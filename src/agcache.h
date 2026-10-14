/* PERF: ordered mirror of a stringmol agent list (nowhead / nexthead).
 * Keeps list order exactly, with O(log n) append, removal and k-th lookup:
 * slots are append-only, removed agents leave a dead slot, and a Fenwick tree over
 * the live flags turns "k-th live agent" into a tree descent. Dead slots are compacted
 * away once they outnumber live ones. The linked list stays the source of truth. */
#ifndef AGCACHE_H_
#define AGCACHE_H_
#include <vector>
#include <cstddef>
#include "SMspp.h"

struct ag_cache {
	std::vector<s_ag *> slot;
	std::vector<char> live;
	std::vector<int> bit;      // 1-based Fenwick tree over live[], bit[0] unused
	int nlive = 0;
	s_ag *tail = NULL;
	bool valid = false;

	void clear(){ slot.clear(); live.clear(); bit.assign(1,0); nlive=0; tail=NULL; valid=true; }

	int prefix(int i) const { int s=0; for(; i>0; i -= i & -i) s += bit[i]; return s; }

	void push(s_ag *a){
		if(bit.empty()) bit.assign(1,0);
		slot.push_back(a); live.push_back(1);
		int i = (int) slot.size();                 // 1-based index of the new slot
		bit.push_back(1 + prefix(i-1) - prefix(i - (i & -i)));
		a->cslot = i - 1;
		nlive++; tail = a;
	}

	bool holds(s_ag *a) const {
		int k = a->cslot;
		return k >= 0 && k < (int) slot.size() && slot[k] == a && live[k];
	}

	/* prev = a->prev as it was before a was unlinked from its list */
	void remove(s_ag *a, s_ag *prev){
		int k = a->cslot;
		live[k] = 0; nlive--;
		for(int i = k+1; i < (int) bit.size(); i += i & -i) bit[i]--;
		if(a == tail) tail = prev;
		if((int) slot.size() > 2*nlive + 4096) compact();
	}

	s_ag *kth(int k) const {                       // 0-based k-th live agent
		int pos = 0, n = (int) bit.size() - 1, step = 1;
		while(step*2 <= n) step *= 2;
		for(; step; step >>= 1)
			if(pos + step <= n && bit[pos+step] <= k){ pos += step; k -= bit[pos]; }
		return slot[pos];                          // pos is 0-based slot of the answer
	}

	s_ag *first() const { return nlive ? kth(0) : NULL; }

	void compact(){
		std::vector<s_ag *> keep;
		keep.reserve(nlive);
		for(size_t i=0;i<slot.size();i++) if(live[i]) keep.push_back(slot[i]);
		s_ag *t = tail;
		clear();
		for(s_ag *a: keep) push(a);
		tail = t;
	}

	void rebuild(s_ag *head){
		clear();
		for(s_ag *p=head; p!=NULL; p=p->next) push(p);
	}
};
#endif
