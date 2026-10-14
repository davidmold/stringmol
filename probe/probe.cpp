// Deterministic binding probe: P(bind) and copier/copied roles for molecule pairs,
// using the same get_sw/get_bprob/set_exec logic as the simulation.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "mt19937-2.h"
#include "randutil.h"
#include "params.h"
#include "alignment.h"
#include "instructions.h"
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"
#include "lodepng.h"

static s_ag *mk(stringPM &A, const char *s){
	s_ag *a = A.make_ag('A');
	a->S = (char*)calloc(A.maxl0,1); strcpy(a->S,s); a->len = strlen(s);
	return a;
}
int main(int argc, char **argv){
	SMspp SP; stringPM A(&SP);
	A.load(argv[1],NULL,0,0);
	std::vector<std::pair<std::string,std::string>> R, P;
	FILE *fp = fopen(argv[2],"r"); char kind[8], name[64], seq[4096];
	while(fscanf(fp,"%7s %63s %4095s",kind,name,seq)==3){
		(kind[0]=='R'?R:P).push_back({name,seq});
	}
	printf("%-10s %-10s %8s %8s %8s  %s\n","replicase","partner","P(r->p)","P(p->r)","mean","copier when bound (r-first / p-first)");
	for(auto &r:R){
		std::vector<std::pair<std::string,std::string>> all = R; all.insert(all.end(),P.begin(),P.end());
		for(auto &p:all){
			s_ag *a=mk(A,r.second.c_str()), *b=mk(A,p.second.c_str());
			align sw1, sw2;
			float p1 = A.get_sw(a,b,&sw1);
			float p2 = A.get_sw(b,a,&sw2);
			const char *c1 = (sw1.s1>=sw1.s2)?"R":"P";   // per set_exec: larger start index is active
			const char *c2 = (sw2.s1>=sw2.s2)?"P":"R";
			printf("%-10s %-10s %8.3f %8.3f %8.3f  %s / %s\n",r.first.c_str(),p.first.c_str(),p1,p2,(p1+p2)/2,c1,c2);
		}
		printf("\n");
	}
	return 0;
}
