// Reaction probe: force-bind two molecules exactly as the simulation would (get_sw + set_exec),
// then execute the reaction to completion with exec_step, mutation off. Reports who copied whom,
// how many steps it took, and every molecule produced.
//
//   react CONF PAIRS [maxsteps]
// PAIRS lines: "R name seq" or "P name seq"; every listed molecule is reacted against every other
// (and itself), with each one taking a turn as the initiator.
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

struct mol { std::string name, seq; };
static std::vector<mol> M;

static s_ag *mk(stringPM &A, const std::string &s){
	s_ag *a = A.make_ag('A');
	a->S = (char*)calloc(A.maxl0,1); strncpy(a->S,s.c_str(),A.maxl0-1); a->len = strlen(a->S);
	return a;
}
static std::string ident(const char *s){
	for(auto &m:M) if(m.seq==s) return m.name;
	return std::string("new(") + std::to_string(strlen(s)) + ")";
}
static void clear_next(stringPM &A){ A.nexthead=NULL; A.cnext.v.clear(); A.cnext.valid=true; }

static void react(stringPM &A, const mol &x, const mol &y, int maxsteps){
	s_ag *a=mk(A,x.seq), *b=mk(A,y.seq);
	align sw;
	float p = A.get_sw(a,b,&sw);
	A.set_exec(a,b,&sw);
	s_ag *act = (a->status==B_ACTIVE)?a:b, *pass = (act==a)?b:a;
	std::string an = (act==a)?x.name:y.name, pn = (act==a)?y.name:x.name;
	std::vector<std::string> products;
	int steps=0; const char *end="ran out of steps";
	while(steps<maxsteps){
		clear_next(A);
		act->next=act->prev=pass->next=pass->prev=NULL;  // as if extracted from nowhead, like the sim
		A.exec_step(act,pass);
		steps++;
		bool act_in=false, pass_in=false;
		for(s_ag *q=A.nexthead;q;q=q->next){
			if(q==act) act_in=true;
			else if(q==pass) pass_in=true;
			else products.push_back(q->S);
		}
		if(!act_in || !pass_in){ end = !act_in ? "active destroyed" : "passive destroyed"; break; }
		if(act->status!=B_ACTIVE){ end="unbound normally"; break; }
	}
	printf("%-8s initiates with %-8s P(bind)=%.3f  copier=%-8s copied=%-8s steps=%-5d %s\n",
		x.name.c_str(), y.name.c_str(), p, an.c_str(), pn.c_str(), steps, end);
	for(auto &s:products) printf("      -> product %-10s %s\n", ident(s.c_str()).c_str(), s.c_str());
	if(products.empty()) printf("      -> no product\n");
}

int main(int argc, char **argv){
	if(argc<3){ fprintf(stderr,"usage: react CONF PAIRS [maxsteps]\n"); return 1; }
	SMspp SP; stringPM A(&SP);
	A.load(argv[1],NULL,0,0);
	A.subrate = 0; A.indelrate = 0;   // no mutation: show the reaction as written
	int maxsteps = argc>3 ? atoi(argv[3]) : 20000;
	FILE *fp=fopen(argv[2],"r"); if(!fp){ perror(argv[2]); return 1; }
	char kind[8], name[64], seq[4096];
	while(fscanf(fp,"%7s %63s %4095s",kind,name,seq)==3) M.push_back({name,seq});
	setvbuf(stdout,NULL,_IONBF,0);
	for(size_t i=0;i<M.size();i++)
		for(size_t j=0;j<M.size();j++) react(A,M[i],M[j],maxsteps);
	return 0;
}
