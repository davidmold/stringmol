/* Copyright (C) 2009-2016 Simon Hickinbotham                           */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL										*/

/* STRINGMOL is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

//#define USE_SDL



#include <stdlib.h>
#include <stdio.h>

//Used for wait time...
#include <time.h>

#ifdef USE_SDL
#ifdef _MSC_VER
#include <SDL.h>
#else
#include <SDL/SDL.h>
#endif /* _MSC_VER */
#endif


#include "randutil.h"
#include "params.h"
#include "hsort.h"
//}

// Stringmol
#include "alignment.h"

// Metabolism
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"

// Writing PNGs
#include "lodepng.h"
#include <iostream>
#include "setupSM.h"

typedef enum td_pic{
	pic_spp,
	pic_len,
	pic_bound
}pictype;










/* ---- slim output helpers ---- */
#include <vector>
#include <stdint.h>
#include <unistd.h>

/* species.tsv: spp <tab> born_step <tab> active_parent <tab> passive_parent <tab> sequence.
 * SMspp prepends, so new species are at the head; write those above the last one logged. */
static void log_new_species(SMspp *SP, FILE *fp, int *lastspp){
	std::vector<l_spp *> fresh;
	for(l_spp *ps=SP->species; ps!=NULL && ps->spp > *lastspp; ps=ps->next) fresh.push_back(ps);
	for(size_t k=fresh.size(); k-->0; ){
		l_spp *ps = fresh[k];
		int pa=-1, pp=-1;
		if(ps->pp){ if(ps->pp->pa) pa=ps->pp->pa->spp; if(ps->pp->pp) pp=ps->pp->pp->spp; }
		fprintf(fp,"%d\t%d\t%d\t%d\t%s\n",ps->spp,ps->tspp,pa,pp,ps->S);
		if(ps->spp > *lastspp) *lastspp = ps->spp;
	}
	fflush(fp);
}

/* grid.snap: repeated frames of [uint32 step][uint32 gx][uint32 gy][uint32 nbytes][zlib data],
 * data = gx*gy uint32 species ids (row-major, y*gx+x), 0 = empty; bound molecules get bit 31. */
static void write_grid_snapshot(FILE *fp, smsprun *run, unsigned int step){
	const int gx=run->gridx, gy=run->gridy;
	std::vector<unsigned char> raw(gx*gy*4,0), z;
	for(int y=0;y<gy;y++) for(int x=0;x<gx;x++){
		s_ag *a = run->grid[x][y]; uint32_t v=0;
		if(a && a->spp){ v = (uint32_t) a->spp->spp; if(a->status!=B_UNBOUND) v |= 0x80000000u; }
		memcpy(&raw[(y*gx+x)*4],&v,4);
	}
	lodepng::compress(z,raw);
	uint32_t h[4] = { step, (uint32_t)gx, (uint32_t)gy, (uint32_t)z.size() };
	fwrite(h,4,4,fp); fwrite(z.data(),1,z.size(),fp); fflush(fp);
}

/* checkpoint = upstream full-state dump (out1_N.conf + RNGstate_N.txt), keeping only the newest few */
static void write_checkpoint(stringPM *A, SMspp *SP, unsigned int step, unsigned int keep, unsigned int every){
	char fn[128];
	sprintf(fn,"out1_%05u.conf",step);
	FILE *fp=fopen(fn,"w"); A->print_conf(fp); fclose(fp);
	if(step >= keep*every){
		unsigned int old = step - keep*every;
		if(old % every == 0 && old>0){
			sprintf(fn,"out1_%05u.conf",old); unlink(fn);
			sprintf(fn,"RNGstate_%u.txt",old); unlink(fn);
		}
	}
}

int main(int argc, char *argv[]) {
//	int smspatial(int argc, char *argv[]) {

	printf("Hello spatial stringmol world\n");

	//int rtype = atoi(argv[1]);

	SMspp		SP;
	stringPM	A(&SP);

	smsprun *run;
	run = NULL;

	/* FIX: upstream never seeded the RNG in this entry point, so every run used the
	 * Mersenne Twister default seed (4357) whatever RANDSEED said. Honour RANDSEED.
	 * Runs made before this fix are reproduced exactly with RANDSEED 4357. */
	A.randseed = init_randseed(argv[2],0);
	smspatial_init(argv[2],&A,&run,1);

	int bt=0,ct=0;
	ct = A.nagents(A.nowhead,-1);
	printf("Initialisation done, number of molecules is %d\n",ct);

	/* Set up SDL if we're using it */
#ifdef USE_SDL

	int **grid;
	grid = (int **)malloc(run->gridx*sizeof(int *));
	for(int i=0;i<run->gridx;i++){
		grid[i] = (int *) malloc(run->gridy*sizeof(int));
		memset(grid[i],0,run->gridy*sizeof(int));
	}

	printf("grid allocated\n");fflush(stdout);

	SDL_Surface *screen;
	SDL_Event sdlEvent;
	if (SDL_Init(SDL_INIT_VIDEO) < 0 ) {
		fprintf(stderr,"*** Unable to init SDL: %s ***\n",SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);
	SDL_WM_SetCaption("Spatial Stringmol","nanopond");
	screen = SDL_SetVideoMode(run->gridy,run->gridx,8,SDL_SWSURFACE);
	if (!screen) {
		fprintf(stderr, "*** Unable to create SDL window: %s ***\n", SDL_GetError());
		exit(1);
	}
	const uintptr_t sdlPitch = screen->pitch;
	printf("SDL window allocated, pitch is %d \n",screen->pitch);fflush(stdout);
#endif /* USE_SDL */


//	while(A.nagents(A.nowhead,-1)){
	/* Output mode. Default is slim: popdy every 100 steps, species log at birth,
	 * compressed grid snapshots every SNAPEVERY steps, full checkpoints every CKPTEVERY
	 * steps (keeping KEEPCKPT). LEGACYOUT 1 restores the upstream every-100-steps dumps. */
	unsigned int legacy_out=0, snapevery=500, ckptevery=50000, keepckpt=3;
	readordef_param_int(argv[2],"LEGACYOUT",&legacy_out,0,0);
	readordef_param_int(argv[2],"SNAPEVERY",&snapevery,500,0);
	readordef_param_int(argv[2],"CKPTEVERY",&ckptevery,50000,0);
	readordef_param_int(argv[2],"KEEPCKPT",&keepckpt,3,0);
	if(!snapevery) snapevery=500;
	if(!ckptevery) ckptevery=50000;
	FILE *splog=NULL, *gridf=NULL; int lastspp=-1;
	if(!legacy_out){
		splog = fopen("species.tsv","a");
		gridf = fopen("grid.snap","ab");
		log_new_species(&SP,splog,&lastspp);
		write_grid_snapshot(gridf,run,A.extit);
	}

	while(A.extit < 1000000){// && A.nagents(A.nowhead,-1)>0){

		smspatial_step(&A,run);
		ct = A.nagents(A.nowhead,-1);
		bt = ct - A.nagents(A.nowhead,B_UNBOUND);
#ifdef DODEBUG
		printf("Nowhead is %d, Nexthead is %d\n",A.nowhead,A.nexthead);
		s_ag *p;
		p=A.nowhead;
		int mno=0;
		while(p!=NULL){
			int x,y;
			find_ag_gridpos(p,run,&x,&y);
			printf("%d, %d, [%d,%d]  status: %d, bound to %d / %d, prev = %d, next = %d\n",++mno,p,x,y,p->status,p->exec,p->pass,p->prev,p->next);
			p = p->next;
		}
#endif

//		if(A.extit == 66396)
//			printf("Pauuuuse\n!");
		A.extit++;
		//if(!(A.extit%100))
		//		printf("Step %d done, number of molecules is %d, nbound = %d\n",(int) A.extit,ct,bt);
		if(!(A.extit%100)){
			printf("Step %d done, number of molecules is %d, nbound = %d\n",(int) A.extit,ct,bt);
			fflush(stdout);
			if(legacy_out){
				/* upstream output: full species list + full state + RNG every 100 steps */
				FILE *fp;char fn[128];
				sprintf(fn,"splist%d.dat",A.extit);
				fp = fopen(fn,"w");
				SP.print_spp_list(fp);
				fclose(fp);
				printsppct(&A,A.extit);
				sprintf(fn,"out1_%05d.conf",A.extit);
				fp = fopen(fn,"w");
				A.print_conf(fp);
				fclose(fp);
			}
			else{
				printsppct(&A,A.extit);              // population counts (popdy001.dat)
				log_new_species(&SP,splog,&lastspp);  // each species once, at birth (species.tsv)
			}
		}
		if(!legacy_out){
			if(!(A.extit%snapevery)) write_grid_snapshot(gridf,run,A.extit);
			if(!(A.extit%ckptevery)) write_checkpoint(&A,&SP,A.extit,keepckpt,ckptevery);
		}

		//if(!(A.extit%10)){
		//	printf("completed iteration %u of simulation\n",A.extit);fflush(stdout);
		//	//printf("grid dimenstions are %d,%d\n",run->gridx,run->gridy);fflush(stdout);
		//}

		//Create the PNG
		std::vector<unsigned char> image;
		image.resize(run->gridx * run->gridy * 4);

#ifdef USE_SDL
    pictype tp = pic_len;
		int x,y,val;
		for(x=0;x<run->gridx;++x){
			for (y=0;y<run->gridy;++y) {
				switch(run->status[x][y]){
				case G_EMPTY:
					val=0;
					break;
				default:
					switch(tp){
					case pic_spp:
						val=run->grid[x][y]->spp->spp * 50;
						break;
					case pic_len:
						val=strlen(run->grid[x][y]->spp->S) * 50;
						break;
					}
				}

				((uint8_t *)screen->pixels)[y + (x * sdlPitch)] = val;//getColor(grid[x][y]);

				if(!(A.extit%10)){
					/*Get the rgb values for writing to PNG*/
					uint8_t r ;
					uint8_t g ;
					uint8_t b ;
					SDL_GetRGB(((uint8_t *)screen->pixels)[y + (x * sdlPitch)], screen->format ,  &r, &g, &b );
					//printf("R is %d G is %d B is %d\n",r,g,b);
					image[4 * run->gridx * y + 4 * x + 0] = r;//255 * !(x & y);
					image[4 * run->gridx * y + 4 * x + 1] = g;//x ^ y;
					image[4 * run->gridx * y + 4 * x + 2] = b;//x | y;
					image[4 * run->gridx * y + 4 * x + 3] = 255;
				}
			}
		}
#endif


#ifdef USE_SDL
		if(!(A.extit%10)){
			char filename[128];
			sprintf(filename,"frame%07d.png",A.extit);
			encodeOneStep(filename, image, run->gridx, run->gridy);
		}

		while (SDL_PollEvent(&sdlEvent)) {
			if (sdlEvent.type == SDL_QUIT) {
				fprintf(stderr,"[QUIT] Quit signal received!\n");

				//TODO: Write the final splist and popdy files here!

				exit(0);
			}
		}

		//printf("sdl grid updated\n",A.extit);fflush(stdout);

        SDL_UpdateRect(screen,0,0,0,0);//run->gridx,run->gridy);

		//printf("sdl rect updated\n",A.extit);fflush(stdout);
#endif

		/* FIX: stop once the population is extinct, after the next 100-step dump so the
		 * final (empty) state is on disk. Upstream ran on to step 1M on an empty grid. */
		if(ct==0 && !(A.extit%100)){
			printf("EXTINCT at step %d\n",(int) A.extit);
			break;
		}
	}

	if(!legacy_out){
		/* final state: full species list (with parent-pair counts) and a last checkpoint */
		FILE *fp = fopen("splist_final.dat","w"); SP.print_spp_list(fp); fclose(fp);
		write_checkpoint(&A,&SP,A.extit,keepckpt+1,ckptevery);
		fclose(splog); fclose(gridf);
	}
	printf("FINISHED smspatial\n");
//	fflush(stdout);
	return 0;
}


