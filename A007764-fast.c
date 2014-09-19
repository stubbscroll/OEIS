/* efficient and parallel version of a program that calculates
   terms of the sequence A007764 (and related sequences such as
   A064297). this program is stand-alone. warning, the code is written entirely
   with optimization in mind, don't expect it to be easy to read and
   understand. programmed by Ruben Spaans in september-october 2012.

   this program was used to calculate n=22,23,24 (as well as verifying
   all lower n correctly), it took around a CPU-month on a machine with
   4 eight-core CPUs and 1 TB RAM (huge thanks to Rune Jensen and Q2S at
   the Norwegian University of Science and Technology for letting me run
   the program on their computer) (n=24 required around 700 GB).
   this program only calculates the answer modulo a large prime (just
   below 2^63). enough runs were done so that the product of all
   primes used exceeded a guaranteed upper bound for the answer, then
   CRT was used to obtain the final anawer. (the CRT bit was done with
   an external program). see the bottom of this file to see how to use the
   program.

   a paper was under preparation, but has been on hold for several years
   (as of september 2014). with luck, it might see the face of the earth.
   don't hold your breath. */

/* find the number of simple paths along the edges of an n*n grid between the
   opposite corners.
   equivalently: find the number of simple paths along the cells of an
   (n+1)*(n+1) grid, from cell (0,0) to cell (n,n). cells are 4-connected.

   algorithm: dynamic programming. state representation is the same as in
   M. Bousquet-Melou, A. J. Guttmann and I. Jensen:
   "Self-avoiding walks crossing a square". when scanning cells row by row
   from left to right, we have a frontier between processed and unprocessed
   cells:

   ******
   ****** <- * denotes scanned cells
   **....    . denotes unscanned cells
   ......

   given a grid of n*n internal cells, the frontier consists of (n+1)
   cell boundaries. each cell boundary can contain a crossing edge or none.
   the crossing edge can belong to the path from (0,0), or it can be an
   "incoming" or "outgoing" edge belonging to a segment having two ends
   going through the frontier. assuming we scan the frontier from left to
   right, the first loose end belonging to a segment with loose ends is
   the incoming edge, and the second loose end of the same segment is the
   outgoint edge. let's use this encoding: 0=no edge, 1=loose end from
   upper left, 2=incoming edge, 3=outgoing edge. a base-4 number with (n+1)
   bits can represent all such boundaries. we can then store all partial
   states in a hashmap with this base-4 number as a key, and the number of
   ways to reach this partial state as the value. notice that we don't need to
   store the current (x,y)-position in the state, as this is implicitly known
   via the i,j iteration variables. whenever we process a given cell we only
   have two hashmaps simultaneously in memory: the one we read from and the one
   we write new values to (for the next iteration).

   we can improve the representation by making a mapping from the base-4 number
   to integers 0..(m-1) where m is the number of base-4 numbers that represent
   legal states. then we can store the values in a regular array and ditch the
   hashmap. this is implemented in the current version.

   further speedup: pthreads with 30 threads (running on a machine with 32
   physical cores), 22x speedup achieved.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>

typedef unsigned long long ull;
typedef long long ll;
typedef unsigned int uint;
/* define the number of bits! 64-bits compared to 16-bits is expected
   to be 4 times faster, using 4 times as much memory */
#define BITS64
//#define BITS16

/* parallel variables */
/* please let MUTEX be a power of two */
/* please tune this value to the target machine */
#define MUTEX (1<<21)

int THREAD=8;

static double gettime() {
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+t.tv_usec/1000000.;
}

#define MAX 32
ull dp[2][MAX][MAX];

/* calculate the number of states (and partial states) for rank/unrank */
/* the sequence dp[1][0][n] is actually A002026 */
static void init(uint max) {
	int i,j,k;
	for(k=0;k<2;k++) for(i=0;i<MAX;i++) for(j=0;j<MAX;j++) dp[k][i][j]=0;
	dp[0][0][0]=1;
	for(i=0;i<MAX-1;i++) for(k=0;k<2;k++) for(j=0;j<MAX-1;j++) if(dp[k][j][i]) {
		dp[k][j][i+1]+=dp[k][j][i];
		if(!k && !j) dp[1][j][i+1]+=dp[k][j][i];
		if(j) dp[k][j-1][i+1]+=dp[k][j][i];
		if(j<MAX-1) dp[k][j+1][i+1]+=dp[k][j][i];
	}
}

/* convert from integer rank to representation in linear time */
static ull unrank(int i,ull r) {
	int j=0;
	ull c0,mask=0;
	while(i--) {
		c0=dp[1][j][i];
		if(r<c0) mask<<=2;
		else {
			r-=c0;
			c0=(!j)?dp[0][0][i]:0;
			if(r<c0) { mask=(mask<<2)+1; goto seen1; }
			else {
				r-=c0;
				c0=(j)?dp[1][j-1][i]:0;
				if(r<c0) mask=(mask<<2)+2,j--;
				else r-=c0,mask=(mask<<2)+3,j++;
			}
		}
	}
	return mask;
seen1:
	while(i--) {
		c0=dp[0][j][i];
		if(r<c0) mask<<=2;
		else {
			r-=c0;
			c0=(j)?dp[0][j-1][i]:0;
			if(r<c0) mask=(mask<<2)+2,j--;
			else r-=c0,mask=(mask<<2)+3,j++;
		}
	}
	return mask;
}

/* convert from represention to integer rank in linear time */
static ull rank(int i,ull mask) {
	int j=0,cur;
	ull r=0;
	while(i--) {
		cur=(mask>>(i<<1))&3;
		if(cur==2) {
			r+=dp[1][j][i];
			j--;
		} else if(cur==3) {
			r+=dp[1][j][i]+(j?dp[1][j-1][i]:dp[0][0][i]);
			j++;
		} else if(cur) { r+=dp[1][j][i]; goto seen1; }
	}
	return r;
seen1:
	while(i--) {
		cur=(mask>>(i<<1))&3;
		if(cur==2) {
			r+=dp[0][j][i];
			j--;
		} else if(cur==3) {
			r+=dp[0][j][i];
			if(j) r+=dp[0][j-1][i];
			j++;
		}
	}
	return r;
}

/* we want to solve huge instances, and wish to use as little memory as
   possible for storing intermediate values. therefore, we use short
   ints to store the intermediate values. we use only half the range so that
   addition+modulo becomes easier. */

#ifdef BITS16
	typedef unsigned short u16;
#elif defined(BITS64)
	typedef ull u16;
#endif

pthread_mutex_t mutex[MUTEX];
ull parans[64];
int id[64];
u16 gmod;

u16 *prev,*cur;
#define ADD(ix,c,mod) { \
	pthread_mutex_lock(&mutex[ix&(MUTEX-1)]);\
	cur[ix]+=c; \
	if(cur[ix]>=mod) cur[ix]-=mod; \
	pthread_mutex_unlock(&mutex[ix&(MUTEX-1)]);\
}
static int swap[16]={0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15};
ull pstart[64],pend[64];

int gi,gj,gn,gm,gw;
void *threadcalc(void *threadid) {
	int i=gi,j=gj,n=gn,m=gm,w=gw,k,left,up,look,ix=*((int *)threadid),l;
	ull z,mask,nz,newmask,o;
	u16 mod=gmod,c,r=0;
	if(i<m-1 && j<n-1) {
		/* regular cell */
		for(z=pstart[ix];z<pend[ix];z++) if((c=prev[z])) {
			mask=unrank(w,z);
			left=(mask>>(j<<1))&3;
			up=(mask>>((j<<1)+2))&3;
			if(left==3 && up==2) {
				/* join, easy case: 32 => 00 */
				nz=rank(w,mask&(~(15ULL<<(j<<1))));
				ADD(nz,c,mod);
			} else if(left==2 && up==2) {
				/* join 22: find mate of right 2, change it from 3 to 2 */
				for(k=j+2,l=1;;k++) {
					o=(mask>>(k<<1))&3;
					if(o==2) l++;
					else if(o==3) {
						l--;
						if(!l) break;
					}
				}
				newmask=mask&(~(15ULL<<(j<<1)));
				newmask=newmask^(1ULL<<(k<<1));
				nz=rank(w,newmask);
				ADD(nz,c,mod);
			} else if(left==3 && up==3) {
				/* join 33: find mate of left 3, change it from 2 to 3 */
				for(k=j-1,l=1;;k--) {
					o=(mask>>(k<<1))&3;
					if(o==3) l++;
					else if(o==2) {
						l--;
						if(!l) break;
					}
				}
				newmask=mask&(~(15ULL<<(j<<1)));
				newmask=newmask^(1ULL<<(k<<1));
				nz=rank(w,newmask);
				ADD(nz,c,mod);
			} else if(left==1 && up) {
				/* we have 12, find up's mate, change it from 3 to 1 */
				for(k=j+2,l=1;;k++) {
					o=(mask>>(k<<1))&3;
					if(o==2) l++;
					else if(o==3) {
						l--;
						if(!l) break;
					}
				}
				nz=rank(w,mask&(~(15ULL<<(j<<1)))&(~(2ULL<<(k<<1))));
				ADD(nz,c,mod);
			} else if(left && up==1) {
				/* we have 31, find left's mate, change it from 2 to 1 */
				for(k=j-1,l=1;;k--) {
					o=(mask>>(k<<1))&3;
					if(o==3) l++;
					else if(o==2) {
						l--;
						if(!l) break;
					}
				}
				newmask=mask&(~(15ULL<<(j<<1)));
				newmask=newmask^(3ULL<<(k<<1));
				nz=rank(w,newmask);
				ADD(nz,c,mod);
			} else if(left) { /* extend single edge, case 1 */
				/* extend down: no change in mask */
				ADD(z,c,mod);
				/* extend to the right: x0 becomes 0x, but only if next
					 cell isn't 2 */
				look=(mask>>((j<<1)+4))&3;
				if(left!=2 || look!=3) {
					nz=rank(w,(mask&(~(15ULL<<(j<<1))))|((ull)swap[left]<<(j<<1)));
					ADD(nz,c,mod);
				}
			} else if(up) { /* extend single edge, case 2 */
				/* extend down: 0x becomes x0 */
				nz=rank(w,(mask&(~(15ULL<<(j<<1))))|((ull)swap[up<<2]<<(j<<1)));
				ADD(nz,c,mod);
				/* extend to the right: no change in mask */
				look=(mask>>((j<<1)+4))&3;
				if(up!=2 || look!=3) ADD(z,c,mod);
			} else if(!(up|left)) { /* no edge */
				/* place nothing */
				ADD(z,c,mod);
				/* place 23 */
				nz=rank(w,mask|(14ULL<<(j<<1)));
				ADD(nz,c,mod);
			} else printf("error uncatched regular %d %d\n",left,up);
		}
	} else if(i<m-1 && j==n-1) {
		/* right column: edges cannot go to the right, mask<<2 */
		for(z=pstart[ix];z<pend[ix];z++) if((c=prev[z])) {
			mask=unrank(w,z);
			left=(mask>>(j<<1))&3;
			up=(mask>>((j<<1)+2))&3;
			if(left==3 && up==3) {
				/* join 33: find mate of left 3, change it from 2 to 3 */
				for(k=j-1,l=1;;k--) {
					o=(mask>>(k<<1))&3;
					if(o==3) l++;
					else if(o==2) {
						l--;
						if(!l) break;
					}
				}
				newmask=mask&(~(15ULL<<(j<<1)));
				newmask=newmask^(1ULL<<(k<<1));
				nz=rank(w,newmask<<2);
				ADD(nz,c,mod);
			} else if(left==1 && up) {
				/* we have 12, find up's mate, change it from 3 to 1 */
				for(k=j+2,l=1;;k++) {
					o=(mask>>(k<<1))&3;
					if(o==2) l++;
					else if(o==3) {
						l--;
						if(!l) break;
					}
				}
				nz=rank(w,(mask&(~(15ULL<<(j<<1)))&(~(2ULL<<(k<<1))))<<2);
				ADD(nz,c,mod);
			} else if(left && up==1) {
				/* we have 31, find left's mate, change it from 2 to 1 */
				for(k=j-1,l=1;;k--) {
					o=(mask>>(k<<1))&3;
					if(o==3) l++;
					else if(o==2) {
						l--;
						if(!l) break;
					}
				}
				newmask=mask&(~(15ULL<<(j<<1)));
				newmask=newmask^(3ULL<<(k<<1));
				nz=rank(w,newmask<<2);
				ADD(nz,c,mod);
			} else if(left) { /* extend single edge, case 1 */
				/* extend down */
				newmask=mask<<2;
				nz=rank(w,newmask);
				ADD(nz,c,mod);
			} else if(up) { /* extend single edge, case 2 */
				/* extend down: 0x becomes x0 */
				nz=rank(w,((mask&(~(15ULL<<(j<<1))))|((ull)swap[up<<2]<<(j<<1)))<<2);
				ADD(nz,c,mod);
			} else if(!(up|left)) { /* no edge */
				/* place nothing */
				newmask=mask<<2;
				nz=rank(w,newmask);
				ADD(nz,c,mod);
			} else printf("error uncatched regular %d %d\n",left,up);
		}
	} else if(i==m-1 && j<n-1) {
		/* lower row_ edges cannot go down */
		for(z=pstart[ix];z<pend[ix];z++) if((c=prev[z])) {
			mask=unrank(w,z);
			left=(mask>>(j<<1))&3;
			up=(mask>>((j<<1)+2))&3;
			if(left==3 && up==2) {
				/* join, easy case: 32 => 00 */
				nz=rank(w,mask&(~(15ULL<<(j<<1))));
				ADD(nz,c,mod);
			} else if(left==2 && up==2) {
				/* join 22: find mate of right 2, change it from 3 to 2 */
				for(k=j+2,l=1;;k++) {
					o=(mask>>(k<<1))&3;
					if(o==2) l++;
					else if(o==3) {
						l--;
						if(!l) break;
					}
				}
				newmask=mask&(~(15ULL<<(j<<1)));
				newmask=newmask^(1ULL<<(k<<1));
				nz=rank(w,newmask);
				ADD(nz,c,mod);
			} else if(left==3 && up==3) {
				/* join 33: find mate of left 3, change it from 2 to 3 */
				for(k=j-1,l=1;;k--) {
					o=(mask>>(k<<1))&3;
					if(o==3) l++;
					else if(o==2) {
						l--;
						if(!l) break;
					}
				}
				newmask=mask&(~(15ULL<<(j<<1)));
				newmask=newmask^(1ULL<<(k<<1));
				nz=rank(w,newmask);
				ADD(nz,c,mod);
			} else if(left==1 && up) {
				/* we have 12, find up's mate, change it from 3 to 1 */
				for(k=j+2,l=1;;k++) {
					o=(mask>>(k<<1))&3;
					if(o==2) l++;
					else if(o==3) {
						l--;
						if(!l) break;
					}
				}
				nz=rank(w,mask&(~(15ULL<<(j<<1)))&(~(2ULL<<(k<<1))));
				ADD(nz,c,mod);
			} else if(left && up==1) {
				/* we have 31, find left's mate, change it from 2 to 1 */
				for(k=j-1,l=1;;k--) {
					o=(mask>>(k<<1))&3;
					if(o==3) l++;
					else if(o==2) {
						l--;
						if(!l) break;
					}
				}
				newmask=mask&(~(15ULL<<(j<<1)));
				newmask=newmask^(3ULL<<(k<<1));
				nz=rank(w,newmask);
				ADD(nz,c,mod);
			} else if(left) { /* extend single edge, case 1 */
				/* extend to the right: x0 becomes 0x, but only if next
					 cell isn't 2 */
				look=(mask>>((j<<1)+4))&3;
				if(left!=2 || look!=3) {
					nz=rank(w,(mask&(~(15ULL<<(j<<1))))|((ull)swap[left]<<(j<<1)));
					ADD(nz,c,mod);
				}
			} else if(up) { /* extend single edge, case 2 */
				/* extend to the right: no change in mask */
				look=(mask>>((j<<1)+4))&3;
				if(up!=2 || look!=3) ADD(z,c,mod);
			} else if(!(up|left)) { /* no edge */
				/* place nothing */
				ADD(z,c,mod);
			} else printf("error uncatched regular %d %d\n",left,up);
		}
	} else if(i==m-1 && j==n-1) {
		/* lower right corner, just take the sum */
		for(z=pstart[ix];z<pend[ix];z++) if((c=prev[z])) {
			r+=c;
			if(r>=mod) r-=mod;
		}
		parans[ix]=r;
	}
	pthread_exit(NULL);
	return NULL;
}

void *threaddelcur(void *threadid) {
	int ix=*((int *)threadid);
	memset(cur+pstart[ix],0,sizeof(u16)*(pend[ix]-pstart[ix]));
	pthread_exit(NULL);
	return NULL;
}

pthread_t threads[64];

static u16 calc(int n,int m,u16 mod) {
	ull num=dp[1][0][n+1],inc;
	int i,j,w=n+1,rc,k;
	u16 r=0,*t;
	if(!(prev=malloc(sizeof(u16)*num))) { puts("out of memory"); exit(1); }
	if(!(cur=malloc(sizeof(u16)*num))) { puts("out of memory"); exit(1); }
	memset(prev,0,sizeof(u16)*num);
	memset(cur,0,sizeof(u16)*num);
	/* force starting edge to go down. the exact number of paths is twice the
	   answer we get. */
	prev[rank(w,1)]=1;
	if(m!=n) prev[rank(w,4)]=1;	/* cannot use symmetry if non-square */
	gn=n; gm=m; gw=w; gmod=mod;
	for(k=0;k<THREAD;k++) id[k]=k;
	inc=num/THREAD;
	for(k=0;k<THREAD;k++) {
		pstart[k]=k*inc;
		pend[k]=(k+1)*inc;
	}
	pend[THREAD-1]=num;
	/* we are not worried about creating 2 million mutexes and creating/tearing
	   down threads n^2 times, as the time for doing this is insignificant
	   compared to the calculation for large n */
	for(k=0;k<MUTEX;k++) pthread_mutex_init(&mutex[k],0);
	for(i=0;i<m;i++) for(j=0;j<n;j++) {
		fprintf(stderr,"[%d %d] ",i,j);fflush(stderr);
		if(i==0 && j==0) continue;
		gi=i; gj=j;
		for(k=0;k<THREAD;k++) {
			rc=pthread_create(&threads[k],NULL,threadcalc,id+k);
			if(rc) puts("error creating threads");
		}
		for(k=0;k<THREAD;k++) pthread_join(threads[k],NULL);
		t=prev; prev=cur; cur=t;
		if(i!=m-1 || j!=n-1) {
			/* delete array with threads */
			for(k=0;k<THREAD;k++) {
				rc=pthread_create(&threads[k],NULL,threaddelcur,id+k);
				if(rc) puts("error creating threads");
			}
			for(k=0;k<THREAD;k++) pthread_join(threads[k],NULL);
		}
	}
	fprintf(stderr,"\n");fflush(stderr);
	free(cur);
	free(prev);
	for(i=0;i<THREAD;i++) r=(r+parans[i])%mod;
	if(n==m) r=(r*2)%mod;
	return r;
}
#undef ADD

/* m,n: dimensions of rectangle with n*m internal cells */
/* when solving A007764, invoke with n+1,n+1 */
static void solve(int n,int m,ull prime) {
	int t;
	double start=gettime();
	if(n>m) t=n,n=m,m=t;
	/* require n<=m, frontier width is n+1 */
	printf("answer mod %llu for %d x %d: %llu\n",prime,n-1,m-1,calc(n,m,prime));
	start=gettime()-start;
	printf("time elapsed %.3f\n",start);
}

ull primes[15]={
	9223372036854775783,9223372036854775643,9223372036854775549,
	9223372036854775507,9223372036854775433,9223372036854775421,
	9223372036854775417,9223372036854775399,9223372036854775351,
	9223372036854775337,9223372036854775291,9223372036854775279
};

int main(int argc,char **argv) {
	int n,m,ix;
	if(argc<5) {
		puts("usage: sq64 n m ix thr");
		puts("where n*m is the problem size and ix is an index into a list of prime numbers");
		puts("(between 0 and 11 inclusive) and thr is the number of threads");
		return 0;
	}
	init(29);
	sscanf(argv[1],"%d",&n);
	sscanf(argv[2],"%d",&m);
	sscanf(argv[3],"%d",&ix);
	sscanf(argv[4],"%d",&THREAD);
	solve(n+1,m+1,primes[ix]);
	return 0;
}
