/* A109224 minimum number of moves to solve second panex puzzle of order n */
/* rather efficient implementation of dumb bfs. state is as follows. for order
   n, the number of possible blocks per position are given as:
   2*n+1 - 2*n+1
   2*n+1 2*n+1 2*n+1
   2*n-1 2*n-1 2*n-1
   ...
   11 11 11
   9 9 9
   7 7 7
   5 5 5
   3 3 3
   in each cell, the following ids are used:
   0 is space, 1-2 are bottommost blocks, 3-4 are second row of blocks from
   bottom, ..., 2n-1 and 2n are topmost blocks. cannot store blocks on the
   minus.
   store each position packed as a multiradix number (need simple bignum
   arithmetic for that, as we quickly exceed 64 bits for the interesting
   orders). for n=10, there are 21 different things on the board.
   bfs tries to do mostly linear memory access: read from previous iteration,
   try all moves and store them as current iteration. don't do any duplication
   checks until we're done with this iteration, or memory is full (which
   triggers "premature" check before proceeding with iteration). then sort
   current iteration, and do linear scan and remove all positions from current
   iteration that exist in thw two previous iterations, as well as duplicates
   within itself.
   throw enough ram and time at the program and see how far we can solve.
   some results:
   n=4: 0.109 s, 0 mb, visited 57 387 states
   n=5: 1.794 s, 0 mb, visited 1 088 951 states
   n=6: 39.811 s, 2 mb, visited 20 674 713 states
   n=7: 1092 s, 20 mb, visited 389 941 755 states
	 n=8: 24105 s, 152 mb, visited 7 343 745 359 states
   TODO possible improvements:
   - there is only one win state, so bidirectional search will work (extreme
     savings are expected if the branching factor is high, which it probably
     isn't)
   - check if some states never occur. if yes, see if we can make a tighter
     encoding. if we can make a tight encoding, we can achieve massive memory
     savings by doing vbyte encoding on top. might get massive improvement
     just with regular permutation rank and vbyte. disadvantage: the bfs
     routine with delayed duplicate checking and all that becomes a nightmare
   - there are two symmetrical ways of solving. throw away one of them to halve
     memory usage
   - peek in snapdragon's solver which currently trounces this program
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define MEM 12000000000LL
#define MAXY 20
#define MAXST 64

int n;                /* number of blocks per tower */
int packedsize;       /* number of 32-bit uints needed for state */

uint32_t *b;          /* memory area */
long long blen;       /* number of elements in area */
long long bblen;      /* number of states that fit in area */

/* bfs variables */
/* we read from iteration -1 (prev1) and push to current */

long long prev2;      /* start of iteration -2 */
long long prev1;      /* end of iteration -2, start of iteration -1 */
long long cur;        /* end of iteration -1, start of current iteration */
long long cure;       /* end of current iteration (processed) */
long long curp;       /* end of current iteration (unprocessed) */
char wewon;

/* state variables */

char state[3][MAXY];

double gettime() {
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+t.tv_usec/1000000.;
}

/* debug routines */
void printpackedstate(long long at) {
	int i;
	at*=packedsize;
	printf("[");
	for(i=0;i<packedsize;i++) printf("%08x ",b[at+i]);
	puts("]");
}
void printstate() {
	int i,j;
	for(j=0;j<=n;j++) {
		for(i=0;i<3;i++) {
			if(!state[i][j]) putchar(' ');
			else if(state[i][j]&1) putchar(state[i][j]/2+'a');
			else putchar(state[i][j]/2+'A'-1);
		}
		putchar('\n');
	}
}
void printpackedrange(long long a,long long b) {
	while(a<b) printf("%I64d: ",a),printpackedstate(a++);
}

void findstatesize() {
	double b=0;
	int i;
	for(i=1;i<=n;i++) b+=log(i*2+1)*(3+(i==n)*2);
	packedsize=((int)(b/log(2)))/32+1;
}

/* big arithmetic! big numbers are little endian */

/* a=b+c where c is small (a and b can be the same) */
void addst(uint32_t *a,uint32_t *b,uint32_t c) {
	int i;
	for(i=0;i<packedsize;i++) a[i]=b[i];
	a[0]+=c;
	if(a[0]<c) for(i=1;i<n;i++) {
		a[i]++;
		if(a[i]) break;
	}
}

/* helper function for mul */
void addat(uint32_t *a,uint32_t b,int len) {
	if(*a+b<b) {
		*a+=b; a++; len--;
		if(!len) return;
		(*a)++;
		while(!(*a)) {
			a++; len--;
			if(!len) break;
			(*a)++;
		}
	} else *a+=b;
}

/* a=b*c where c is small (a and b can be the same) */
void mulst(uint32_t *a,uint32_t *b,uint32_t c) {
	uint32_t t[MAXST];
	uint64_t r;
	int i;
	for(i=0;i<packedsize;i++) t[i]=0;
	for(i=0;i<packedsize;i++) {
		r=(uint64_t)b[i]*c;
		addat(t+i,r&0xFFFFFFFFu,packedsize-i);
		if(r>=(1ULL<<32) && i<packedsize-1) addat(t+i+1,r>>32,packedsize-i-1);
	}
	for(i=0;i<packedsize;i++) a[i]=t[i];
}

/* return a%b where b is small */
uint32_t modst(uint32_t *a,uint32_t b) {
	uint32_t r=0;
	int i;
	for(i=packedsize-1;i>=0;i--) r=(((uint64_t)r<<32)+a[i])%b;
	return r;
}

/* a=b/c where c is small */
void divst(uint32_t *a,uint32_t *b,uint32_t c) {
	uint64_t carry=0;
	uint32_t t[MAXST];
	int i;
	for(i=packedsize-1;i>=0;i--) {
		t[i]=(b[i]+carry)/c;
		carry=((b[i]+carry)%c)<<32;
	}
	for(i=0;i<packedsize;i++) a[i]=t[i];
}

int multable[MAXY];

/* pack state and store it starting at b[pos*packedsize] */
void pack(int pos) {
	uint32_t a[MAXST];
	int i,j;
	for(i=0;i<packedsize;i++) a[i]=0;
	for(j=0;j<=n;j++) for(i=0;i<3;i++) {
		if(i==1 && !j) continue;
		if(i+j) mulst(a,a,multable[j]);
		addst(a,a,state[i][j]);
	}
	for(i=0;i<packedsize;i++) b[pos*packedsize+i]=a[i];
}

/* unpack from b[pos*packedsize] to state */
void unpack(int pos) {
	uint32_t a[MAXST];
	int i,j;
	for(i=0;i<packedsize;i++) a[i]=b[pos*packedsize+i];
	for(j=n;j>=0;j--) for(i=2;i>=0;i--) {
		if(i==1 && !j) continue;
		state[i][j]=modst(a,multable[j]);
		if(i+j<1) break;
		divst(a,a,multable[j]);
	}
}

void init() {
	int i,j;
	findstatesize();
	/* memory stuff */
	blen=MEM/(sizeof(uint32_t));
	if(!(b=malloc(blen))) printf("error"),exit(1);
	bblen=blen/packedsize;
	prev1=prev2=0; cur=1;
	/* stuff */
	multable[0]=2*n+1;
	for(i=0;i<n;i++) multable[i+1]=(n-i)*2+1;
	wewon=0;
	/* put start state into queue */
	for(i=0;i<3;i++) for(j=0;j<=n;j++) state[i][j]=0;
	for(i=0;i<n;i++) state[0][i+1]=(n-i-1)*2+1,state[2][i+1]=(n-i-1)*2+2;
	pack(prev1);
}

/* sort by packed value (most significant first), qsort interface */
/* potential hotspot */
int compposq(const void *A,const void *B) {
	const uint32_t *a=A,*b=B;
	int i;
	for(i=packedsize-1;i>=0;i--) {
		if(a[i]<b[i]) return -1;
		if(a[i]>b[i]) return 1;
	}
	return 0;
}

/* sort by packed value (most significant first), normal interface */
int comppos(long long a,long long c) {
	return compposq(b+a*packedsize,b+c*packedsize);
}

/* a <- c */
void copypos(long long a,long long c) {
	int i;
	a*=packedsize; c*=packedsize;
	for(i=0;i<packedsize;i++) b[a++]=b[c++];
}

long long peak; /* maximal memory usage */

/* sort and compress all positions from cure to curp, return new curp */
long long sortandcompress(long long cure,long long curp) {
	long long i,j;
	qsort(b+cure*packedsize,curp-cure,packedsize*sizeof(uint32_t),compposq);
	for(i=j=cure+1;i<curp;i++) if(comppos(i-1,i)) {
		if(i!=j) copypos(j,i);
		j++;
	}
	return j;
}

/* remove duplicates from cs,ce against as,ae and bs,be
   return new end pos for ce */
long long removeduplicates2(long long as,long long ae,long long bs,long long be,long long cs,long long ce) {
	long long cto=cs;
	while(cs<ce) {
		while(as<ae && comppos(as,cs)<0) as++;
		while(bs<be && comppos(bs,cs)<0) bs++;
		if(as<ae && !comppos(as,cs)) goto skip;
		if(bs<be && !comppos(bs,cs)) goto skip;
		if(cs!=cto) copypos(cto,cs);
		cto++;
	skip:
		cs++;
	}
	return cto;
}

/* copy block from froms to frome, to to. required. to <= froms <= frome */
void copymem(long long to,long long froms,long long frome) {
	if(to<froms && frome-froms>0)
		memmove(b+to*packedsize,b+froms*packedsize,(frome-froms)*packedsize*sizeof(uint32_t));
}

void savepos() {
	int i;
	/* won? */
	if(wewon) return;
	for(i=0;i<n;i++) if(state[0][i+1]!=(n-i)*2) goto not;
	for(i=0;i<n;i++) if(state[2][i+1]!=(n-i)*2-1) goto not;
	wewon=1; return;
not:
	if(curp==bblen) {
		/* WARNING, this part of the code is actually not tested, as i had enough ram
		   to avoid repacking */
		/* repack */
		curp=sortandcompress(cure,curp);
		curp=removeduplicates2(prev2,prev1,prev1,cur,cure,curp);
		/* TODO this second sort+compress can be replaced by linear time merge */
		if(cur<cure) cure=curp=sortandcompress(cur,curp);
		cure=curp;
	}
	if(curp==bblen) puts("out of memory"),exit(0);
	pack(curp++);
	if(peak<curp) peak=curp;
}

void tryallmoves() {
	int from,to,fromy,toy,maxy;
	for(from=0;from<3;from++) for(to=0;to<3;to++) if(from!=to) {
		for(fromy=0;fromy<=n;fromy++) if(state[from][fromy]) break;
		if(fromy>n) break;
		/* enforce max y-coordinate restriction */
		maxy=n-(state[from][fromy]-1)/2;
		for(toy=0;toy<=maxy;toy++) if(state[to][toy]) break;
		toy--;
		if(toy<0) continue;
		if(to==1 && !toy) continue;
		state[to][toy]=state[from][fromy]; state[from][fromy]=0;
		savepos();
		state[from][fromy]=state[to][toy]; state[to][toy]=0;
	}
}

/* breadth-first search, in memory, delayed duplicate detection */
long long bfsmemddd() {
	long long iter=1;
	long long at;
	long long tot=1;
	double tid=gettime();
	peak=0;
	while(prev1-cur) {
		if(iter%1000==0) fprintf(stderr,"  n=%d: iter %I64d, q %I64d tot %I64d, time %.3f\n",n,iter,cur-prev1,tot,gettime()-tid);
		for(curp=cure=cur,at=prev1;at<cur;at++) {
			unpack(at);
			tryallmoves();
			if(wewon) {
				fprintf(stderr,"  solved in %.3f seconds, %.0f mb, visited %I64d\n",gettime()-tid,peak*packedsize*sizeof(uint32_t)/1048576.,tot);
				return iter;
			}
		}
		curp=sortandcompress(cure,curp);
		curp=removeduplicates2(prev2,prev1,prev1,cur,cure,curp);
		/* TODO this second sort+compress can be replaced by linear time merge */
		if(cur<cure) cure=curp=sortandcompress(cur,curp);
		cure=curp;
		/* get rid of prev2 and set pointers accordingly */
		copymem(prev2,prev1,cure);
		curp-=prev1;
		prev1=cur-prev1;
		cur=curp;
		tot+=cur-prev1;
		iter++;
	}
	puts("SANITY ERROR, no solution found");
	return -1;
}

void shutdown() {
	free(b);
}

long long solve(int _n) {
	long long r;
	n=_n;
	init();
	r=bfsmemddd();
	shutdown();
	return r;
}

int main() {
	int i;
	for(i=1;i<20;i++) printf("%d %I64d\n",i,solve(i));
	return 0;
}
