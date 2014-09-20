/* A000002 kolakoski sequence
   a(n) is length of a-th run, a(1)=1, only 1, 2 in sequence */

#include <stdio.h>
#include <stdlib.h>

/* 10 gigs, 80*10^9 terms. please reduce if this is too large.
   takes 214 seconds to precalculate on my i7-2600k */
#define MAX 10000000000ULL

typedef long long ll;
typedef unsigned char uchar;

uchar *a;
ll n;

#define READBIT(n) ((a[(n)>>3]&(1<<((n)&7)))>>((n)&7))
#define SETBIT(n) a[(n)>>3]|=(1<<((n)&7))

/* calculate sequence before outputting */
void precalc() {
	ll i,j;
	int flip=0,k;
	if(!(a=calloc(MAX,1))) puts("out of memory, try reducing MAX"),exit(0);
	n=MAX*8;
	SETBIT(2);
	for(i=j=2;j<n;i++) {
		k=READBIT(i)+1;
		if(flip) j+=k;
		else while(k-- && j<n) SETBIT(j),j++;
		flip^=1;
	}
}

int main() {
	ll i;
	precalc();
	for(i=1;i<n;i++) printf("%lld %d\n",i,READBIT(i)+1);
	return 0;
}
