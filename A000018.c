/* A000018 a(n) is the number of positive integers <= 2^n on the
   form x^2+16y^2.
   algorithm: stupid brute force, try all x,y. this is currently not good
   enough to generate all the terms that OEIS has, this program stops at
   term 29 with 1.6 GB RAM. a(n+1) requires roughly twice the amount of
   RAM as a(n). */

#include <stdio.h>
#include <stdlib.h>

typedef long long ll;

#define MAX 200000000
ll a[MAX];
int n;

int compl(const void *A,const void *B) {
	const ll *a=A,*b=B;
	if(*a<*b) return -1;
	if(*a>*b) return 1;
	return 0;
}

ll calc(ll N) {
	ll r=0;
	int x,y,i;
	n=0;
	for(x=0;(ll)x*x<=N;x++) {
		for(y=0;(ll)x*x+16LL*y*y<=N;y++) {
			if(n==MAX) puts("out of memory"),exit(0);
			a[n++]=(ll)x*x+16LL*y*y;
		}
	}
	qsort(a,n,sizeof(long long),compl);
	for(i=1;i<n;i++) if(a[i]!=a[i-1]) r++;
	return r;
}

int main() {
	int i;
	for(i=0;i<63;i++) printf("%d %lld\n",i,calc(1LL<<i));
	return 0;
}
