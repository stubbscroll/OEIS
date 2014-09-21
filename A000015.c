/* A000015 smallest prime power >= n
   for n, find the smallest p^k>=n where p is prime and k>=1 */

#include <stdio.h>
#include <math.h>

/* naive primality test */
int isprime(int n) {
	int i;
	if(n<4) return n>1;
	if(!(n&1)) return 0;
	for(i=3;(long long)i*i<=n;i++) if(n%i==0) return 0;
	return 1;
}

/* evaluate n^k */
long long intpow(int n,int k) {
	long long r=1;
	int i;
	for(i=0;i<k;i++) r*=n;
	return r;
}

/* find integer part of k-th root of n */
int introot(int n,int k) {
	int b=pow(n,1./k+1e-9);
	/* don't trust floating point math! adjust if too low */
	while(intpow(b+1,k)<=n) b++;
	/* adjust if too high */
	while(intpow(b-1,k)>=n) b--;
	return b;
}

long long smallestpower(int n) {
	long long r,t;
	int k,b;
	/* special case: p^0 */
	if(n==1) return 1;
	/* find smallest prime >= n */
	for(r=n;!isprime(r);r++);
	/* try each power until 2^power>=smallest so far */
	for(k=2;;k++) {
		if((1<<k)>=r) break;
		b=introot(n,k);
		while(!isprime(b)) b++;
		while(n>intpow(b,k)) {
			b++;
			while(!isprime(b)) b++;
		}
		t=intpow(b,k);
		if(r>t) r=t;
	}
	return r;
}

int main() {
	int i;
	for(i=1;i>-1;i++) printf("%d %lld\n",i,smallestpower(i));
	return 0;
}
