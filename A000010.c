/* A000010 euler's totient (phi) function
   algorithm: trial division and a formula equivalent to the standard one
   which is n*Product_{distinct primes p dividing n} (1-1/p) */

#include <stdio.h>

int phi(int n) {
	int r=1,i;
	for(i=2;(long long)i*i<=n;i++) if(n%i==0) {
		for(r*=i-1,n/=i;n%i==0;n/=i,r*=i);
	}
	if(n>1) r*=n-1;
	return r;
}

int main() {
	int i;
	for(i=1;i>-1;i++) printf("%d %d\n",i,phi(i));
	return 0;
}
