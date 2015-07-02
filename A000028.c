/* A000028: numbers with odd number of ones in binary expansions of exponents
   in prime factorization
   algorithm: trial division */

#include <stdio.h>

typedef long long ll;

int calc(ll n) {
	int i,j,r=0;
	for(i=2;(ll)i*i<=n;i++) if(n%i==0) {
		for(j=1,n/=i;n%i==0;n/=i,j++);
		for(;j;j>>=1) if(j&1) r++;
	}
	if(n>1) r++;
	return r;
}

int main() {
	long long i=1,n=1;
	for(;n>-1;n++) if(calc(n)&1) printf("%I64d %I64d\n",i++,n);
	return 0;
}
