/* A000005 number of divisors */

#include <stdio.h>

typedef long long ll;

/* calculate number of divisors by finding exponents in prime factorization.
   algorithm: trial division */
ll numdiv(ll n) {
	ll r=1;
	int i,j;
	for(i=2;(ll)i*i<=n;i++) if(n%i==0) {
		for(j=2,n/=i;n%i==0;n/=i,j++);
		r*=j;
	}
	if(n>1) r*=2;
	return r;
}

int main() {
	ll i;
	for(i=1;i>-1;i++) printf("%lld %lld\n",i,numdiv(i));
	return 0;
}
