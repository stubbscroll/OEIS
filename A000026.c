/* A000026: mosaic numbers
   algorithm: n=prod p_i^q_j, a(n) = prod p_i*q_j, using trial division */

#include <stdio.h>

typedef long long ll;

ll calc(ll n) {
	ll r=1;
	int i,j;
	for(i=2;(ll)i*i<=n;i++) if(n%i==0) {
		for(j=1,n/=i;n%i==0;n/=i,j++);
		r*=i*j;
	}
	if(n>1) r*=n;
	return r;
}


int main() {
	long long i;
	for(i=1;i>-1;i++) printf("%I64d %I64d\n",i,calc(i));
	return 0;
}
