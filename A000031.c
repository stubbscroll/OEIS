/* A000031 number of n-bead necklaces with 2 colours with turning over not
   allowed
   algorithm: formula Sum_{ d divides n } (phi(2d)*2^(n/d))/n
   where the divisors and phi are found using trial division */

#include <stdio.h>
#include <gmp.h>

#define MAX 1000000000

int phi(int n) {
	int r=1,i;
	for(i=2;(long long)i*i<=n;i++) if(n%i==0) {
		for(r*=i-1,n/=i;n%i==0;n/=i,r*=i);
	}
	if(n>1) r*=n-1;
	return r;
}

void calc(int d,int n,mpz_t r) {
	mpz_t t;
	mpz_init(t);
	mpz_ui_pow_ui(t,2,n/d);
	mpz_mul_si(t,t,phi(d));
	mpz_add(r,r,t);
	mpz_clear(t);
}

void calcterm(int n,mpz_t r) {
	int i;
	mpz_set_ui(r,0);
	for(i=1;(long long)i*i<=n;i++) if(n%i==0) {
		calc(i,n,r);
		if(i!=n/i) calc(n/i,n,r);
	}
	mpz_fdiv_q_ui(r,r,n);
}

int main() {
	int i;
	mpz_t r;
	mpz_init(r);
	puts("0 1"); /* algorithm doesn't like n=0, fake it */
	for(i=1;i<MAX;i++) {
		calcterm(i,r);
		gmp_printf("%d %Zd\n",i,r);
	}
	mpz_clear(r);
	return 0;
}
