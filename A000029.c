/* A000029 number of n-bead necklaces with 2 colours, allowing turning over
   algorithm: a(n) = Sum_{ d divides n } phi(d)*2^(n/d)/(2*n) +
   either 2^((n-1)/2) if n odd or 2^(n/2-1)+2^(n/2-2)-A000013(n/2-1)/2 if even
   (taken from oeis with "tweak"), rather naive implementation */

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
	mpz_mul_si(t,t,phi(2*d));
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
	mpz_fdiv_q_ui(r,r,2*n);
}

void calcterm2(int n,mpz_t r) {
	mpz_t z;
	calcterm(n,r);
	mpz_init(z);
	if(n&1) {
		mpz_ui_pow_ui(z,2,n/2);
		mpz_add(r,r,z);
	} else {
		mpz_ui_pow_ui(z,2,n/2-1);
		mpz_add(r,r,z);
		if(n>=4) {
			mpz_ui_pow_ui(z,2,n/2-2);
			mpz_add(r,r,z);
		}
		mpz_set_ui(z,0);
		calcterm(n/2,z);
		mpz_fdiv_q_ui(z,z,2);
		mpz_sub(r,r,z);
	}
	mpz_clear(z);
}

int main() {
	int i;
	mpz_t r;
	mpz_init(r);
	puts("0 1"); /* algorithm doesn't like n=0, fake it */
	for(i=1;i<MAX;i++) {
		calcterm2(i,r);
		gmp_printf("%d %Zd\n",i,r);
	}
	mpz_clear(r);
	return 0;
}
