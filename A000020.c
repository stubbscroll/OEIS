/* A000020 number of primitive polynomials of degree n over GF(2),
   except for the erroneous a(1)=2.
   algorithm: phi(2^n-1)/n.
   naive factorization is used which becomes very slow after 60-70 terms or so.
*/

#include <stdio.h>
#include <gmp.h>

/* calculate euler's totient (phi) for bigint */
void bigphi(mpz_t n,mpz_t r) {
	mpz_t i,ii,m,t;
	mpz_init(i); mpz_init(ii); mpz_init(m); mpz_init(t);
	mpz_set_ui(r,1);
	if(!mpz_probab_prime_p(n,200)) {
		for(mpz_set_ui(i,2);;mpz_add_ui(i,i,1)) {
			mpz_mul(ii,i,i);
			if(mpz_cmp(ii,n)>0) break;
			mpz_fdiv_r(m,n,i);
			if(!mpz_cmp_ui(m,0)) {
				mpz_sub_ui(t,i,1);
				mpz_mul(r,r,t);
				mpz_fdiv_q(n,n,i);
				while(1) {
					mpz_fdiv_qr(t,m,n,i);
					if(mpz_cmp_ui(m,0)) break;
					mpz_set(n,t);
					mpz_mul(r,r,i);
				}
				/* prune if remaining n is prime */
				if(mpz_probab_prime_p(n,200)) break;
			}
		}
	}
	if(mpz_cmp_ui(n,1)>0) {
		mpz_sub_ui(n,n,1);
		mpz_mul(r,r,n);
	}
	mpz_clear(i); mpz_clear(ii); mpz_clear(m); mpz_clear(t);
}

void calc(int n,mpz_t r) {
	mpz_t n2,phi;
	mpz_init(n2);
	mpz_init(phi);
	mpz_ui_pow_ui(n2,2,n);
	mpz_sub_ui(n2,n2,1);
	bigphi(n2,phi);
	mpz_fdiv_q_ui(r,phi,n);
	mpz_clear(phi);
	mpz_clear(n2);
}

int main() {
	int i;
	mpz_t r;
	mpz_init(r);
	puts("1 2"); /* wrong term */
	for(i=2;;i++) {
		calc(i,r);
		gmp_printf("%d %Zd\n",i,r);
	}
	mpz_clear(r);
	return 0;
}
