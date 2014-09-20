/* A000051 2^n+1
   algorithm: a(0)=2, a(n)=2*a(n-1)-1 */

#include <stdio.h>
#include <gmp.h>

int main() {
	int i;
	mpz_t r;
	mpz_init_set_si(r,2);
	for(i=0;i>-1;i++) {
		gmp_printf("%d %Zd\n",i,r);
		mpz_mul_si(r,r,2);
		mpz_sub_ui(r,r,1);
	}
	mpz_clear(r);
	return 0;
}
