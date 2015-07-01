/* A000023 E.g.f.: exp(-2*x)/(1-x), whatever that is
   algorithm: a(0)=1, a(n)=a(n-1)*n+(-2)^n */

#include <stdio.h>
#include <gmp.h>

int main() {
	int i;
	mpz_t a,p2;
	mpz_init_set_si(a,1);
	mpz_init_set_si(p2,2);
	gmp_printf("0 %Zd\n",a);
	for(i=1;i>-1;i++) {
		mpz_mul_si(a,a,i);
		if(i&1) mpz_sub(a,a,p2);
		else mpz_add(a,a,p2);
		gmp_printf("%d %Zd\n",i,a);
		mpz_mul_si(p2,p2,2);
	}
	mpz_clear(p2);
	mpz_clear(a);
	return 0;
}
