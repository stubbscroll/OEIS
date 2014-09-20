/* A000032 lucas numbers */

#include <stdio.h>
#include <gmp.h>

int main() {
	mpz_t a,b;
	int i;
	mpz_init_set_si(a,2);
	mpz_init_set_si(b,1);
	gmp_printf("0 %Zd\n",a);
	gmp_printf("1 %Zd\n",b);
	for(i=2;i>-1;i+=2) {
		mpz_add(a,a,b);
		mpz_add(b,b,a);
		gmp_printf("%d %Zd\n",i,a);
		gmp_printf("%d %Zd\n",i+1,b);
	}
	return 0;
}
