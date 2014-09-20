/* A000044 dying rabbits */

#include <stdio.h>
#include <gmp.h>

mpz_t last[14];

int main() {
	int i,j;
	mpz_init_set_si(last[0],1);
	mpz_init_set_si(last[1],1);
	mpz_init_set_si(last[2],1);
	mpz_init_set_si(last[3],2);
	mpz_init_set_si(last[4],3);
	mpz_init_set_si(last[5],5);
	mpz_init_set_si(last[6],8);
	mpz_init_set_si(last[7],13);
	mpz_init_set_si(last[8],21);
	mpz_init_set_si(last[9],34);
	mpz_init_set_si(last[10],55);
	mpz_init_set_si(last[11],89);
	mpz_init_set_si(last[12],144);
	mpz_init(last[13]);
	for(i=0;i<13;i++) gmp_printf("%d %Zd\n",i,last[i]);
	for(i=13;i>-1;i++) {
		mpz_set(last[13],last[12]);
		mpz_add(last[13],last[13],last[11]);
		mpz_sub(last[13],last[13],last[0]);
		for(j=0;j<13;j++) mpz_set(last[j],last[j+1]);
		gmp_printf("%d %Zd\n",i,last[13]);
	}
	return 0;
}
