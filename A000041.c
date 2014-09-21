/* A000041: number of unordered partitions of n into integers
   formula: a(n) = a(n-1) + a(n-2) - a(n-5) - a(n-7) + a(n-12) + ...
   where i in a(n-i) is from the sequence n(3n-1)/2 with n=+-1, +-2, +-3, ...
   and the signs follow the sequence ++--++--...
   source: third formula at http://oeis.org/A000041
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define MAX 1000001

/* dynamically growing an array of mpz_t with realloc somehow breaks all
   numbers in the array. using static array instead */
mpz_t a[MAX];
int size;

void calcterm(int n,mpz_t r) {
	int i,j,k;
	if(!n) {
		mpz_set_si(a[0],1);
		mpz_set(r,a[0]);
		return;
	}
	mpz_set_si(a[n],0);
	for(i=1;;i++) {
		for(j=1;j>-2;j-=2) {
			k=(i*j)*(3*i*j-1)/2;
			if(k>n) goto out;
			if(i&1) mpz_add(a[n],a[n],a[n-k]);
			else mpz_sub(a[n],a[n],a[n-k]);
		}
	}
out:
	mpz_set(r,a[n]);
}

int main() {
	mpz_t r;
	int i;
	mpz_init(r);
	for(i=0;i<MAX;i++) mpz_init(a[i]);
	for(i=0;i<MAX;i++) {
		calcterm(i,r);
		gmp_printf("%d %Zd\n",i,r);
	}
	mpz_clear(r);
	for(i=0;i<MAX;i++) mpz_clear(a[i]);
	return 0;
}
