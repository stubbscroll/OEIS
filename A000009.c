/* A000009 partition into odd numbers, or distinct numbers
   algorithm: a(n) = P(n)-P(n-2)-P(n-4)+P(n-10)+P(n-14)+...
   +(-1)^m P(n-2p_m)+..., where P(n) is the partition function (A000041)
   and p_m =m(3m-1)/2 is the m-th generalized pentagonal number (A001318)
   (second last formula on OEIS page). this approach takes linear memory
   and O((n^1.5)^1.5) time (ignoring the size of the numbers) */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define MAX 1000001

mpz_t p[MAX];
int size;

void calcterm(int n,mpz_t r) {
	int i,j,k;
	if(!n) {
		mpz_set_si(p[0],1);
		mpz_set(r,p[0]);
		return;
	}
	/* first, calculate regular partition number */
	mpz_set_si(p[n],0);
	for(i=1;;i++) {
		for(j=1;j>-2;j-=2) {
			k=(i*j)*(3*i*j-1)/2;
			if(k>n) goto out;
			if(i&1) mpz_add(p[n],p[n],p[n-k]);
			else mpz_sub(p[n],p[n],p[n-k]);
		}
	}
out:
	/* now calculate partitions into odd numbers */
	mpz_set(r,p[n]);
	for(i=1;;i++) {
		for(j=1;j>-2;j-=2) {
			k=(i*j)*(3*i*j-1);
			if(k>n) return;
			if(i&1) mpz_sub(r,r,p[n-k]);
			else mpz_add(r,r,p[n-k]);
		}
	}
}

int main() {
	mpz_t r;
	int i;
	mpz_init(r);
	for(i=0;i<MAX;i++) mpz_init(p[i]);
	for(i=0;i<MAX;i++) {
		calcterm(i,r);
		gmp_printf("%d %Zd\n",i,r);
	}
	mpz_clear(r);
	for(i=0;i<MAX;i++) mpz_clear(p[i]);
	return 0;
}
