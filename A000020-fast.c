/* A000020 number of primitive polynomials of degree n over GF(2),
   except for the erroneous a(1)=2.
   formula: phi(2^n-1)/n.
   algorithm: start with N=2^n-1. we have that when d|n, 2^d-1|2^n-1.
   for all d|n, factorize gcd(N,2^d-1) and let N=gcd(N,2^d-1).
   whenever a huge number needs to be factorized, bring out the heavy
   machinery (pollard-rho implemented so far).
   then process the list of prime factors and calculate phi.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

/* limit for trial division */
#define SMALL 1000000

/* limit (number of digits) for pollard-rho. TODO adjust. the routine takes
   a long while on as little as 33 digits */
/* we haven't implemented anything more powerful than pollard-rho yet,
   so run it on everything */
//#define POLLARDLIM 36
#define POLLARDLIM 10000

#define MAXF 10000
mpz_t factors[MAXF];
int fn;

#define ERROR puts("error, too many factors. increase MAXF and recompile"),exit(0)

/* start of factorization routines! */

/* pollard-rho! when a factor is found, return it in a. the routine will run
   until a factor is found. it is the caller's responsibility to call this
   function on composite values that aren't too huge. */
void pollardrho(mpz_t n,mpz_t a) {
	int c,j;
	mpz_t max,range,k,x1,x2,product,temp,g;
	mpz_init(max); mpz_init(range); mpz_init(k); mpz_init(x1);
	mpz_init(x2); mpz_init(product); mpz_init(temp); mpz_init(g);
	/* try x^2+c for increasing c */
	for(c=1;;c++) {
		/* birthday! probability 0.5 to find a factor p after 1.179*sqrt(p)
		   iterations. set the max number of iterations slightly higher before
		   trying another c */
		mpz_mul_ui(max,n,2);
		mpz_root(max,max,4);
//		gmp_fprintf(stderr,"pollard rho n=%Zd max=%Zd c=%d, %d digits\n",n,max,c,mpz_sizeinbase(n,10));
		mpz_set_ui(range,1);
		mpz_set_ui(x1,2);
		mpz_set_ui(x2,4); mpz_add_ui(x2,x2,c); mpz_mod(x2,x2,n);
		mpz_set_ui(product,1);
		for(j=0;mpz_cmp_ui(max,0);mpz_sub_ui(max,max,1)) {
			for(mpz_set_ui(k,0);mpz_cmp(k,range)<0;mpz_add_ui(k,k,1)) {
				mpz_powm_ui(x2,x2,2,n);
				mpz_add_ui(x2,x2,c);
				if(mpz_cmp(x2,n)>-1) mpz_sub(x2,x2,n);
				mpz_sub(temp,x1,x2);
				mpz_abs(temp,temp);
				mpz_mul(product,product,temp); mpz_mod(product,product,n);
				j++;
				if(!(j&7)) {
					mpz_gcd(g,product,n);
					if(!mpz_cmp(g,n)) goto failed;
					if(mpz_cmp_si(g,1)>0) {
						mpz_set(a,g);
						goto done;
					}
					mpz_set_ui(product,1);
				}
			}
			mpz_set(x1,x2);
			mpz_mul_2exp(range,range,1);
			for(mpz_set_ui(k,0);mpz_cmp(k,range)<0;mpz_add_ui(k,k,1)) {
				mpz_powm_ui(x2,x2,2,n);
				mpz_add_ui(x2,x2,c);
				if(mpz_cmp(x2,n)>-1) mpz_sub(x2,x2,n);
			}
		}
	failed:;
	}
done:
	mpz_clear(max); mpz_clear(range); mpz_clear(k); mpz_clear(x1);
	mpz_clear(x2); mpz_clear(product); mpz_clear(temp); mpz_clear(g);
}

/* "black fox"-factorization for n with no prime factors smaller than SMALL.
   pollard-rho and TODO depending on number of digits */
void factorize2(mpz_t n) {
	int d;
	mpz_t m;
	if(mpz_probab_prime_p(n,200)) {
		if(fn==MAXF) ERROR;
		mpz_set(factors[fn++],n);
		return;
	}
	d=mpz_sizeinbase(n,10);
	if(d<=POLLARDLIM) {
		mpz_init(m);
		pollardrho(n,m);
		factorize2(m);
		mpz_divexact(m,n,m);
		factorize2(m);
		mpz_clear(m);
		return;
	} else {
		gmp_printf("TODO large number (%d digits)\n",d);
//		gmp_printf("todo heavier machinery needed to factorize %Zd (%d digits)\n",n,d);
		if(fn==MAXF) ERROR;
		mpz_set(factors[fn++],n);
		return;
	}
}

/* "black box" factorization that lumps all factors in factors[] in no
   particular order */
void factorize(mpz_t in) {
	int i;
	mpz_t n,ii;
	/* in prime? */
	if(mpz_probab_prime_p(in,200)) {
		mpz_set(factors[fn++],in);
		return;
	}
	/* don't modify in */
	mpz_init_set(n,in);
	mpz_init(ii);
	/* first, trial division with small numbers */
	for(i=3;i<SMALL;i+=2) {
		mpz_set_ui(ii,i);
		mpz_mul_ui(ii,ii,i);
		if(mpz_cmp(ii,n)>0) break;
		if(mpz_divisible_ui_p(n,i)) {
			mpz_divexact_ui(n,n,i);
			if(fn==MAXF) ERROR;
			mpz_set_ui(factors[fn++],i);
			while(mpz_divisible_ui_p(n,i)) {
				mpz_divexact_ui(n,n,i);
				if(fn==MAXF) ERROR;
				mpz_set_ui(factors[fn++],i);
			}
			if(mpz_probab_prime_p(n,200)) break;
		}
	}
	if(mpz_probab_prime_p(n,200)) {
		if(fn==MAXF) ERROR;
		mpz_set(factors[fn++],n);
	}
	else if(mpz_cmp_ui(n,1)>0) factorize2(n);
	mpz_clear(n); mpz_clear(ii);
}

void phi(mpz_t r) {
	int i,j;
	mpz_t t;
	mpz_init(t);
	/* sort factors. get away with bubble sort for now */
	for(i=0;i<fn-1;i++) for(j=0;j<fn-1;j++) if(mpz_cmp(factors[j],factors[j+1])>0) {
		mpz_set(t,factors[j]);
		mpz_set(factors[j],factors[j+1]);
		mpz_set(factors[j+1],t);
	}
	/* evaluate phi */
	mpz_set_ui(r,1);
	if(fn) {
		mpz_set(r,factors[0]);
		mpz_sub_ui(r,r,1);
	}
	for(i=1;i<fn;i++) if(!mpz_cmp(factors[i-1],factors[i])) mpz_mul(r,r,factors[i]);
	else mpz_set(t,factors[i]),mpz_sub_ui(t,t,1),mpz_mul(r,r,t);
	mpz_clear(t);
}

void calc(int n,mpz_t r) {
	int d;
	mpz_t N,D,rem;
	mpz_init(N);
	mpz_init(D);
	mpz_init(rem);
	/* obtain N=2^n-1 */
	mpz_ui_pow_ui(N,2,n);
	mpz_sub_ui(N,N,1);
	/* for each d dividing N, factorize 2^d-1 and remove factors from N */
	fn=0;
	for(d=2;d*d<=n;d++) if(n%d==0) {
		/* obtain D=2^d-1 */
		mpz_ui_pow_ui(D,2,d);
		mpz_sub_ui(D,D,1);
		/* take gcd because of overlapping divisors */
		mpz_gcd(D,N,D);
		if(mpz_cmp_ui(D,1)>0) {
			factorize(D);
			mpz_divexact(N,N,D);
		}
		/* obtain D=2^(n/d)-1 and do the same as above */
		if(d!=n/d) {
			mpz_ui_pow_ui(D,2,n/d);
			mpz_sub_ui(D,D,1);
			mpz_gcd(D,N,D);
			if(mpz_cmp_ui(D,1)>0) {
				factorize(D);
				mpz_divexact(N,N,D);
			}
		}
	}
	/* factorize what remains */
	if(mpz_cmp_ui(N,1)>0) factorize(N);
	phi(r);
	mpz_divexact_ui(r,r,n);
	mpz_clear(N);
	mpz_clear(D);
	mpz_init(rem);
}

int main() {
	int i;
	mpz_t r;
	mpz_init(r);
	for(i=0;i<MAXF;i++) mpz_init(factors[i]);
	puts("1 2"); /* wrong term */
	for(i=2;;i++) {
		calc(i,r);
		gmp_printf("%d %Zd\n",i,r);
	}
	mpz_clear(r);
	for(i=0;i<MAXF;i++) mpz_clear(factors[i]);
	return 0;
}
