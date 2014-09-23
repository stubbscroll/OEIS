/* A000020 number of primitive polynomials of degree n over GF(2),
   except for the erroneous a(1)=2.
   formula: phi(2^n-1)/n.
   algorithm: start with N=2^n-1. we have that when d|n, 2^d-1|2^n-1.
   for all d|n, factorize gcd(N,2^d-1) and let N=gcd(N,2^d-1).
   we try the following factorization algorithms in order:
   - trial division
   - pollard-rho
   - pollard p-1
   - elliptic curve method. SLOW and in dire need of optimization or
     discovering what's wrong. alpertron destroys my implementation
     speed-wise
   - (more to come later, qs will definitely be added)
   then process the list of prime factors and calculate phi.
   this version finds 256 of 400 terms, but it had to run overnight.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

/* limit for trial division */
#define SMALL 1000000

/* limit (number of digits) for pollard-rho. TODO adjust. the routine takes
   a long while on as little as 33 digits */
#define POLLARDLIM 32

#define MAXF 10000
mpz_t factors[MAXF];
int fn;

#define ERROR puts("error, too many factors. increase MAXF and recompile"),exit(0)

gmp_randstate_t gmprand;

/* start of factorization routines! */

int isprime(int n) {
	int i;
	if(n<4) return n>1;
	if(!(n&1)) return 0;
	for(i=3;i*i<=n;i++) if(n%i==0) return 0;
	return 1;
}

/* start of elliptic curve method! implementation based on "factorization and
   primality testing" (bressaud) */

/* helpful parameters for ECM (B1, number of curves)
   http://www.alpertron.com.ar/ECM.HTM
*/

void sub_2i(mpz_t r,mpz_t s,mpz_t n,mpz_t a,mpz_t b,mpz_t X2,mpz_t Z2) {
	mpz_t t,u,v;
	mpz_init(t); mpz_init(u); mpz_init(v);
	/* calculate X2 */
	mpz_mul(t,r,r); mpz_mul(u,s,s); mpz_mul(u,u,a);
	mpz_sub(t,t,u); mpz_mod(t,t,n);
	mpz_mul(t,t,t);
	mpz_pow_ui(u,s,3); mpz_mul(u,u,r); mpz_mul(u,u,b); mpz_mul_ui(u,u,8);
	mpz_sub(t,t,u); mpz_mod(v,t,n);
	/* calculate Z2 */
	mpz_pow_ui(t,r,3); mpz_mul(u,a,r); mpz_mul(u,u,s); mpz_mul(u,u,s);
	mpz_add(t,t,u);
	mpz_pow_ui(u,s,3); mpz_mul(u,u,b); mpz_add(t,t,u); mpz_mod(t,t,n);
	mpz_mul(t,s,t); mpz_mul_ui(t,t,4); mpz_mod(Z2,t,n); mpz_set(X2,v);
	mpz_clear(t); mpz_clear(u); mpz_clear(v);
}

void sub_2i_plus(mpz_t r,mpz_t s,mpz_t u,mpz_t v,mpz_t X,mpz_t Z,mpz_t n,mpz_t a,mpz_t b,mpz_t U1,mpz_t U2) {
	mpz_t t,t1,t2,u1;
	mpz_init(t); mpz_init(t1); mpz_init(t2); mpz_init(u1);
	/* obtain U1 */
	mpz_mul(t1,r,u); mpz_mul(t,a,s); mpz_mul(t,t,v); mpz_sub(t1,t1,t); mpz_mod(t1,t1,n);
	mpz_mul(t2,r,v); mpz_mul(t,s,u); mpz_add(t2,t2,t); mpz_mul(t2,t2,v);
	mpz_mul(t2,t2,s); mpz_mul(t2,t2,b); mpz_mod(t2,t2,n);
	mpz_mul(t1,t1,t1); mpz_mul_ui(t2,t2,4); mpz_sub(t1,t1,t2); mpz_mul(t1,t1,Z);
	mpz_mod(u1,t1,n);
	/* obtain U2 */
	mpz_mul(t,u,s); mpz_mul(t1,r,v); mpz_sub(t,t,t1); mpz_mod(t,t,n);
	mpz_mul(t,t,t); mpz_mul(t,t,X); mpz_mod(U2,t,n);
	mpz_set(U1,u1);
	mpz_clear(t); mpz_clear(t1); mpz_clear(t2); mpz_clear(u1);
}

/* multiply (X,Z) with k in y^2=x^3+ax+b mod n */
void nextvalues(mpz_t X,mpz_t Z,int p,mpz_t n,mpz_t a,mpz_t b) {
	mpz_t X1,Z1,X2,Z2,U1,U2,t;
	int len=0,c[32],i;
	mpz_init(U1); mpz_init(U2); mpz_init(t);
	mpz_init_set(X1,X); mpz_init_set(Z1,Z); mpz_init(X2); mpz_init(Z2);
	while(p) c[len++]=p&1,p>>=1;
	sub_2i(X,Z,n,a,b,X2,Z2);
	/* len-1 or len-2? */
	for(i=len-2;i>=0;i--) {
		sub_2i_plus(X1,Z1,X2,Z2,X,Z,n,a,b,U1,U2);
		if(!c[i]) {
			sub_2i(X1,Z1,n,a,b,t,Z1);
			mpz_set(X1,t); mpz_set(X2,U1); mpz_set(Z2,U2);
		} else {
			sub_2i(X2,Z2,n,a,b,t,Z2);
			mpz_set(X2,t); mpz_set(X1,U1); mpz_set(Z1,U2);
		}
	}
	mpz_set(X,X1); mpz_set(Z,Z1);
	mpz_clear(U1); mpz_clear(U2); mpz_clear(t);
	mpz_clear(X1); mpz_clear(Z1); mpz_clear(X2); mpz_clear(Z2);
}

/* elliptic curve method! returns 1 and factor in out if factor is found,
   otherwise return 0, maxb is max prime, maxc is number of curves to test */
int ECM(mpz_t n,mpz_t out,int maxb,int maxc) {
	int r=0,p;
	mpz_t X,Y,Z,a,b,g,t;
	mpz_init(X); mpz_init(Y); mpz_init(Z); mpz_init(a); mpz_init(b); mpz_init(g); mpz_init(t);
	gmp_printf("start ecm on %Zd with b1 %d curves %d\n",n,maxb,maxc);
	while(maxc--) {
	newcurve:
		/* pick random curve */
		mpz_urandomm(X,gmprand,n); mpz_urandomm(Y,gmprand,n); mpz_urandomm(a,gmprand,n);
		/* determine b=(y^2-x^3-ax)%n */
		mpz_mul(b,Y,Y); mpz_pow_ui(t,X,3); mpz_sub(b,b,t);
		mpz_mul(t,a,X); mpz_sub(b,b,t); mpz_mod(b,b,n);
//		gmp_printf("curve a %Zd\n      b %Zd\n      x %Zd\n      y %Zd\n",a,b,X,Y);
		/* check gcd(4a^3+27b^2,n) */
		mpz_pow_ui(g,a,3); mpz_mul_si(g,g,4); mpz_mul(t,b,b);
		mpz_mul_si(t,t,27); mpz_add(g,g,t); mpz_gcd(g,g,n);
		if(!mpz_cmp(g,n)) goto newcurve;
		if(mpz_cmp_ui(g,1)>0) {
			/* factor found! */
			r=1;
			mpz_set(a,g);
			goto end;
		}
		mpz_set_ui(Z,1);
		for(p=2;p<=maxb;p++) {
			nextvalues(X,Z,p,n,a,b);
			if(!(p&15)) {
				mpz_gcd(g,Z,n);
				if(!mpz_cmp(g,n)) goto fail;
				if(mpz_cmp_ui(g,1)>0) {
					r=1;
					mpz_set(out,g);
					goto end;
				}
			}
		}
		/* loop done, check again */
		mpz_gcd(g,Z,n);
		if(mpz_cmp_ui(g,1)>0) {
			r=1;
			mpz_set(out,g);
			goto end;
		}
	fail:;
	}
end:
	mpz_clear(X); mpz_clear(Z); mpz_clear(a); mpz_clear(b); mpz_clear(g); mpz_clear(t);
	return r;
}

/* pollard p-1 with stage 2 with large prime. hope that n has a factor p such
   that p-1 is B1-smooth with at most one factor larger than B1 and not larger
   than B2. return 0 if it fails, otherwise return 1 and return factor in a. */
int pollardp1(mpz_t n,mpz_t a,int B1,int B2,int maxc) {
	int b,c,r=0,q;
	mpz_t m,g;
	mpz_init(m); mpz_init(g);
	for(c=2;c<2+maxc;c++) {
		mpz_set_ui(m,c);
		for(b=2;b<=B1;b++) {
			mpz_powm_ui(m,m,b,n);
			/* check for factor periodically */
			if(!(b&1023)) {
				mpz_sub_ui(g,m,1);
				mpz_gcd(g,g,n);
				if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) {
					r=1;
					mpz_set(a,g);
					goto end;
				}
			}
		}
		/* we're done, check again */
		mpz_sub_ui(g,m,1);
		mpz_gcd(g,g,n);
		if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) {
			r=1;
			mpz_set(a,g);
			goto end;
		}
		/* stage 2! for each prime p between B1 and B2, check p*m */
		if(!(b&1)) b++;
		for(q=0;b<=B2;b+=2) if(isprime(b)) {
			mpz_powm_ui(m,m,b-q,n); q=b;
			mpz_sub_ui(g,m,1);
			mpz_gcd(g,g,n);
			if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) {
				r=1;
				mpz_set(a,g);
				goto end;
			}
		}
	}
end:
	mpz_clear(m); mpz_clear(g);
	return r;
}

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
		mpz_set_ui(range,1);
		mpz_set_ui(x1,2);
		mpz_set_ui(x2,4); mpz_add_ui(x2,x2,c); mpz_mod(x2,x2,n);
		mpz_set_ui(product,1);
		for(j=0;mpz_cmp_ui(max,0);mpz_sub_ui(max,max,1)) {
			for(mpz_set_ui(k,0);mpz_cmp(k,range)<0;mpz_add_ui(k,k,1)) {
				mpz_powm_ui(x2,x2,2,n); mpz_add_ui(x2,x2,c);
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
				mpz_powm_ui(x2,x2,2,n); mpz_add_ui(x2,x2,c);
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
	}
	/* always try pollard p-1 */
	mpz_init(m);
	if(pollardp1(n,m,500000,500000,2)) {
		factorize2(m);
		mpz_divexact(m,n,m);
		factorize2(m);
		mpz_clear(m);
		return;
	}
	/* ecm */
	if(ECM(n,m,2000,25)) goto found;
	if(ECM(n,m,11000,90)) goto found;
	if(ECM(n,m,50000,300)) goto found;
	if(ECM(n,m,250000,300)) goto found;
	goto wrong;
found:
	factorize2(m);
	mpz_divexact(m,n,m);
	factorize2(m);
	mpz_clear(m);
	return;
wrong:
//	gmp_printf("TODO large number (%d digits)\n",d);
	gmp_printf("todo %Zd (%d digits)\n",n,d);
	if(fn==MAXF) ERROR;
	mpz_set(factors[fn++],n);
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
	gmp_randinit_mt(gmprand);
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
