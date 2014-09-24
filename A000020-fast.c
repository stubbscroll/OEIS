/* A000020 number of primitive polynomials of degree n over GF(2),
   except for the erroneous a(1)=2.
   formula: phi(2^n-1)/n.
   algorithm: start with N=2^n-1. we have that when d|n, 2^d-1|2^n-1.
   for all d|n, factorize gcd(N,2^d-1) and let N=gcd(N,2^d-1).
   if n is prime, then all prime factors must be of the form 2kn+1 for
   integer k>=1 and must also be 1 or 7 mod 8.
   the following factorization algorithms are tried in order:
   - trial division
   - pollard-rho
   - pollard p-1
   - elliptic curve method (speed is roughly 66% of alpertron, i think)
   - (more to come later, qs will definitely be added)
   - (maybe also add a check for primes of the form 2kn+1)
   then process the list of prime factors and calculate phi.
   this is very tunable. ecm is much slower than the other methods, so
   pollard-rho and p-1 could be given more time.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

/* misc routines, such as sieve */

#define MAXP 100000000
typedef long long ll;

char sieve[MAXP];
int prime[5761456];
int primes;

void createsieve() {
	int i,j,q;
	memset(sieve,1,sizeof(sieve));
	q=sqrt(MAXP);
	for(sieve[0]=sieve[1]=0,i=2;i<=q;i++) if(sieve[i]) for(j=i*i;j<MAXP;j+=i) sieve[j]=0;
}

void genprimes() {
	int i;
	for(primes=i=0;i<MAXP;i++) if(sieve[i]) prime[primes++]=i;
}

/* start of ECM routines */

/* helpful parameters for ECM (B1, number of curves)
   http://www.alpertron.com.ar/ECM.HTM
*/

/* global variables for ECM */
struct {
	mpz_t n,C;
	/* temp variables */
	mpz_t t,t1,t2,t3,t4;
	mpz_t U,V,T,W;
} ecm;

void addh(mpz_t X1,mpz_t Z1,mpz_t X2,mpz_t Z2,mpz_t X,mpz_t Z,mpz_t Xo,mpz_t Zo) {
	/* Xo = Z*(X1*X2 - Z1*Z2)^2 */
	mpz_mul(ecm.t1,X1,X2); mpz_submul(ecm.t1,Z1,Z2); mpz_mul(ecm.t1,ecm.t1,ecm.t1);
	mpz_mul(ecm.t1,ecm.t1,Z);
	/* Zo = X*(X1*Z2 - X2*Z1)^2 */
	mpz_mul(ecm.t2,X1,Z2); mpz_submul(ecm.t2,X2,Z1); mpz_mul(ecm.t2,ecm.t2,ecm.t2);
	mpz_mul(Zo,ecm.t2,X); mpz_mod(Zo,Zo,ecm.n);
	mpz_mod(Xo,ecm.t1,ecm.n);
}

/* double X,Y and put result in X2,Z2 */
void doubleh(mpz_t X,mpz_t Z,mpz_t X2,mpz_t Z2) {
	/* t1=X*X, t2=Z*Z */
	mpz_mul(ecm.t1,X,X); mpz_mul(ecm.t2,Z,Z);
	/* Z2 = 4*Z*(X*X*X + C*X*X*Z + X*Z*Z), build inner stuff in t3 */
	mpz_mul(ecm.t3,ecm.t1,X); mpz_mul(ecm.t4,ecm.t1,ecm.C); mpz_mul(ecm.t4,ecm.t4,Z);
	mpz_add(ecm.t3,ecm.t3,ecm.t4); mpz_mul(ecm.t4,ecm.t2,X);
	mpz_add(ecm.t3,ecm.t3,ecm.t4);
	mpz_mul(ecm.t3,Z,ecm.t3); mpz_mul_ui(Z2,ecm.t3,4); mpz_mod(Z2,Z2,ecm.n);
	/* X2=(X*X-Z*Z)^2 */
	mpz_sub(X2,ecm.t1,ecm.t2); mpz_mul(X2,X2,X2); mpz_mod(X2,X2,ecm.n);
}

void multiply(mpz_t X,mpz_t Z,int p,mpz_t X2,mpz_t Z2) {
	int b;
	if(p<2) puts("error not implemented"),exit(0);
	if(p==2) return doubleh(X,Z,X2,Z2);
	mpz_set(ecm.U,X); mpz_set(ecm.V,Z);
	doubleh(X,Z,ecm.T,ecm.W);
	for(b=30;;b--) if(p&(1<<b)) break;
	for(b--;b>=0;b--) {
		if(p&(1<<b)) {
			addh(ecm.T,ecm.W,ecm.U,ecm.V,X,Z,ecm.U,ecm.V);
			doubleh(ecm.T,ecm.W,ecm.T,ecm.W);
		} else {
			addh(ecm.U,ecm.V,ecm.T,ecm.W,X,Z,ecm.T,ecm.W);
			doubleh(ecm.U,ecm.V,ecm.U,ecm.V);
		}
	}
	if(p&1) return addh(ecm.U,ecm.V,ecm.T,ecm.W,X,Z,X2,Z2);
	doubleh(ecm.U,ecm.V,X2,Z2);
}

#define D 100

/* faster ECM from "prime numbers - a computational perspective" (crandall,
   pomerance), algorithm 7.4.4 */
/* return 1 and factor in out if factor is found, otherwise return 0.
   B1 is max prime (must be even), maxc is number of curves to test */
int ECM(mpz_t n,mpz_t out,int B1,int maxc) {
	long long q;
	int r=0,B2=100*B1,i,B,j,delta;
	mpz_t sigma,u,v,Qx,Qz,g,t,t1;
	mpz_t Sx[D],Sz[D],beta[D],Tx,Tz,Rx,Rz,alpha;
	mpz_init(t); mpz_init(t1);
	mpz_init(sigma); mpz_init(u); mpz_init(v); mpz_init(Qx); mpz_init(Qz); mpz_init(g);
	/* global variables */
	mpz_init_set(ecm.n,n); mpz_init(ecm.C);
	mpz_init(ecm.U); mpz_init(ecm.V); mpz_init(ecm.T); mpz_init(ecm.W);
	mpz_init(ecm.t); mpz_init(ecm.t1); mpz_init(ecm.t2); mpz_init(ecm.t3);
	while(maxc--) {
		/* choose random curve */
		do {
			mpz_urandomm(sigma,gmprand,n); /* sigma between 6 and n-1 */
		} while(mpz_cmp_si(sigma,6)<0);
		/* u=(sigma^2-5) mod n */
		mpz_mul(u,sigma,sigma); mpz_sub_ui(u,u,5); mpz_mod(u,u,n);
		/* v=4*sigma mod n */
		mpz_mul_si(v,sigma,4); mpz_mod(v,v,n);
		/* C=((u-v)^3)(3u+v)/(4u^3v)-2) mod n */
		mpz_sub(t,v,u); mpz_powm_ui(t,t,3,n);
		mpz_mul_ui(t1,u,3); mpz_add(t1,t1,v); mpz_mul(t,t,t1); mpz_mod(t,t,n);
		mpz_powm_ui(t1,u,3,n); mpz_mul(t1,t1,v); mpz_mul_si(t1,t1,4); mpz_mod(t1,t1,n);
		mpz_invert(t1,t1,n); mpz_mul(t,t,t1); mpz_sub_ui(t,t,2); mpz_mod(ecm.C,t,n);
		/* Qx=u^3 mod n, Qy=v^3 mod n */
		mpz_powm_ui(Qx,u,3,n); mpz_powm_ui(Qz,v,3,n);
		/* perform stage 1 */
		for(j=0;j<primes && prime[j]<B1;j++) {
			for(q=1;q<=B1;q*=prime[j]) multiply(Qx,Qz,prime[j],Qx,Qz);
		}
		mpz_gcd(g,Qz,n);
		if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) {
			r=1;
			mpz_set(out,g);
			goto end;
		}
		/* perform stage 2 */
		for(i=0;i<D;i++) mpz_init(Sx[i]),mpz_init(Sz[i]),mpz_init(beta[i]);
		mpz_init(Tx); mpz_init(Tz); mpz_init(Rx); mpz_init(Rz); mpz_init(alpha);
		doubleh(Qx,Qz,Sx[0],Sz[0]);
		doubleh(Sx[0],Sz[0],Sx[1],Sz[1]);
		for(i=1;i<=D;i++) {
			if(i>2) addh(Sx[i-2],Sz[i-2],Sx[0],Sz[0],Sx[i-3],Sz[i-3],Sx[i-1],Sz[i-1]);
			mpz_mul(beta[i-1],Sx[i-1],Sz[i-1]); mpz_mod(beta[i-1],beta[i-1],n);
		}
		mpz_set_ui(g,1);
		B=B1-1;
		multiply(Qx,Qz,B-2*D,Tx,Tz);
		multiply(Qx,Qz,B,Rx,Rz);
		for(i=B;i<B2;i+=2*D) {
			mpz_mul(alpha,Rx,Rz); mpz_mod(alpha,alpha,n);
			for(;j<primes && prime[j]<=i+2*D;j++) {
				delta=(prime[j]-i)/2-1;
				if(delta>=D || delta<0) printf("error p %d i %d delta %d error\n",prime[j],i,delta),exit(0);
				/* g = g*( (Rx-Sx[delta])*(Rz+Sz[delta])-alpha+beta[delta]) mod n */
				mpz_sub(t,Rx,Sx[delta]); mpz_add(t1,Rz,Sz[delta]);
				mpz_mul(t,t1,t1); mpz_sub(t,t,alpha); mpz_add(t,t,beta[delta]);
				mpz_mul(g,g,t); mpz_mod(g,g,n);
			}
			addh(Rx,Rz,Sx[D-1],Sz[D-1],Tx,Tz,t,t1);
			mpz_set(Tx,Rx); mpz_set(Tz,Rz); mpz_set(Rx,t); mpz_set(Rz,t1);
		}
		mpz_gcd(g,g,n);
		if(mpz_cmp_ui(g,1)>0 && mpz_cmp(g,n)<0) r=1,mpz_set(out,g);
		for(i=0;i<D;i++) mpz_clear(Sx[i]),mpz_clear(Sz[i]),mpz_clear(beta[i]);
		mpz_clear(Tx); mpz_clear(Tz); mpz_clear(Rx); mpz_clear(Rz); mpz_clear(alpha);
		if(r) goto end;
	}
end:
	mpz_clear(ecm.n); mpz_clear(ecm.C);
	mpz_clear(ecm.t); mpz_clear(ecm.t1); mpz_clear(ecm.t2); mpz_clear(ecm.t3);
	mpz_clear(ecm.U); mpz_clear(ecm.V); mpz_clear(ecm.T); mpz_clear(ecm.W);
	mpz_clear(sigma); mpz_clear(u); mpz_clear(v); mpz_clear(Qx); mpz_clear(Qz); mpz_clear(g);
	mpz_clear(t); mpz_clear(t1);
	return r;
}
#undef D

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
		for(q=0;b<=B2;b+=2) if(sieve[b]) {
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
	/* always try pollard p-1 and ecm interleaved */
	mpz_init(m);
	if(pollardp1(n,m,500000,500000,2)) goto found;
	if(ECM(n,m,2000,25)) goto found;
	if(ECM(n,m,11000,90)) goto found;
	if(ECM(n,m,50000,300)) goto found;
	if(pollardp1(n,m,1000000,5000000,2)) goto found;
	if(ECM(n,m,250000,700)) goto found;
	if(ECM(n,m,1000000,1800)) goto found;
	/* TODO strategy:
	   for "small" numbers, run qs immediately (<50 digits?)
	   otherwise, run ECM before QS, the larger the number the more ECM. */
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
	createsieve();
	genprimes();
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
