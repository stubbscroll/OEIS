/* A000020 number of primitive polynomials of degree n over GF(2),
   except for the erroneous a(1)=2.
   formula: phi(2^n-1)/n.
   algorithm: start with N=2^n-1. we have that when d|n, 2^d-1|2^n-1.
   for all d|n, factorize gcd(N,2^d-1) and let N=gcd(N,2^d-1) after each step.
   if n is prime, then all prime factors must be of the form 2kn+1 for
   integer k>=1 and must also be 1 or 7 mod 8.
   the following factorization algorithms are tried in order:
   - trial division
   - pollard-rho
   - pollard p-1
   - elliptic curve method (speed is roughly 66% of alpertron, i think)
   - quadratic sieve (simple variant with one large prime)
   * (maybe also add a check for primes of the form 2kn+1)
   then process the list of prime factors and calculate phi.
   this is very tunable. ecm is much slower than the other methods, so i guess
   pollard-rho and p-1 could be given more time.
   good up to a(292), didn't calculate a(293) because it needs to factorize a
   89 digit composite.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

/* limit for trial division */
#define SMALL 100000000

/* limit (number of digits) for pollard-rho. lowered because QS is fast */
#define POLLARDLIM 30

#define MAXF 10000
mpz_t factors[MAXF];
int fn;

#define ERROR puts("error, too many factors. increase MAXF and recompile"),exit(0)

gmp_randstate_t gmprand;

/* start of factorization routines! */

/* misc routines, such as sieve */

#define MAXP 120000000
typedef long long ll;

char sieve[MAXP];
int prime[6841700];
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

/* start of quadratic sieve! */

void mpz_set_ull(mpz_t m,unsigned long long a) {
	mpz_import(m,1,-1,sizeof(unsigned long long),0,0,&a);
}

/* WARNING, don't call on m > maxull */
/* TODO check for m > maxull */
ll mpz_get_ull(mpz_t m) {
	long long r;
	mpz_export(&r,0,-1,sizeof(long long),0,0,m);
	return r;
}

unsigned int powmod(unsigned int n,unsigned int k,unsigned int mod) {
	unsigned int ans=1;
	while(k) {
		if(k&1) ans=(unsigned long long)ans*n%mod;
		k>>=1;
		n=(unsigned long long)n*n%mod;
	}
	return ans;
}

int rand31() {
	return (rand()&32767)+((rand()&32767)<<15)+((rand()&1)<<30);
}

/* find square root of a modulo p (p prime) */
int sqrtmod(int a,int p) {
	int p8,alpha,i;
	int x,c,s,n,b,J,r2a,r,ret;
	mpz_t za,zp;
	if(p==2) return a&1;
	mpz_init_set_si(za,a);
	mpz_init_set_si(zp,p);
	if(mpz_jacobi(za,zp)!=1) { /* no square root */
		ret=0;
		goto end;
	}
	p8=p&7;
	if(p8==3 || p8==5 || p8==7) {
		if((p8&3)==3) {
			ret=powmod(a,(p+1)/4,p);
			goto end;
		}
		x=powmod(a,(p+3)/8,p);
		c=(ll)x*x%p;
		ret=c==a?x:(ll)x*powmod(2,(p-1)/4,p)%p;
		goto end;
	}
	alpha=0;
	s=p-1;
	while(!(s&1)) s>>=1,alpha++;
	r=powmod(a,(s+1)/2,p);
	r2a=(ll)r*powmod(a,(s+1)/2-1,p)%p;
	do {
		n=rand31()%(p-2)+2;
		mpz_set_si(za,n);
	} while(mpz_jacobi(za,zp)!=-1);
	b=powmod(n,s,p);
	J=0;
	for(i=0;i<alpha-1;i++) {
		c=powmod(b,2*J,p);
		c=(ll)r2a*c%p;
		c=powmod(c,1<<(alpha-i-2),p);
		if(c==p-1) J+=1<<i;
	}
	ret=(ll)r*powmod(b,J,p)%p;
end:
	mpz_clear(zp);
	mpz_clear(za);
	return ret;
}

/* number of extra relations in linear algebra */
#define EXTRAREL 20
/* cache block, experiment with this */
#define BLOCKSIZE 30000

/* structures for partial relations */
struct partial2_s {
	mpz_t rel;               /* smooth number with one large prime */
	unsigned long long *v;   /* exponent vector mod 2 (only factor base primes) */
};

struct partial_s {
	long long p;             /* large prime */
	int p2ix;                /* index into partial array */
};

/* global variables for quadratic sieve */
struct {
	/* basic info */
	mpz_t n;                 /* number to factorize */
	/* factor base */
	int B1;                  /* upper bound for primes in factor base */
	long long B2;            /* upper bound for large prime */
	mpz_t Big2;              /* large prime bound as bigint */
	int LGSLACK;             /* base-2 log slack to accept for trial division */
	int *p;                  /* primes in factor base */
	int *a;                  /* root of n mod p */
	int *lg;                 /* integer logarithm base 2 */
	int fn;                  /* number of primes in factor base */
	/* full relations */
	short sieve[BLOCKSIZE];  /* the sieve */
	mpz_t *rel;              /* smooth numbers, this includes the "first half" of
	                            merged partial relations (actually not the smooth
	                            number itself, but the x resulting in x^2-n) */
	int rn;                  /* number of relations */
	int frn;                 /* number of full relations */
	/* partial relations */
	struct partial_s *part;  /* partial relations, can be sorted, points into
	                            part2 where most of the actual data is */
	struct partial2_s *part2;/* partial relations: smooth number (with large
	                            prime) and exponent vector */
	int prn;                 /* number of partial relations */
	int maxprn;              /* number of partial relations allocated */
	/* merged partial relations consist of 2 almost smooth numbers. the first is
	   in *rel[i], the other is *rel2[i-f] where f is the number of full
	   relations. */
	mpz_t *rel2;             /* smooth numbers, second half of merged partial
	                            relations (actually, the x giving x^2-n) */
	mpz_t *lp;               /* common large prime of merged relation */
	/* hash table for counting occurrences of large primes */
	long long *hkey;         /* keys (large prime) in hash table */
	/* TODO consider making hcount uchar if memory use becomes significant */
	int *hcount;             /* value (count) in hash table */
	long long maxhash;       /* number of elements in hash table */
	int mcnt;                /* number of relations formed from partials */
	/* matrix */
	unsigned long long **m;
	/* final assembly */
	int *ev;                 /* cumulative exponent vector */
	/* statz */
	int false,sieves;        /* failed trial divisions, blocks sieved */
} qs;

/* find index in hash table for element key. if it doesn't exist, create it */
long long QSgethash(long long key) {
	long long ix=key%qs.maxhash;
	while(qs.hkey[ix] && qs.hkey[ix]!=key) {
		ix++;
		if(ix==qs.maxhash) ix=0;
	}
	if(!qs.hkey[ix]) qs.hkey[ix]=key;
	return ix;
}

int QScompp(const void *A,const void *B) {
	const struct partial_s *a=A,*b=B;
	if(a->p<b->p) return -1;
	if(a->p>b->p) return 1;
	return 0;
}

#define SETBITVAL(a,b,v) qs.m[(a)][(b)>>6]=(qs.m[(a)][(b)>>6]&~(1ULL<<((b)&63)))|((v>0)<<((b)&63))
#define SETBIT(a,b) qs.m[(a)][(b)>>6]|=(1ULL<<((b)&63))
#define CLRBIT(a,b) qs.m[(a)][(b)>>6]&=~(1ULL<<((b)&63))
#define XORBIT(a,b) qs.m[(a)][(b)>>6]^=(1ULL<<((b)&63))
#define ISSET(a,b) (qs.m[(a)][(b)>>6]&(1ULL<<((b)&63)))

/* gaussian elimination mod 2 on bitmasks, A is n*m, b is n*o */
/* a is a malloced array of pointers, each a[i] is of size
   sizeof(unsigned long long)*(m+o+63)/64 */
/* TODO optimize this later */
/* return 0: no solutions, 1: one solution, 2: free variables */
int bitgauss64(int n,int m,int o) {
	int i,j,k,z=m+o,c=0,fri=0,bz=(z+63)>>6;
	unsigned long long t;
	/* process each column */
	for(i=0;i<m;i++) {
		/* TODO check words instead of bits */
		for(j=c;j<n;j++) if(ISSET(j,i)) break;
		if(j==n) { fri=1; continue; }
		/* swap? */
		if(j>c)  for(k=0;k<bz;k++) {
			t=qs.m[j][k],qs.m[j][k]=qs.m[c][k],qs.m[c][k]=t;
		}
		/* subtract multiples of this row */
		for(j=0;j<n;j++) if(j!=c && ISSET(j,i)) {
			for(k=0;k<bz;k++) qs.m[j][k]^=qs.m[c][k];
		}
		c++;
	}
	/* detect no solution: rows with 0=b and b!=0 */
	for(i=0;i<n;i++) {
		/* TODO make bit-efficient solution later */
		for(j=0;j<m;j++) if(ISSET(i,j)) goto ok;
		for(;j<z;j++) if(ISSET(i,j)) return 0;
	ok:;
	}
	return 1+fri;
}

int QSgenfactorbase() {
	int i,j;
	mpz_t t;
	mpz_init(t);
	qs.fn=0;
	/* generate factor base! p must be a quadratic residue of n */
	for(qs.fn=2,i=1;i<primes && prime[i]<qs.B1;i++) {
		mpz_set_ui(t,prime[i]);
		if(mpz_jacobi(qs.n,t)>0) qs.fn++;
	}
	/* allocate memory! */
	if(!(qs.p=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(qs.a=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(qs.lg=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(qs.rel=malloc(sizeof(mpz_t)*(qs.fn+EXTRAREL)))) puts("out of memory"),exit(1);
	for(i=0;i<qs.fn+EXTRAREL;i++) mpz_init(qs.rel[i]);
	if(!(qs.m=malloc(sizeof(unsigned long long *)*(qs.fn)))) puts("out of memory"),exit(1);
	for(i=0;i<qs.fn;i++) if(!(qs.m[i]=calloc(((qs.fn+EXTRAREL+63)/64),sizeof(unsigned long long)))) puts("out of memory"),exit(1);
	qs.fn=0;
	qs.p[qs.fn++]=-1;
	qs.p[qs.fn]=2; qs.lg[qs.fn]=1; qs.a[qs.fn++]=1;
	for(i=1;i<primes && prime[i]<qs.B1;i++) {
		mpz_set_ui(t,prime[i]);
		j=mpz_jacobi(qs.n,t);
		/* in the exceedingly rare event that the prime divides n: return factor */
		if(!j) return prime[i];
		if(j>0) {
			qs.p[qs.fn]=prime[i];
			/* find square root a = n (mod p) */
			qs.a[qs.fn]=sqrtmod(mpz_fdiv_ui(qs.n,prime[i]),prime[i]);
			qs.lg[qs.fn]=log(prime[i])/log(2)+.5;
			qs.fn++;
		}
	}
	printf("  factor base bound %d primes %d\n",qs.B1,qs.fn);
	mpz_clear(t);
	return 0;
}

/* dir: 1 if positive x^2-n, -1 if negative */
void QStrialdiv(mpz_t start,int dir) {
	mpz_t t,t1,t2;
	long long big,ix,iy;
	int lim,i,j,lo,hi,mid,k;
	mpz_init_set(t,start); mpz_init(t1); mpz_init(t2);
	/* find lowest value in interval */
	if(dir<0) mpz_add_ui(t,t,BLOCKSIZE-1);
	/* take lg(t^2-n) */
	mpz_mul(t,t,t);
	mpz_sub(t,t,qs.n);
	mpz_abs(t,t);
	lim=0.5+log(mpz_get_d(t))/log(2);
	for(i=0;i<BLOCKSIZE;i++) if(qs.sieve[i]>=lim-qs.LGSLACK) {
		mpz_set(t,start);
		mpz_add_ui(t,t,i);
		mpz_mul(t,t,t);
		mpz_sub(t,t,qs.n);
		mpz_set(t2,t);
		/* put -1 in exponent vector if negative */
		if(mpz_cmp_ui(t,0)<0) SETBIT(0,qs.rn);
		mpz_abs(t,t);
		for(j=1;j<qs.fn;j++) {
			/* break if remainder < next prime^2 */
			if(qs.p[j]<46340) {
				if(mpz_cmp_ui(t,qs.p[j]*qs.p[j])<0) break;
			} else {
				mpz_set_ui(t1,qs.p[j]);
				mpz_mul(t1,t1,t1);
				if(mpz_cmp(t,t1)<0) break;
			}
			while(!mpz_fdiv_ui(t,qs.p[j])) {
				mpz_fdiv_q_ui(t,t,qs.p[j]);
				/* put factor in exponent vector */
				XORBIT(j,qs.rn);
			}
		}
		if(mpz_cmp_ui(t,1)>0 && mpz_cmp_ui(t,qs.p[qs.fn-1])<0) {
			/* find index of remaining factor in list */
			lo=j; hi=qs.fn; k=mpz_get_ui(t);
			while(lo<hi) {
				mid=lo+(hi-lo)/2;
				if(qs.p[mid]<k) lo=mid+1;
				else hi=mid;
			}
			if(k!=qs.p[lo]) printf("sanity error, remainder %d found %d\n",k,qs.p[lo]);
			/* put factor in exponent vector */
			if(lo<0 || lo>=qs.fn) printf("lo %d out of bounds\n",lo);
			XORBIT(lo,qs.rn);
			mpz_set_ui(t,1);
		} else if(mpz_cmp_ui(t,1)>0 && mpz_cmp(t,qs.Big2)<1) {
			/* the remainder is a large prime within our B2 bounds */
			if(qs.prn==qs.maxprn) puts("ERROR, exhausted large prime storage"),exit(1);
			big=mpz_get_ull(t);
			qs.part[qs.prn].p2ix=qs.prn;
			qs.part[qs.prn].p=big;
			mpz_init_set(qs.part2[qs.prn].rel,start);
			mpz_add_ui(qs.part2[qs.prn].rel,start,i);
			if(!(qs.part2[qs.prn].v=calloc((qs.fn+63)/64,sizeof(unsigned long long)))) puts("out of memory"),exit(1);
			for(j=0;j<qs.fn;j++) if(ISSET(j,qs.rn)) qs.part2[qs.prn].v[j>>6]|=(1ULL<<(j&63));
			qs.prn++;
			/* update count */
			if(++qs.hcount[QSgethash(big)]>1) qs.mcnt++;
			/* don't set t to 1 to invoke fail routine on purpose */
		}
		if(!mpz_cmp_ui(t,1)) {
			mpz_set(qs.rel[qs.rn],start);
			mpz_add_ui(qs.rel[qs.rn],qs.rel[qs.rn],i);
			qs.rn++;
		} else {
			/* factorization failed, clear vector */
			for(j=0;j<qs.fn;j++) CLRBIT(j,qs.rn);
			qs.false++;
		}
		/* break if we have enough relations */
		if(qs.rn+qs.mcnt==qs.fn+EXTRAREL) {
			/* if we have enough partial relations, merge them and end sieving */
			if(qs.mcnt) {
				/* sort on largest prime */
				qsort(qs.part,qs.prn,sizeof(qs.part[0]),QScompp);
				/* allocate storage for merged smooth numbers */
				if(!(qs.rel2=malloc(qs.mcnt*sizeof(mpz_t)))) puts("out of memory"),exit(1);
				if(!(qs.lp=malloc(qs.mcnt*sizeof(mpz_t)))) puts("out of memory"),exit(1);
				for(i=0;i<qs.mcnt;i++) mpz_init(qs.rel2[i]),mpz_init(qs.lp[i]);
				i=qs.rn;
				/* for each prime with count>=2: merge 0-1, 0-2, ..., 0-(cnt-1) */
				for(ix=0;ix<qs.prn;) {
					for(iy=ix+1;iy<qs.prn && qs.part[ix].p==qs.part[iy].p;iy++) {
						/* merge ix and iy */
						mpz_set(qs.rel[i],qs.part2[qs.part[ix].p2ix].rel);
						mpz_set(qs.rel2[i-qs.rn],qs.part2[qs.part[iy].p2ix].rel);
						mpz_set_ull(qs.lp[i-qs.rn],qs.part[ix].p);
						for(j=0;j<qs.fn;j++) if(qs.part2[qs.part[ix].p2ix].v[j>>6]&(1ULL<<(j&63)))
							XORBIT(j,i);
						for(j=0;j<qs.fn;j++) if(qs.part2[qs.part[iy].p2ix].v[j>>6]&(1ULL<<(j&63)))
							XORBIT(j,i);
						i++;
					}
					ix=iy;
				}
			}
			/* adjust the meaning of some variables... */
			qs.frn=qs.rn;
			qs.rn+=qs.mcnt;
			break;
		}
	}
	mpz_clear(t); mpz_clear(t1); mpz_clear(t2);
}

void QSsieve() {
	int *pfrontp,*pfrontm,*pbackp,*pbackm;
	int i;
	mpz_t xfront,xback,t;
	if(!(pfrontp=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(pfrontm=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(pbackp=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	if(!(pbackm=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	mpz_init(xfront); mpz_init(xback); mpz_init(t);
	/* start sieving from ceil(sqrt(n)), in both directions */
	mpz_sqrt(xfront,qs.n);
	mpz_add_ui(xfront,xfront,1);
	mpz_set(xback,xfront);
	mpz_sub_ui(xback,xback,BLOCKSIZE);
	qs.rn=0;
	/* find offsets for front and back sieves, both roots */
	for(i=1;i<qs.fn;i++) {
		pfrontp[i]=qs.a[i]-mpz_fdiv_ui(xfront,qs.p[i]);
		while(pfrontp[i]<0) pfrontp[i]+=qs.p[i];
		pfrontm[i]=-qs.a[i]-mpz_fdiv_ui(xfront,qs.p[i]);
		while(pfrontm[i]<0) pfrontm[i]+=qs.p[i];
		pbackp[i]=pfrontp[i]-qs.p[i];
		pbackm[i]=pfrontm[i]-qs.p[i];
		if(pfrontp[i]<0) printf("error");
		if(pfrontm[i]<0) printf("error");
	}
	/* sieve! */
	qs.false=qs.sieves=0;
	printf("  ");
	do {
		/* forward */
		for(i=0;i<BLOCKSIZE;i++) qs.sieve[i]=0;
		i=1;
		if(qs.p[i]==2) {
			while(pfrontp[i]<BLOCKSIZE) qs.sieve[pfrontp[i]]+=qs.lg[i],pfrontp[i]+=qs.p[i];
			pfrontp[i]-=BLOCKSIZE;
			i++;
		}
		for(;i<qs.fn;i++) {
			while(pfrontp[i]<BLOCKSIZE) qs.sieve[pfrontp[i]]+=qs.lg[i],pfrontp[i]+=qs.p[i];
			pfrontp[i]-=BLOCKSIZE;
			while(pfrontm[i]<BLOCKSIZE) qs.sieve[pfrontm[i]]+=qs.lg[i],pfrontm[i]+=qs.p[i];
			pfrontm[i]-=BLOCKSIZE;
		}
		/* find smooth numbers (positive) */
		if(qs.rn+qs.mcnt<qs.fn+EXTRAREL) QStrialdiv(xfront,1);
		mpz_add_ui(xfront,xfront,BLOCKSIZE);
		/* backward */
		if(mpz_cmp_ui(xback,0)<0) continue;
		for(i=0;i<BLOCKSIZE;i++) qs.sieve[i]=0;
		i=1;
		if(qs.p[i]==2) {
			pbackp[i]+=BLOCKSIZE;
			while(pbackp[i]>=0) qs.sieve[pbackp[i]]+=qs.lg[i],pbackp[i]-=qs.p[i];
			i++;
		}
		for(;i<qs.fn;i++) {
			pbackp[i]+=BLOCKSIZE;
			while(pbackp[i]>=0) qs.sieve[pbackp[i]]+=qs.lg[i],pbackp[i]-=qs.p[i];
			pbackm[i]+=BLOCKSIZE;
			while(pbackm[i]>=0) qs.sieve[pbackm[i]]+=qs.lg[i],pbackm[i]-=qs.p[i];
		}
		/* find smooth numbers (negative) */
		if(qs.rn<qs.fn+EXTRAREL) QStrialdiv(xback,-1);
		mpz_sub_ui(xback,xback,BLOCKSIZE);
		qs.sieves++;
		if(qs.sieves%1000000==0) printf("[%d] ",qs.rn);
	} while(qs.rn+qs.mcnt<qs.fn+EXTRAREL);
	printf("%d fr %d pr %d mr %d fail trial division %d sieve blocks\n",qs.rn,qs.prn,qs.rn-qs.frn,qs.false,qs.sieves);
	free(pfrontp); free(pfrontm); free(pbackp); free(pbackm);
	mpz_clear(xfront); mpz_clear(xback); mpz_clear(t);
}

/* build final exponent vector by trial division of x^2-n */
void QSbuildev(mpz_t x) {
	mpz_t t,t1;
	int i,lo,hi,mid,k;
	mpz_init_set(t,x); mpz_init(t1);
	mpz_mul(t,t,t);
	mpz_sub(t,t,qs.n);
	if(mpz_cmp_ui(t,0)<0) {
		qs.ev[0]++;
		mpz_abs(t,t);
	}
	for(i=1;i<qs.fn;i++) {
		if(qs.p[i]<46340) {
			if(mpz_cmp_ui(t,qs.p[i]*qs.p[i])<0) break;
		} else {
			mpz_set_ui(t1,qs.p[i]);
			mpz_mul(t1,t1,t1);
			if(mpz_cmp(t,t1)<0) break;
		}
		while(!mpz_fdiv_ui(t,qs.p[i])) {
			mpz_fdiv_q_ui(t,t,qs.p[i]);
			qs.ev[i]++;
		}
	}
	if(mpz_cmp_ui(t,1)>0 && mpz_cmp_ui(t,qs.p[qs.fn-1])<0) {
		/* find index of remaining factor in list */
		lo=i; hi=qs.fn; k=mpz_get_ui(t);
		while(lo<hi) {
			mid=lo+(hi-lo)/2;
			if(qs.p[mid]<k) lo=mid+1;
			else hi=mid;
		}
		if(k!=qs.p[lo]) printf("sanity error, remainder %d found %d\n",k,qs.p[lo]);
		mpz_set_ui(t,1);
		qs.ev[lo]++;
	}
	if(mpz_cmp(t,qs.Big2)>=0) {
		gmp_printf("sanity error, remaining factor %Zd\n",t);
		exit(1);
	}
	mpz_clear(t); mpz_clear(t1);
}

/* examine each solution vector and try to get factor */
int QSroot(mpz_t a) {
	int i,j,k,r=0,f,tried=0;
	char *freevar,*v;
	mpz_t x,y,t;
	mpz_init(x); mpz_init(y); mpz_init(t);
	if(!(qs.ev=malloc(sizeof(int)*qs.fn))) puts("out of memory"),exit(1);
	/* find all free variables. variable i is free if there is no row having
	   its first 1-element in column i */
	if(!(freevar=malloc(qs.rn))) puts("out of memory"),exit(1);
	if(!(v=malloc(qs.rn))) puts("out of memory"),exit(1);
	for(i=0;i<qs.rn;i++) freevar[i]=1;
	for(i=0;i<qs.fn;i++) {
		for(j=0;j<qs.rn;j++) if(ISSET(i,j)) {
			freevar[j]=0;
			break;
		}
	}
	for(f=0;f<qs.rn;f++) if(freevar[f]) {
		tried++;
		/* set free variable i to 1 and the others to 0 */
		for(i=0;i<qs.rn;i++) v[i]=i==f;
		/* solution vector by back-substitution! set the first 1-element to
		   the xor of the others */
		for(i=qs.fn-1;i>=0;i--) {
			for(j=0;j<qs.rn;j++) if(ISSET(i,j)) goto ok;
			continue;
		ok:
			for(k=j++;j<qs.rn;j++) if(ISSET(i,j) && v[j]) v[k]^=1;
		}
		/* v[i]=1 means that i-th relation is part of the solution */
		/* take square root of left side, the product of x^2 */
		mpz_set_ui(x,1);
		for(i=0;i<qs.rn;i++) if(v[i]) {
			mpz_mul(x,x,qs.rel[i]),mpz_mod(x,x,qs.n);
			if(i>=qs.frn) mpz_mul(x,x,qs.rel2[i-qs.frn]);
		}
		/* take square root of right side, the product of (x^2-n) */
		for(i=0;i<qs.fn;i++) qs.ev[i]=0;
		/* we didn't want to spend lots of memory storing the factorization of
		   each x^2-n, so trial divide again */
		for(i=0;i<qs.rn;i++) if(v[i]) QSbuildev(qs.rel[i]);
		for(i=0;i<qs.rn-qs.frn;i++) if(v[i+qs.frn]) QSbuildev(qs.rel2[i]);
		mpz_set_ui(y,1);
		/* multiply half the exponents */
		for(i=0;i<qs.fn;i++) for(j=0;j<qs.ev[i];j+=2) {
			mpz_mul_si(y,y,qs.p[i]);
			mpz_mod(y,y,qs.n);
		}
		/* multiply in the large primes */
		for(i=0;i<qs.rn-qs.frn;i++) if(v[i+qs.frn]) mpz_mul(y,y,qs.lp[i]),mpz_mod(y,y,qs.n);
		/* get factor */
		mpz_sub(x,x,y);
		mpz_gcd(x,x,qs.n);
		if(mpz_cmp_ui(x,1)>0 && mpz_cmp(x,qs.n)<0) {
			mpz_set(a,x);
			printf("  found factor after %d nullvectors\n",tried);
			r=1;
			break;
		}
	}
	free(v);
	free(freevar);
	free(qs.ev);
	mpz_clear(x); mpz_clear(y); mpz_clear(t);
	return r;
}

/* quadratic sieve! */
/* warning, don't invoke on even integers or powers */
int QS(mpz_t n,mpz_t a) {
	double L=mpz_get_d(n),b1mul,b2mul=200;
	int r,i,d=mpz_sizeinbase(n,10)-50;
	if(d<0) d=0;
	b1mul=0.9-d*0.02;
	if(b1mul<0.4) b1mul=0.4;
	qs.B1=(int)(b1mul*exp(0.5*sqrt(log(L)*log(log(L)))));
	qs.B1+=100;
	qs.B2=(long long)qs.B1*b2mul;
	/* ensure that B2<B1*B1 */
	if(qs.B2>(long long)qs.B1*qs.B1) qs.B2=(long long)qs.B1*qs.B1-1;
	mpz_init(qs.Big2);
	mpz_set_ull(qs.Big2,qs.B2);
	/* the following formula of LGSLACK found by experimentation */
	/* additional slack for large prime */
	qs.LGSLACK=22+d*0.4;
	/* allocate hash table. TODO measure fill degree later */
	qs.maxhash=qs.B2/log(qs.B2)-qs.B1/log(qs.B1);
	if(!(qs.hkey=calloc(qs.maxhash,sizeof(long long)))) puts("out of memory"),exit(1);
	if(!(qs.hcount=calloc(qs.maxhash,sizeof(int)))) puts("out of memory"),exit(1);
	qs.mcnt=0;
	mpz_init_set(qs.n,n);
	if((i=QSgenfactorbase())) {
		/* factor base prime divides n */
		mpz_set_si(a,i);
		r=1;
		goto done;
	}
	qs.maxprn=qs.fn*20+10000; /* max number of partial relations */
	qs.prn=0;
	if(!(qs.part=malloc(qs.maxprn*sizeof(struct partial_s)))) puts("out of memory"),exit(1);
	if(!(qs.part2=malloc(qs.maxprn*sizeof(struct partial2_s)))) puts("out of memory"),exit(1);
	QSsieve();
	bitgauss64(qs.fn,qs.rn,0);
	r=QSroot(a);
	/* dealloc stuff */
	if(qs.mcnt) {
		for(i=0;i<qs.mcnt;i++) mpz_clear(qs.rel2[i]),mpz_clear(qs.lp[i]);
		free(qs.rel2); free(qs.lp);
	}
	for(i=0;i<qs.prn;i++) free(qs.part2[i].v);
	for(i=0;i<qs.prn;i++) mpz_clear(qs.part2[i].rel);
	free(qs.part2);
	free(qs.part);
done:
	free(qs.hcount);
	free(qs.hkey);
	for(i=0;i<qs.fn+EXTRAREL;i++) mpz_clear(qs.rel[i]);
	free(qs.p); free(qs.a); free(qs.lg); free(qs.rel);
	for(i=0;i<qs.fn;i++) free(qs.m[i]);
	free(qs.m);
	mpz_clear(qs.Big2);
	mpz_clear(qs.n);
	return r;
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
	printf("  ECM %d %d\n",B1,maxc);
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
   choose algorithm depending on number of digits */
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
	/* try pollard p-1 and ecm and qs */
	mpz_init(m);
	if(d<=42 && QS(n,m)) goto found;
	if(pollardp1(n,m,500000,500000,2)) goto found;
	if(ECM(n,m,2000,25)) goto found;
	if(d<=50 && QS(n,m)) goto found;
	if(ECM(n,m,11000,90)) goto found;
	if(d<=58 && QS(n,m)) goto found;
	if(ECM(n,m,50000,300)) goto found;
	if(d<=68 && QS(n,m)) goto found;
	if(pollardp1(n,m,1000000,5000000,2)) goto found;
	if(ECM(n,m,250000,700)) goto found;
	if(d<=80 && QS(n,m)) goto found;
	if(ECM(n,m,1000000,1800)) goto found;
	printf("last resort, try QS with %d digits\n",d);
	if(QS(n,m)) goto found;
	goto wrong;
found:
	factorize2(m);
	mpz_divexact(m,n,m);
	factorize2(m);
	mpz_clear(m);
	return;
wrong:
//	gmp_printf("TODO large number (%d digits)\n",d);
	gmp_printf("%Zd d\n",n,d);
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
	for(d=2;d<n;d++) if(n%d==0) {
		/* obtain D=2^d-1 */
		mpz_ui_pow_ui(D,2,d);
		mpz_sub_ui(D,D,1);
		/* take gcd because of overlapping divisors */
		mpz_gcd(D,N,D);
		if(mpz_cmp_ui(D,1)>0) {
			factorize(D);
			mpz_divexact(N,N,D);
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
