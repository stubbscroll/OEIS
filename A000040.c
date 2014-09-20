/* A000040 the prime numbers!
   algorithm: trial division for small n, otherwise deterministic
   miller-rabin */

#include <stdio.h>

#define SMALL 10000

/* check primality by trial division */
int isprime(int n) {
	int i;
	if(n<4) return n>1;
	if(n%2==0) return 0;
	for(i=3;i*i<=n;i++) if(n%i==0) return 0;
	return 1;
}

/* start of miller rabin for unsigned int */

/* evaluate n^k mod mod */
unsigned int powmod(unsigned int n,unsigned int k,unsigned int mod) {
	unsigned int ans=1;
	while(k) {
		if(k&1) ans=(unsigned long long)ans*n%mod;
		k>>=1;
		n=(unsigned long long)n*n%mod;
	}
	return ans;
}

/* subroutine for miller-rabin, only works for odd n */
int witness(unsigned int n,unsigned int a) {
	int s=1,r;
	unsigned int d=(n-1)>>1,x;
	while(!(d&1)) s++,d>>=1;
	x=powmod(a,d,n);
	if(x==1 || x==n-1) return 1;
	for(r=1;r<s;r++) {
		x=(unsigned long long)x*x%n;
		if(x==1) return 0;
		if(x==n-1) return 1;
	}
	return 0;
}

/* deterministic miller-rabin, return 1 if n is prime */
int millerrabin(unsigned int n) {
	if(n<4) return n>1;
	if(!(n&1)) return 0;
	if(n<1373653) return witness(n,2) && witness(n,3);
	if(n<9080191) return witness(n,31) && witness(n,73);
	if(n<25326001) return witness(n,2) && witness(n,3) && witness(n,5);
	return witness(n,2) && witness(n,7) && witness(n,61);
}

/* start of miller rabin for unsigned long long */

/* evaluate a*b%mod. NOT tested for values close to 2^64 */
unsigned long long ullmulmod(unsigned long long a,unsigned long long b,unsigned long long mod) {
	unsigned long long r=0,t;
	if(a>=mod) a%=mod;
	if(b>=mod) b%=mod;
	if(a>b) t=a,a=b,b=t;
	while(b) {
		if(b&1) {
			r+=a;
			if(r>=mod || r<a) r-=mod;
		}
		a=(a<<1)%mod;
		b>>=1;
	}
	return r;
}

unsigned long long ullpowmod(unsigned long long n,unsigned long long k,unsigned long long mod) {
	unsigned long long ans=1;
	while(k) {
		if(k&1) ans=ullmulmod(ans,n,mod);
		k>>=1;
		n=ullmulmod(n,n,mod);
	}
	return ans;
}

/* subroutine for miller-rabin, only works for odd n */
int ullwitness(unsigned long long n,unsigned long long a) {
	unsigned long long d=(n-1)>>1,x;
	int s=1,r;
	while(!(d&1)) s++,d>>=1;
	x=ullpowmod(a,d,n);
	if(x==1 || x==n-1) return 1;
	for(r=1;r<s;r++) {
		x=ullmulmod(x,x,n);
		if(x==1) return 0;
		if(x==n-1) return 1;
	}
	return 0;
}

/* deterministic miller-rabin for 64-bit numbers, return 1 if n is prime */
/* it seems the routine breaks down if n>=2^63 */
int ullmillerrabin(unsigned long long n) {
	if(n<4294967296LU) return millerrabin(n);
	if(!(n&1)) return 0;
	if(n<4759123141LL)
		return ullwitness(n,2) && ullwitness(n,7) && ullwitness(n,61);
	if(n<2152302898747LL)
		return ullwitness(n,2) && ullwitness(n,3) && ullwitness(n,5) &&
		       ullwitness(n,7) && ullwitness(n,11);
	if(n<3474749660383LL)
		return ullwitness(n,2) && ullwitness(n,3) && ullwitness(n,5) &&
		       ullwitness(n,7) && ullwitness(n,11) && ullwitness(n,13);
	if(n<341550071728321LL)
		return ullwitness(n,2) && ullwitness(n,3) && ullwitness(n,5) &&
		       ullwitness(n,7) && ullwitness(n,11) && ullwitness(n,13) &&
		       ullwitness(n,17);
	return ullwitness(n,2) && ullwitness(n,325) && ullwitness(n,9375) &&
	       ullwitness(n,28178) && ullwitness(n,450775) && ullwitness(n,9780504) &&
	       ullwitness(n,1795265022);
}

int main() {
	unsigned long long n,i=1;
	for(n=2;n<SMALL;n++) if(isprime((int)n)) printf("%lld %lld\n",i++,n);
	if(n%2==0) n++;
	for(;n<(1ULL<<32);n+=2) {
		if(millerrabin((unsigned int)n)) printf("%lld %lld\n",i++,n);
	}
	for(;n>10;n+=2) {
		if(ullmillerrabin(n)) printf("%lld %lld\n",i++,n);
	}
	return 0;
}
