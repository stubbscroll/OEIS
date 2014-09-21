/* A000010 euler's totient (phi) function
   calculate all values using backtracking in time
   O(total number of (non-distinct) prime factors in all numbers)
   timing for 10 million numbers including output to file to SSD:
   regular version: 17.140 s
   this version:     3.180 s
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

typedef long long ll;

#define MAX 100000000
int phivalues[MAX];
char sieve[MAX];
int prime[5761455]; /* number of primes less than MAX */
int primes;

void createsieve() {
	int i,j,q;
	memset(sieve,1,sizeof(sieve));
	q=sqrt(MAX);
	for(sieve[0]=sieve[1]=0,i=2;i<=q;i++) if(sieve[i])
		for(j=i*i;j<MAX;j+=i) sieve[j]=0;
}

void genprimes() {
	int i;
	for(primes=i=0;i<MAX;i++) if(sieve[i]) prime[primes++]=i;
}

/* generate all phi up to MAX, needs all prime[]<MAX.
   this routine is faster than calculating each phi(i) individually */
void phigenbtr(int n,int t,int p) {
	int n2,t2;
	phivalues[n]=t;
	while(p<primes && (ll)n*prime[p]<MAX) {
		n2=n*prime[p];
		t2=t*(prime[p]-1);
		phigenbtr(n2,t2,p+1);
		while((ll)n2*prime[p]<MAX) {
			n2*=prime[p];
			t2*=prime[p];
			phigenbtr(n2,t2,p+1);
		}
		p++;
	}
}

void phigen() {
	phigenbtr(1,1,0);
}

int main() {
	int i;
	createsieve();
	genprimes();
	phigen();
	for(i=1;i<MAX;i++) printf("%d %d\n",i,phivalues[i]);
	return 0;
}
