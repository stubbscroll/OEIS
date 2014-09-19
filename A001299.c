/* A001299 number of ways for coin change problem with 1,5,10,25.
   algorithm: the standard dp algorithm.
   calculates up to a(n)<=maxlonglong (n<=41050344) */

#include <stdio.h>
#include <stdlib.h>

#define MAX 41050345
long long dp[MAX];

int val[]={5,10,25,0};

int main() {
	int i,j;
	for(i=0;i<MAX;i++) dp[i]=1;
	for(j=0;val[j];j++) for(i=val[j];i<MAX;i++) dp[i]+=dp[i-val[j]];
	for(i=0;i<MAX;i++) printf("%d %lld\n",i,dp[i]);
	return 0;
}
