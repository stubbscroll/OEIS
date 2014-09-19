/* A000008 number of ways for coin change problem with 1,2,5,10.
   algorithm: the standard dp algorithm.
   calculates up to a(n)<=maxlonglong (n<=17688056) */

#include <stdio.h>
#include <stdlib.h>

#define MAX 17688057
long long dp[MAX];

int val[]={2,5,10,0};

int main() {
	int i,j;
	for(i=0;i<MAX;i++) dp[i]=1;
	for(j=0;val[j];j++) for(i=val[j];i<MAX;i++) dp[i]+=dp[i-val[j]];
	for(i=0;i<MAX;i++) printf("%d %lld\n",i,dp[i]);
	return 0;
}
