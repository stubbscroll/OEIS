/* A067675 number of fixed polyominoes with n cells.
   the algorithm as designed kind of needs to calculate all values in one
   swoop, rather than calculating in-between spitting out terms.
   calculating 1000 terms takes 3 hours 49 minutes on my i7-2600k. */

#include <gmp.h>
#include <stdio.h>
#include <string.h>

/* DP state: [prev/cur][leftmost][rightmost][width][cells used]
   leftmost: 1 if the leftmost part of the current row is to the right
             of the leftmost part of the previous row (that is, we have
             passed the leftmost point)
   rightmost: similar
   width: width of last row
*/

#define MAX 1001
mpz_t dp[2][2][2][MAX][MAX];
int dp2[2][2][2][2][MAX][MAX];

mpz_t table[MAX];

/* precalculate helper values for main algorithm so that we can do innermost
   loop in O(1) with table lookup. precalculation is O(n^3), less than the
   main algorithm */
void precalc() {
	int i,j,k,m,left,right,forceleft,forceright;
	memset(dp2,0,sizeof(dp2));
	for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=1;k<MAX;k++) for(m=1;m<MAX;m++) {
		/* given state i,j,k,m where i,j are current values of leftmost and
		   rightmost, k is the width of the last row and m is the width of the
		   current row, try all ways of placing the current row and count the
		   number of times each of the 4 cases occur (new values of leftmost/
		   rightmost) */
		if(i+j>1 && m>k) break;
		for(left=1-m;left<k;left++) {
			if(left<0 && i) continue;
			right=left+m;
			if(right>k && j) continue;
			forceleft=left>0?1:i;
			forceright=right<k?1:j;
			dp2[forceleft][forceright][i][j][k][m]++;
		}
	}
}

/* calculate a(1) ... a(max) with DP
   time complexity: O(n^4):
   - do/while loop is executed n times
   - loop over n^2 states containing another n loop that
     generates states for the next iteration
   the bound is not necessarily tight, many states never occur.
   (can probably be reduced to O(n^3.5) with some symmetry and making sure
    that we only directly generate polyominoes with height<=width and do some
    very careful counting.)
   and then there's the complexity from the size of the terms, which i haven't
   thought about yet. i guess it depends on things like the complexity of
   gmp's multiply function among other things.
   disclaimer, there might be much better approaches.
*/
int calcterms() {
	int cur=1,prev=0,done,i,j,k,l,m;
	precalc();
  for(i=1;i<MAX;i++) mpz_init_set_si(table[i],0);
	/* initialize dp array */
	for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<2;k++)
		for(l=0;l<MAX;l++) for(m=0;m<MAX;m++)
			mpz_init_set_si(dp[i][j][k][l][m],0);
	/* initial dp with a top row of every possible length */
	for(i=1;i<=MAX;i++) mpz_init_set_si(dp[prev][0][0][i][i],1);
	do {
		done=1;
		/* initialize cur */
		for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<MAX;k++)
			for(l=0;l<MAX;l++) mpz_set_si(dp[cur][i][j][k][l],0);
		for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=1;k<MAX;k++)
			for(l=1;l<MAX;l++) if(mpz_cmp_si(dp[prev][i][j][k][l],0)) {
				done=0;
				/* finish off the current polyomino and put it in table */
				mpz_add(table[l],table[l],dp[prev][i][j][k][l]);
				/* try all widths for next row */
				for(m=1;m<=MAX-l;m++) {
					if(i+j>1 && m>k) break;
					if(dp2[0][0][i][j][k][m]) mpz_addmul_ui(dp[cur][0][0][m][l+m],dp[prev][i][j][k][l],dp2[0][0][i][j][k][m]);
					if(dp2[0][1][i][j][k][m]) mpz_addmul_ui(dp[cur][0][1][m][l+m],dp[prev][i][j][k][l],dp2[0][1][i][j][k][m]);
					if(dp2[1][0][i][j][k][m]) mpz_addmul_ui(dp[cur][1][0][m][l+m],dp[prev][i][j][k][l],dp2[1][0][i][j][k][m]);
					if(dp2[1][1][i][j][k][m]) mpz_addmul_ui(dp[cur][1][1][m][l+m],dp[prev][i][j][k][l],dp2[1][1][i][j][k][m]);
				}
			}
		cur=prev;
		prev=1-cur;
	} while(!done);
	/* clean up */
	for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<2;k++)
		for(l=0;l<MAX;l++) for(m=0;m<MAX;m++) mpz_clear(dp[i][j][k][l][m]);
	return 1;
}

int main() {
	int i;
	calcterms();
	for(i=1;i<MAX;i++) gmp_printf("%d %Zd\n",i,table[i]);
	return 0;
}
