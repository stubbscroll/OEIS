/* A064297 triangle of selv-avoiding rook paths joining opposite corners of
   n*k board.
   for an explanation of the algorithm, see my implementation of A007764 at
   https://github.com/stubbscroll/OEIS/blob/master/A007764-fast.c
   this implementation is very similar and is quite optimized, but does not
   use parallelism.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

typedef unsigned long long ull;

void mpz_set_ull(mpz_t b,ull a) {
  mpz_import(b,1,1, sizeof(a),0,0,&a);
}

/* define the number of bits! 64-bits compared to 16-bits is expected
   to be 4 times faster, using 4 times as much memory */
#define BITS64

#define MAX 32
ull dp[2][MAX][MAX];

/* calculate the number of states (and partial states) for rank/unrank */
/* the sequence dp[1][0][n] is actually A002026 */
void init() {
	int i,j,k;
	for(k=0;k<2;k++) for(i=0;i<MAX;i++) for(j=0;j<MAX;j++) dp[k][i][j]=0;
	dp[0][0][0]=1;
	for(i=0;i<MAX-1;i++) for(k=0;k<2;k++) for(j=0;j<MAX-1;j++) if(dp[k][j][i]) {
		dp[k][j][i+1]+=dp[k][j][i];
		if(!k && !j) dp[1][j][i+1]+=dp[k][j][i];
		if(j) dp[k][j-1][i+1]+=dp[k][j][i];
		if(j<MAX-1) dp[k][j+1][i+1]+=dp[k][j][i];
	}
}

/* convert from integer rank to representation in linear time */
ull unrank(int i,ull r) {
	int j=0;
	ull c0,mask=0;
	while(i--) {
		c0=dp[1][j][i];
		if(r<c0) mask<<=2;
		else {
			r-=c0;
			c0=(!j)?dp[0][0][i]:0;
			if(r<c0) { mask=(mask<<2)+1; goto seen1; }
			else {
				r-=c0;
				c0=(j)?dp[1][j-1][i]:0;
				if(r<c0) mask=(mask<<2)+2,j--;
				else r-=c0,mask=(mask<<2)+3,j++;
			}
		}
	}
	return mask;
seen1:
	while(i--) {
		c0=dp[0][j][i];
		if(r<c0) mask<<=2;
		else {
			r-=c0;
			c0=(j)?dp[0][j-1][i]:0;
			if(r<c0) mask=(mask<<2)+2,j--;
			else r-=c0,mask=(mask<<2)+3,j++;
		}
	}
	return mask;
}

/* convert from represention to integer rank in linear time */
ull rank(int i,ull mask) {
	int j=0,cur;
	ull r=0;
	while(i--) {
		cur=(mask>>(i<<1))&3;
		if(cur==2) {
			r+=dp[1][j][i];
			j--;
		} else if(cur==3) {
			r+=dp[1][j][i]+(j?dp[1][j-1][i]:dp[0][0][i]);
			j++;
		} else if(cur) { r+=dp[1][j][i]; goto seen1; }
	}
	return r;
seen1:
	while(i--) {
		cur=(mask>>(i<<1))&3;
		if(cur==2) {
			r+=dp[0][j][i];
			j--;
		} else if(cur==3) {
			r+=dp[0][j][i];
			if(j) r+=dp[0][j-1][i];
			j++;
		}
	}
	return r;
}

#ifdef BITS16
	typedef unsigned short u16;
	static u16 maxprime=32767;
#elif defined(BITS64)
	typedef ull u16;
	static u16 maxprime=9223372036854775783;
#endif

u16 *prev,*cur;

#define ADD(ix,c,mod) { cur[ix]+=c; if(cur[ix]>=mod) cur[ix]-=mod; }

int swap[16]={0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15};

u16 calc(int n,int m,u16 mod) {
	ull num=dp[1][0][n+1],z,mask,nz,newmask;
	int i,j,left,up,k,w=n+1,l,o,look;
	u16 c,r=0,*t;
	if(!(prev=malloc(sizeof(u16)*num))) { puts("out of memory"); exit(1); }
	if(!(cur=malloc(sizeof(u16)*num))) { puts("out of memory"); exit(1); }
	memset(prev,0,sizeof(u16)*num);
	memset(cur,0,sizeof(u16)*num);
	/* force first edge to go down. the rest of the paths can be obtained by flipping
	   along the diagonal. hence, count the first case only and multiply by 2 */
	prev[rank(w,1)]=1;
	if(m!=n) prev[rank(w,4)]=1;	/* cannot use symmetry if non-square */
	for(i=0;i<m;i++) for(j=0;j<n;j++) {
		if(i==0 && j==0) continue;
		else if(i<m-1 && j<n-1) {
			/* regular cell */
			for(z=0;z<num;z++) if((c=prev[z])) {
				mask=unrank(w,z);
				left=(mask>>(j<<1))&3;
				up=(mask>>((j<<1)+2))&3;
				if(left==3 && up==2) {
					/* join, easy case: 32 => 00 */
					nz=rank(w,mask&(~(15ULL<<(j<<1))));
					ADD(nz,c,mod);
				} else if(left==2 && up==2) {
					/* join 22: find mate of right 2, change it from 3 to 2 */
					for(k=j+2,l=1;;k++) {
						o=(mask>>(k<<1))&3;
						if(o==2) l++;
						else if(o==3) {
							l--;
							if(!l) break;
						}
					}
					newmask=mask&(~(15ULL<<(j<<1)));
					newmask=newmask^(1ULL<<(k<<1));
					nz=rank(w,newmask);
					ADD(nz,c,mod);
				} else if(left==3 && up==3) {
					/* join 33: find mate of left 3, change it from 2 to 3 */
					for(k=j-1,l=1;;k--) {
						o=(mask>>(k<<1))&3;
						if(o==3) l++;
						else if(o==2) {
							l--;
							if(!l) break;
						}
					}
					newmask=mask&(~(15ULL<<(j<<1)));
					newmask=newmask^(1ULL<<(k<<1));
					nz=rank(w,newmask);
					ADD(nz,c,mod);
				} else if(left==1 && up) {
					/* we have 12, find up's mate, change it from 3 to 1 */
					for(k=j+2,l=1;;k++) {
						o=(mask>>(k<<1))&3;
						if(o==2) l++;
						else if(o==3) {
							l--;
							if(!l) break;
						}
					}
					nz=rank(w,mask&(~(15ULL<<(j<<1)))&(~(2ULL<<(k<<1))));
					ADD(nz,c,mod);
				} else if(left && up==1) {
					/* we have 31, find left's mate, change it from 2 to 1 */
					for(k=j-1,l=1;;k--) {
						o=(mask>>(k<<1))&3;
						if(o==3) l++;
						else if(o==2) {
							l--;
							if(!l) break;
						}
					}
					newmask=mask&(~(15ULL<<(j<<1)));
					newmask=newmask^(3ULL<<(k<<1));
					nz=rank(w,newmask);
					ADD(nz,c,mod);
				} else if(left) { /* extend single edge, case 1 */
					/* extend down: no change in mask */
					ADD(z,c,mod);
					/* extend to the right: x0 becomes 0x, but only if next
					   cell isn't 2 */
					look=(mask>>((j<<1)+4))&3;
					if(left!=2 || look!=3) {
						nz=rank(w,(mask&(~(15ULL<<(j<<1))))|((ull)swap[left]<<(j<<1)));
						ADD(nz,c,mod);
					}
				} else if(up) { /* extend single edge, case 2 */
					/* extend down: 0x becomes x0 */
					nz=rank(w,(mask&(~(15ULL<<(j<<1))))|((ull)swap[up<<2]<<(j<<1)));
					ADD(nz,c,mod);
					/* extend to the right: no change in mask */
					look=(mask>>((j<<1)+4))&3;
					if(up!=2 || look!=3) ADD(z,c,mod);
				} else if(!(up|left)) { /* no edge */
					/* place nothing */
					ADD(z,c,mod);
					/* place 23 */
					nz=rank(w,mask|(14ULL<<(j<<1)));
					ADD(nz,c,mod);
				} else printf("error uncatched regular %d %d\n",left,up);
			}
		} else if(i<m-1 && j==n-1) {
			/* right column: edges cannot go to the right, mask<<2 */
			for(z=0;z<num;z++) if((c=prev[z])) {
				mask=unrank(w,z);
				left=(mask>>(j<<1))&3;
				up=(mask>>((j<<1)+2))&3;
				if(left==3 && up==3) {
					/* join 33: find mate of left 3, change it from 2 to 3 */
					for(k=j-1,l=1;;k--) {
						o=(mask>>(k<<1))&3;
						if(o==3) l++;
						else if(o==2) {
							l--;
							if(!l) break;
						}
					}
					newmask=mask&(~(15ULL<<(j<<1)));
					newmask=newmask^(1ULL<<(k<<1));
					nz=rank(w,newmask<<2);
					ADD(nz,c,mod);
				} else if(left==1 && up) {
					/* we have 12, find up's mate, change it from 3 to 1 */
					for(k=j+2,l=1;;k++) {
						o=(mask>>(k<<1))&3;
						if(o==2) l++;
						else if(o==3) {
							l--;
							if(!l) break;
						}
					}
					nz=rank(w,(mask&(~(15ULL<<(j<<1)))&(~(2ULL<<(k<<1))))<<2);
					ADD(nz,c,mod);
				} else if(left && up==1) {
					/* we have 31, find left's mate, change it from 2 to 1 */
					for(k=j-1,l=1;;k--) {
						o=(mask>>(k<<1))&3;
						if(o==3) l++;
						else if(o==2) {
							l--;
							if(!l) break;
						}
					}
					newmask=mask&(~(15ULL<<(j<<1)));
					newmask=newmask^(3ULL<<(k<<1));
					nz=rank(w,newmask<<2);
					ADD(nz,c,mod);
				} else if(left) { /* extend single edge, case 1 */
					/* extend down */
					newmask=mask<<2;
					nz=rank(w,newmask);
					ADD(nz,c,mod);
				} else if(up) { /* extend single edge, case 2 */
					/* extend down: 0x becomes x0 */
					nz=rank(w,((mask&(~(15ULL<<(j<<1))))|((ull)swap[up<<2]<<(j<<1)))<<2);
					ADD(nz,c,mod);
				} else if(!(up|left)) { /* no edge */
					/* place nothing */
					newmask=mask<<2;
					nz=rank(w,newmask);
					ADD(nz,c,mod);
				} else printf("error uncatched regular %d %d\n",left,up);
			}
		} else if(i==m-1 && j<n-1) {
			/* lower row_ edges cannot go down */
			for(z=0;z<num;z++) if((c=prev[z])) {
				mask=unrank(w,z);
				left=(mask>>(j<<1))&3;
				up=(mask>>((j<<1)+2))&3;
				if(left==3 && up==2) {
					/* join, easy case: 32 => 00 */
					nz=rank(w,mask&(~(15ULL<<(j<<1))));
					ADD(nz,c,mod);
				} else if(left==2 && up==2) {
					/* join 22: find mate of right 2, change it from 3 to 2 */
					for(k=j+2,l=1;;k++) {
						o=(mask>>(k<<1))&3;
						if(o==2) l++;
						else if(o==3) {
							l--;
							if(!l) break;
						}
					}
					newmask=mask&(~(15ULL<<(j<<1)));
					newmask=newmask^(1ULL<<(k<<1));
					nz=rank(w,newmask);
					ADD(nz,c,mod);
				} else if(left==3 && up==3) {
					/* join 33: find mate of left 3, change it from 2 to 3 */
					for(k=j-1,l=1;;k--) {
						o=(mask>>(k<<1))&3;
						if(o==3) l++;
						else if(o==2) {
							l--;
							if(!l) break;
						}
					}
					newmask=mask&(~(15ULL<<(j<<1)));
					newmask=newmask^(1ULL<<(k<<1));
					nz=rank(w,newmask);
					ADD(nz,c,mod);
				} else if(left==1 && up) {
					/* we have 12, find up's mate, change it from 3 to 1 */
					for(k=j+2,l=1;;k++) {
						o=(mask>>(k<<1))&3;
						if(o==2) l++;
						else if(o==3) {
							l--;
							if(!l) break;
						}
					}
					nz=rank(w,mask&(~(15ULL<<(j<<1)))&(~(2ULL<<(k<<1))));
					ADD(nz,c,mod);
				} else if(left && up==1) {
					/* we have 31, find left's mate, change it from 2 to 1 */
					for(k=j-1,l=1;;k--) {
						o=(mask>>(k<<1))&3;
						if(o==3) l++;
						else if(o==2) {
							l--;
							if(!l) break;
						}
					}
					newmask=mask&(~(15ULL<<(j<<1)));
					newmask=newmask^(3ULL<<(k<<1));
					nz=rank(w,newmask);
					ADD(nz,c,mod);
				} else if(left) { /* extend single edge, case 1 */
					/* extend to the right: x0 becomes 0x, but only if next
					   cell isn't 2 */
					look=(mask>>((j<<1)+4))&3;
					if(left!=2 || look!=3) {
						nz=rank(w,(mask&(~(15ULL<<(j<<1))))|((ull)swap[left]<<(j<<1)));
						ADD(nz,c,mod);
					}
				} else if(up) { /* extend single edge, case 2 */
					/* extend to the right: no change in mask */
					look=(mask>>((j<<1)+4))&3;
					if(up!=2 || look!=3) ADD(z,c,mod);
				} else if(!(up|left)) { /* no edge */
					/* place nothing */
					ADD(z,c,mod);
				} else printf("error uncatched regular %d %d\n",left,up);
			}
		} else if(i==m-1 && j==n-1) {
			/* lower right corner, just take the sum */
			for(z=0;z<num;z++) if((c=prev[z])) {
				r+=c;
				if(r>=mod) r-=mod;
			}
		}
		t=prev; prev=cur; cur=t;
		memset(cur,0,sizeof(u16)*num);
	}
	free(cur);
	free(prev);
	if(n==m) r=(r*2)%mod;
	return r;
}
#undef ADD

/* solve a set of modular equations using chinese remainder theorem */
/* n: number of equations of the form
   x = a[i] mod b[i] */
void chinese(int n,ull *a,ull *b,mpz_t x) {
	int i;
	mpz_t N,temp,t2,n2,bi;
	mpz_init(temp);
	mpz_init(n2);
	mpz_init(bi);
	mpz_init(t2);
	mpz_set_si(x,0);
	mpz_init_set_si(N,1);
	for(i=0;i<n;i++) {
		mpz_set_ull(temp,b[i]);
		mpz_mul(N,N,temp);
	}
	for(i=0;i<n;i++) {
		mpz_set_ull(temp,a[i]);
		mpz_set_ull(bi,b[i]);
		mpz_tdiv_q(n2,N,bi);
		mpz_mul(temp,temp,n2);
		mpz_tdiv_r(t2,n2,bi);
		mpz_invert(t2,t2,bi);
		mpz_mul(temp,temp,t2);
		mpz_add(x,x,temp);
		mpz_tdiv_r(x,x,N);
	}
	mpz_clear(N);
	mpz_clear(temp);
	mpz_clear(n2);
	mpz_clear(bi);
	mpz_clear(t2);
}

/* m,n: dimensions of rectangle with n*m internal cells */
/* when solving A007764, invoke with n+1,n+1 */
void solve(int n,int m,mpz_t ans) {
	ull prime=maxprime;
	int o=0,t,tries;
	static ull res[500],prim[500];
	mpz_t mul,temp;
	/* upper bound proved by Bousquet-Melou et al */
	/* TODO find a better upper bound, it will improve run time by a lot */
	double bound=pow(1.782,(n+1)*(m+1));
	if(n>m) t=n,n=m,m=t;
	if(n==1) { mpz_set_si(ans,1); return; }
	if(n>=MAX-2) puts("grid is too large"),exit(1);
	mpz_init_set_si(mul,1);
	mpz_init(temp);
	/* require n<=m, frontier width is n+1 */
	do {
		while(1) {
			if(n<10) tries=30;
			else tries=200;
			mpz_set_ull(temp,prime);
			/* we want to be REALLY sure */
			if(mpz_probab_prime_p(temp,tries)) break;
			prime-=2;
		}
		res[o]=calc(n,m,prime);
		prim[o++]=prime;
		mpz_mul(mul,mul,temp);
		/* if mul > upper bound, then break */
		prime-=2;
	} while(mpz_cmp_d(mul,bound)<0);
	/* glue together result with chinese remainder theorem */
	chinese(o,res,prim,mul);
	mpz_set(ans,mul);
	mpz_clear(mul);
	mpz_clear(temp);
}

int main() {
	int height,width,i=1;
	mpz_t r;
	mpz_init(r);
	init();
	for(height=0;;height++) for(width=0;width<=height;width++) {
		solve(height+1,width+1,r);
		gmp_printf("%d %Zd\n",i++,r);
	}
	mpz_clear(r);
	return 0;
}
