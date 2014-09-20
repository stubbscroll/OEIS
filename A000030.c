/* A000030 initial digit of n */

#include <stdio.h>

int main() {
	int i,j;
	for(i=0;i>-1;i++) {
		j=i;
		while(j>9) j/=10;
		printf("%d %d\n",i,j);
	}
	return 0;
}
