/* A000042 unary representation of natural numbers */

#include <stdio.h>

int main() {
	int i,j;
	for(i=1;i>-1;i++) {
		printf("%d ",i);
		for(j=0;j<i;j++) putchar('1');
		putchar('\n');
	}
	return 0;
}
