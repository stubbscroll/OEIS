/* A000034 sequence of repeating 1, 2 */

#include <stdio.h>

int main() {
	int i;
	for(i=0;i>-1;i++) printf("%d %d\n",i,i%2+1);
	return 0;
}
