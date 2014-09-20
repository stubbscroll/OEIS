/* A000038 the sequence that begins with 2 and continues with 0 */

#include <stdio.h>

int main() {
	int i;
	for(i=0;i>-1;i++) printf("%d %d\n",i,(i<1)*2);
	return 0;
}
