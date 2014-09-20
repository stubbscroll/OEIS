/* A000037 non-squares */

#include <stdio.h>

int main() {
	long long next=1,i,a;
	for(i=a=1;i>-1;i++,a++) {
		if(a==next*next) a++,next++;
		printf("%lld %lld\n",i,a);
	}
	return 0;
}
