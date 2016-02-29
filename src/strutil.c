#include "strutil.h"

#ifndef HAVE_WCTYPE_H
# include <ctype.h>
#else
# include <wctype.h>
#endif
#include <string.h>

void de_space(char *s) {
	int n=strlen(s);
	int i,j=0;
	char ch;
	for(i=0;i<n;i++) {
		ch=s[i];
		if(!isspace(ch)) {
			s[j]=ch;
			j++;
		}
	}
	s[j]=0;
}


void strupr(char *s) {
	int i=0;
	while(s[i]) {
		if(islower(s[i])) {
			s[i]-=32;
		}
		i++;
	}
}

void strlwr(char *s) {
	int i=0;
	while(s[i]) {
		if(isupper(s[i])) {
			s[i]+=32;
		}
		i++;
	}
}
