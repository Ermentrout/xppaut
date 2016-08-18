#include "strutil.h"

#ifndef HAVE_WCTYPE_H
# include <ctype.h>
#else
# include <wctype.h>
#endif
#include <string.h>


/**
 * Removes whitespace from string s.
 *
 * @param s the string clean up.
 */
void de_space(char *s) {
	int j=0;
	char ch;
	for(int i=0; i<strlen(s); i++) {
		ch=s[i];
		if(!isspace(ch)) {
			s[j]=ch;
			j++;
		}
	}
	s[j]=0;
}


/**
 * Copies string s2 into s1 up to position len.
 *
 * @param s1 the string to copy into.
 * @param s2 the string to copy from.
 * @param len the length of the copied string.
 */
void memmov(char *s1, const char *s2, int len) {
	for(int i=0; i<len; i++) {
		s1[i]=s2[i];
	}
}


/**
 * Copies string s2 into s1 backwards from position len.
 *
 * @param s1 the string to copy into.
 * @param s2 the string to copy from.
 * @param len the length of the copied string.
 */
void movmem(char *s1, const char *s2, int len) {
	for(int i=len-1; i>=0; i--) {
		s1[i]=s2[i];
	}
}


/**
 * Cuts string target at the first divergence from string sother.
 *
 * @param target string to compare against.
 * @param sother string to compare with.
 */
void stringintersect(char *target, const char *sother) {
	int j = 0;
	int m = strlen(target);
	int n = strlen(sother);
	if (n < m) {
		m = n;
	}
	while (j < m) {
		if (target[j] != sother[j]) {
			break;
		}
		j++;
	}
	target[j] = '\0';
}


/**
 * Returns whether pre is a prefix of s.
 *
 * @param pre the prefix to search for.
 * @param s the string to search in.
 * @return non-zero if pre is a prefix of s.
 */
int strprefix(const char *pre, const char *s) {
	int n = strlen(pre);

	if (strlen(s) < n) {
		return 0;
	}
	for (int i = 0; i < n; i++) {
		if (pre[i] != s[i]) {
			return 0;
		}
	}
	return 1;
}


/**
 * Turns string s into lowercase.
 *
 * @param s the string to turn lowercase.
 */
void strlwr(char *s) {
	int i=0;
	while(s[i]) {
		if(isupper(s[i])) {
			s[i]+=32;
		}
		i++;
	}
}


/**
 * Turns string s into uppercase.
 *
 * @param s the string to turn uppercase.
 */
void strupr(char *s) {
	int i=0;
	while(s[i]) {
		if(islower(s[i])) {
			s[i]-=32;
		}
		i++;
	}
}
