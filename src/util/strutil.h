#ifndef XPPAUT_STRUTIL_H
#define XPPAUT_STRUTIL_H

/* --- Functions --- */
void de_space(char *s);
void memmov(char *s1, const char *s2, int len);
void movmem(char *s1, const char *s2, int len);
void stringintersect(char *target, const char *sother);
int strprefix(const char *pre, const char *s);
void strupr(char *s);
void strlwr(char *s);

#endif /* XPPAUT_STRUTIL_H */
