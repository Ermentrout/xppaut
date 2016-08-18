#include "timeutil.h"

#include <stdlib.h>
#include <sys/time.h>
#include <time.h>


/* --- Functions --- */
int gettimenow(void) {
	struct timeval now;
	gettimeofday(&now,NULL);
	return now.tv_usec;
}


void waitasec(int msec) {
	struct timespec rgt;
	rgt.tv_sec = msec / 1000;
	rgt.tv_nsec = (msec % 1000) * 1000000l;
	nanosleep(&rgt, NULL);
}
