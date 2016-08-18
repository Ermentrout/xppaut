#ifndef XPPAUT_READ_DIR_H
#define XPPAUT_READ_DIR_H

#include "load_eqn.h"

/* --- Types --- */
typedef struct {
	char **dirnames,**filenames;
	int nfiles,ndirs;
} FILEINFO;

/* --- Data --- */
extern char cur_dir[XPP_MAX_NAME];
extern FILEINFO my_ff;

/* --- Functions --- */
int change_directory(const char *path);
void free_finfo(FILEINFO *ff);
int get_fileinfo_tab(const char *wild, const char *direct, FILEINFO *ff,const char *wild2);
int get_fileinfo(const char *wild, const char *direct, FILEINFO *ff);
int get_directory(char *direct);

#endif /* XPPAUT_READ_DIR_H */
