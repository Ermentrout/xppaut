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
void free_finfo(FILEINFO *ff);
int cmpstringp(const void *p1, const void *p2);
int get_fileinfo_tab(char *wild, char *direct, FILEINFO *ff,char *wild2);
int get_fileinfo(char *wild, char *direct, FILEINFO *ff);
int fil_count(char *direct, int *ndir, int *nfil, char *wild, int *mld, int *mlf);
int change_directory(char *path);
int get_directory(char *direct);
int IsDirectory(char *root, char *path);
void MakeFullPath(char *root, char *filename, char *pathname);
int wild_match(char *string, char *pattern);

#endif /* XPPAUT_READ_DIR_H */
