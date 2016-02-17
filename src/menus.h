#ifndef XPPAUT_MENUS_H
#define XPPAUT_MENUS_H

#include <X11/Xlib.h>

/* --- Macros --- */
#define MAIN_MENU 0
#define FILE_MENU 1
#define NUM_MENU 2
#define MAIN_ENTRIES 20
#define FILE_ENTRIES 16
#define NUM_ENTRIES 18

/* --- Types --- */
typedef struct {
  Window base;
  Window title;
  Window w[25];
  char key[25];
  char **names;
  char **hints;
  int n;
  int visible;
} MENUDEF;

/* --- Data --- */
extern char *main_menu[];
extern char *num_menu[];
extern char *fileon_menu[];
extern char *fileoff_menu[];
extern char *main_hint[];
extern char *file_hint[];
extern char *num_hint[];
extern char *auto_hint[];
extern char *no_hint[];
extern char *aaxes_hint[];
extern char *afile_hint[];
extern char *arun_hint[];
extern char *aspecial_hint[];
extern char *edrh_hint[];

#endif /* XPPAUT_MENUS_H */
