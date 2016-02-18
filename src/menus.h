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
extern char *fileon_menu[];
extern char *fileoff_menu[];
extern char *main_menu[];
extern char *num_menu[];

extern char *aaxes_hint[];
extern char *adj_hint[];
extern char *afile_hint[];
extern char *arun_hint[];
extern char *aspecial_hint[];
extern char *auto_hint[];
extern char *browse_hint[];
extern char *bvp_hint[];
extern char *cmap_hint[];
extern char *color_hint[];
extern char *edit_hint[];
extern char *edrh_hint[];
extern char *file_hint[];
extern char *flow_hint[];
extern char *frz_hint[];
extern char *graf_hint[];
extern char *half_hint[];
extern char *ic_hint[];
extern char *kin_hint[];
extern char *main_hint[];
extern char *map_hint[];
extern char *meth_hint[];
extern char *no_hint[];
extern char *num_hint[];
extern char *null_hint[];
extern char *null_freeze[];
extern char *phas_hint[];
extern char *sing_hint[];
extern char *tab_hint[];
extern char *text_hint[];
extern char *stoch_hint[];
extern char *view_hint[];
extern char *wind_hint[];

#endif /* XPPAUT_MENUS_H */
