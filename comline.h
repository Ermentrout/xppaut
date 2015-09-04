#ifndef _comline_h_
#define _comline_h_


typedef struct {
  char *name;
  char *does;
  unsigned int use;
} INTERN_SET;


typedef struct {
   char *name;
   struct SET_NAME * next;
} SET_NAME;



int is_set_name(SET_NAME *set, char *nam);
SET_NAME *add_set(SET_NAME *set, char *nam);
SET_NAME *rm_set(SET_NAME *set, char *nam);
void do_comline(int argc, char **argv);
int if_needed_select_sets(void);
int if_needed_load_set(void);
int if_needed_load_par(void);
int if_needed_load_ic(void);
int if_needed_load_ext_options(void);
int parse_it(char *com);


#endif
