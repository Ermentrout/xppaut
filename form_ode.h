
#ifndef _form_ode_h
#define _form_ode_h

#include "xpplim.h"
#include "newpars.h"
#include <stdio.h>

#define MAXVNAM 33
#define MAXLINES 5000


/*void break_up_list(char *rhs);
void compile_em();
void free_varinfo();
void remove_blanks(char *s1);
void read_a_line(FILE *fp,char *s);

void subsk(char *big,char *new,int k,int flag);
void free_comments();

void add_comment(char *s);
void init_varinfo();
void add_varinfo(int type,char *lhs,char *rhs,int nargs,char args[MAXARG][NAMLEN+1]);
void stor_internopts(char *s1);
*/

typedef struct {
  char *name,*value;} FIXINFO;
  
  
  
int make_eqn(void);
void strip_saveqn(void);
int disc(char *string);
void dump_src(void);
void dump_comments(void);
void format_list(char **s, int n);
int get_a_filename(char *filename, char *wild);
void list_em(char *wild);
int read_eqn(void);
int get_eqn(FILE *fptr);
int compiler(char *bob, FILE *fptr);
void list_upar(void);
void welcome(void);
void show_syms(void);
void take_apart(char *bob, double *value, char *name);
char *get_first(char *string, char *src);
char *get_next(char *src);
void find_ker(char *string, int *alt);
void pos_prn(char *s, int x, int y);
void clrscr(void);
int getuch(void);
int getchi(void);
int if_include_file(char *old, char *nf);
int if_end_include(char *old);
int do_new_parser(FILE *fp, char *first, int nnn);
void create_plot_list(void);
void add_only(char *s);
void break_up_list(char *rhs);
int find_the_name(char list[1949][33], int n, char *name);
void compile_em(void);
int formula_or_number(char *expr, double *z);
void strpiece(char *dest, char *src, int i0, int ie);
int parse_a_string(char *s1, VAR_INFO *v);
void init_varinfo(void);
void add_varinfo(int type, char *lhs, char *rhs, int nargs, char args[20][10 +1]);
void free_varinfo(void);
int extract_ode(char *s1, int *ie, int i1);
int strparse(char *s1, char *s2, int i0, int *i1);
int extract_args(char *s1, int i0, int *ie, int *narg, char args[20][10 +1]);
int find_char(char *s1, char *s2, int i0, int *i1);
int next_nonspace(char *s1, int i0, int *i1);
void remove_blanks(char *s1);
void read_a_line(FILE *fp, char *s);
int search_array(char *old, char *new, int *i1, int *i2, int *flag);
int check_if_ic(char *big);
int not_ker(char *s, int i);
int is_comment(char *s);
void subsk(char *big, char *new, int k, int flag);
void keep_orig_comments(void);
void default_comments(void);
void free_comments(void);
void new_comment(FILE *f);
void add_comment(char *s);

/* for parsing par, init with whitespace correctly */
char* new_string2(char* old, int length);
void advance_past_first_word(char** sptr);
char* get_next2(char** tokens_ptr);
void strcpy_trim(char* dest, char* source);
void strncpy_trim(char* dest, char* source, int n);
#endif 
