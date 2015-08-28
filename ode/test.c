#include <stdlib.h>
#include <stdio.h>
#include <string.h>
typedef struct
        {
         char name[10];
         int len;
         int com;
         int arg;
         int pri;
        } SYMBOL;
	

SYMBOL ms[11]=
{  
   {"(",1,999,0,1},      /*  0   */
   {")",1,999,0,2},
   {",",1,999,0,3},
   {"+",1,100,0,4},
   {"-",1,101,0,4},
   {"*",1,102,0,6},
   {"/",1,103,0,6},
   {"^",1,105,0,7},
   {"**",2,105,0,7},
   {"~",1,14,0,6},
   {"START",5,-1,0,0}
};

SYMBOL *m;

double *v;

main()
{
  int i;
  char name[10];
  v=(double *)malloc(11*sizeof(double));
  m=(SYMBOL *)malloc(11*sizeof(SYMBOL));
   for(i=0;i<10;i++){
     strcpy(m[i].name,ms[i].name);
     m[i].len=ms[i].len;
     m[i].com=ms[i].com;
     m[i].arg=ms[i].arg;
     m[i].pri=ms[i].pri;
     v[i]=.05*i;
   }

  for(i=0;i<10;i++)
    printf("%s %d %d %d %d %g \n",m[i].name,m[i].len,m[i].com,m[i].arg,m[i].pri,v[i]);
  /*  m=(SYMBOL *)realloc(m,100*sizeof(SYMBOL));
  v=(double *)realloc(v,100*sizeof(double));  
  for(i=0;i<10;i++)
    printf("%s %d %d %d %d %g\n",m[i].name,m[i].len,m[i].com,m[i].arg,m[i].pri,v[i]);
  */
  for(i=12;i<80;i++){
    sprintf(m[i].name,"var%d",i);
    v[i]=.05*i;
    m[i].len=strlen(m[i].name);
    m[i].com=10000+i;
    m[i].arg=0;
    m[i].pri=17;
  }
  
  printf("here we go \\n");
  for(i=12;i<80;i++)
     printf("%s %d %d %d %d %g\n",m[i].name,m[i].len,m[i].com,m[i].arg,m[i].pri,v[i]);

    
    

}
  
