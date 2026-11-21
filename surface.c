/*  Rotatable surface view of the array plot.

    The array plot shows z = data[index][time] as colored pixels.  This
    window draws the same grid as a 3D surface: x = index, y = time,
    height = value.  Orthographic projection, back-to-front (painter's)
    sorting of the quads, drawn into a pixmap and copied to the window.

    Mouse:  drag with button 1 to rotate, wheel (buttons 4/5) to zoom.
    Keys:   + / -   zoom        [ / ]   shrink / stretch the height
            m       mesh on/off r       reset the view
    Buttons GIF / PS / SVG save the picture (PS is EPS, both are vector).
    The picture is first built as a list of primitives (polygons, lines,
    text) in window pixel coordinates; the screen, PostScript, SVG and GIF
    all draw from that one list, so a saved figure matches the screen.
    Only core Xlib is used (works with X11 and XQuartz).
*/
#include "surface.h"
#include "stdfix.h"
#include "arrayplot.h"
#include "browse.h"
#include "pop_list.h"
#include "color.h"
#include "xpplim.h"
#include "load_eqn.h"
#include "ggets.h"
#include "init_conds.h"
#include "scrngif.h"
#include "lunch-new.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FIRSTCOLOR 30
#define MESH_AUTO_MAX 1500 /* mesh starts off for grids bigger than this many quads */
#define FAST_QUADS 12000   /* above this, decimate the grid while dragging */

extern Display *display;
extern int screen,DCURXs,DCURYs,CURY_OFFs,color_total,COLOR;
extern unsigned int GrFore,GrBack;
extern GC small_gc;
extern BROWSER my_browser;
extern APLOT aplot;

extern int NJMP;
extern double DELTA_T,TRANS;
extern char this_file[XPP_MAX_NAME];

static int s_alive=0;
static Window s_base,s_win,s_bredraw,s_bmesh,s_breset,s_bclose,s_bgif,s_bps,s_bsvg;
static Pixmap s_pix=0;
static GC s_gc;
static int s_w=520,s_h=520;      /* base window size   */
static int s_pw,s_ph,s_top;      /* drawing area size, y offset of drawing area */
static double s_az=30.0,s_el=35.0,s_zoom=1.0,s_zs=1.0;
static int s_mesh=1;
static int s_drag=0,s_lx,s_ly;

/* trig of the current view, set by set_view() */
static double ca,sa,ce,se,SC;
static double CX,CY;

/* ------------------------------------------------------------------ */
/*  The scene: a list of primitives in window pixel coordinates        */
/* ------------------------------------------------------------------ */
#define C_BACK (-1)   /* background (white in saved files) */
#define C_FORE (-2)   /* foreground (black in saved files) */
enum { P_POLY, P_LINES, P_TEXT };
#define A_LEFT 0
#define A_RIGHT 1
#define A_CENTER 2

typedef struct {
  int kind,n,colr,stroke,anchor;
  double x[5],y[5];
  char t[128];
} PRIM;

typedef struct { PRIM *p; int n,cap; } SCENE;

static PRIM *new_prim(SCENE *sc,int kind)
{
  PRIM *p;
  if(sc->n>=sc->cap){
    int nc=sc->cap?2*sc->cap:1024;
    PRIM *np=(PRIM *)realloc(sc->p,sizeof(PRIM)*nc);
    if(np==NULL)return NULL;
    sc->p=np;sc->cap=nc;
  }
  p=&sc->p[sc->n++];
  memset(p,0,sizeof(PRIM));
  p->kind=kind;
  return p;
}

static void add_line(SCENE *sc,double x1,double y1,double x2,double y2)
{
  PRIM *p=new_prim(sc,P_LINES);
  if(p==NULL)return;
  p->n=2;p->x[0]=x1;p->y[0]=y1;p->x[1]=x2;p->y[1]=y2;
}

static void add_text(SCENE *sc,double x,double y,int anchor,const char *s)
{
  PRIM *p=new_prim(sc,P_TEXT);
  if(p==NULL)return;
  p->anchor=anchor;
  snprintf(p->t,sizeof(p->t),"%s",s);
  {  /* keep the text inside the picture (width estimated from the screen font) */
    double w=(double)strlen(p->t)*DCURXs;
    if(anchor==A_LEFT){if(x+w>s_pw-2)x=s_pw-2-w;if(x<2)x=2;}
    else if(anchor==A_RIGHT){if(x>s_pw-2)x=s_pw-2;if(x-w<2)x=2+w;}
    else{if(x+w/2>s_pw-2)x=s_pw-2-w/2;if(x-w/2<2)x=2+w/2;}
  }
  p->x[0]=x;p->y[0]=y;
}

typedef struct { float d; int i,j; } QUAD;

static int qcmp(const void *a,const void *b)
{
  float da=((const QUAD *)a)->d,db=((const QUAD *)b)->d;
  return (da<db)-(da>db);     /* far (large d) first */
}

static void set_view(void)
{
  double a=s_az*M_PI/180.0,e=s_el*M_PI/180.0;
  int m=s_pw<s_ph?s_pw:s_ph;
  ca=cos(a);sa=sin(a);ce=cos(e);se=sin(e);
  SC=s_zoom*0.33*m;
  CX=s_pw/2+14;        /* a little right of center: room for the z axis labels */
  CY=s_ph/2-6;
}

/* model point -> screen position and depth (larger d = farther away) */
static void project(double x,double y,double z,double *sx,double *sy,double *d)
{
  double u=x*ca-y*sa;
  double w=x*sa+y*ca;
  *sx=CX+SC*u;
  *sy=CY-SC*(z*ce+w*se);
  *d=w*ce-z*se;
}

/* a "nice" tick spacing (1, 2, 5 times a power of 10) giving about n ticks */
static double nice_step(double range,int n)
{
  double raw=range/n,mag=pow(10.0,floor(log10(raw))),f=raw/mag;
  return (f<1.5?1.0:f<3.0?2.0:f<7.0?5.0:10.0)*mag;
}

/* vertical z axis with ticks and numbers at floor corner (cx,cy) */
static void add_zaxis(SCENE *sc,double cx,double cy,double zlo,double zhi)
{
  double fz=-0.4*s_zs,tz=0.4*s_zs;
  double x0,y0,x1,y1,d,step,v,first,sx,sy;
  char buf[48];
  int k,nt;
  project(cx,cy,fz,&x0,&y0,&d);
  project(cx,cy,tz,&x1,&y1,&d);
  if(fabs(y0-y1)<30.0)return;          /* looking straight down: no room */
  add_line(sc,x0,y0,x1,y1);
  nt=(int)(fabs(y0-y1)/32.0);          /* about one tick per 32 pixels of axis */
  if(nt<2)nt=2;
  if(nt>6)nt=6;
  step=nice_step(zhi-zlo,nt);
  first=ceil(zlo/step-1e-9)*step;
  for(k=0,v=first;v<=zhi+1e-9*step&&k<12;v+=step,k++){
    double zv=v;
    if(fabs(zv)<1e-12*step)zv=0.0;
    project(cx,cy,((zv-zlo)/(zhi-zlo)-0.5)*0.8*s_zs,&sx,&sy,&d);
    add_line(sc,sx,sy,sx-5.0,sy);
    sprintf(buf,"%g",zv);
    add_text(sc,sx-8.0,sy+4.0,A_RIGHT,buf);
  }
  add_text(sc,x1,y1-8.0,A_CENTER,"z");
}

/* Build the whole picture.  fast: decimate a big grid (used while dragging);
   export: leave out the on-screen help line.  Returns 0 if no data. */
static int build_scene(int fast,int export,SCENE *sc)
{
  int nrows=my_browser.maxrow;
  int ncs=aplot.ncskip<1?1:aplot.ncskip;
  int nc=0,nr=0,i,j,ib,jb,st=1,n;
  int *ci=NULL,*rj=NULL,mc,mr,a,b,nq,k;
  double zlo=aplot.zmin,zhi=aplot.zmax,zn,x,y,z,sx,sy,d;
  double *zv=NULL,*dv=NULL,*px=NULL,*py=NULL;
  QUAD *q=NULL;
  char buf[300];
  double fz,cxs[4]={-1,1,1,-1},cys[4]={-1,-1,1,1},fsx[4],fsy[4],fd[4];
  int left=0,axis_back=0;
  double tlo,thi;
  PRIM *p;

  sc->p=NULL;sc->n=0;sc->cap=0;
  set_view();

  /* same grid the array plot uses */
  for(i=0;i<aplot.nacross/ncs;i++){
    ib=aplot.index0+i*ncs;
    if(ib>=my_browser.maxcol||ib<0)break;
    nc++;
  }
  for(j=0;j<aplot.ndown;j++){
    jb=aplot.nstart+aplot.nskip*j;
    if(jb>=nrows||jb<0)break;
    nr++;
  }
  if(nrows<=2||aplot.plotdef==0||nc<2||nr<2){
    add_text(sc,10,20,A_LEFT,"Nothing to plot yet: run a simulation and set up the array plot");
    return 0;
  }
  if(zhi<=zlo)zhi=zlo+1.0;

  nq=(nc-1)*(nr-1);
  if(fast&&nq>FAST_QUADS)st=(int)ceil(sqrt((double)nq/FAST_QUADS));

  ci=(int *)malloc(sizeof(int)*(nc+1));
  rj=(int *)malloc(sizeof(int)*(nr+1));
  if(ci==NULL||rj==NULL){free(ci);free(rj);return 0;}
  mc=0;for(i=0;i<nc;i+=st)ci[mc++]=i;
  if(ci[mc-1]!=nc-1)ci[mc++]=nc-1;
  mr=0;for(j=0;j<nr;j+=st)rj[mr++]=j;
  if(rj[mr-1]!=nr-1)rj[mr++]=nr-1;

  px=(double *)malloc(sizeof(double)*mc*mr);
  py=(double *)malloc(sizeof(double)*mc*mr);
  zv=(double *)malloc(sizeof(double)*mc*mr);
  dv=(double *)malloc(sizeof(double)*mc*mr);
  q=(QUAD *)malloc(sizeof(QUAD)*(mc-1)*(mr-1));
  if(px==NULL||py==NULL||zv==NULL||dv==NULL||q==NULL){
    free(ci);free(rj);free(px);free(py);free(zv);free(dv);free(q);return 0;
  }

  /* project every grid vertex once */
  for(a=0;a<mc;a++){
    ib=aplot.index0+ci[a]*ncs;
    x=2.0*ci[a]/(double)(nc-1)-1.0;
    for(b=0;b<mr;b++){
      jb=aplot.nstart+aplot.nskip*rj[b];
      y=1.0-2.0*rj[b]/(double)(nr-1);   /* time runs toward -y, so the top view matches the array plot */
      zn=(my_browser.data[ib][jb]-zlo)/(zhi-zlo);
      if(zn<0.0)zn=0.0;
      if(zn>1.0)zn=1.0;
      z=(zn-0.5)*0.8*s_zs;
      project(x,y,z,&sx,&sy,&d);
      k=a*mr+b;
      px[k]=sx;py[k]=sy;zv[k]=zn;dv[k]=d;
    }
  }

  /* floor outline, floor labels */
  fz=-0.4*s_zs;
  p=new_prim(sc,P_LINES);
  for(k=0;k<4;k++){
    project(cxs[k],cys[k],fz,&fsx[k],&fsy[k],&fd[k]);
    if(p){p->x[k]=fsx[k];p->y[k]=fsy[k];}
  }
  if(p){p->x[4]=p->x[0];p->y[4]=p->y[0];p->n=5;}
  /* z axis at the left-most floor corner; behind the surface if that corner is far */
  for(k=1;k<4;k++)if(fsx[k]<fsx[left])left=k;
  axis_back=(fd[left]>0.0);
  if(axis_back)add_zaxis(sc,cxs[left],cys[left],zlo,zhi);

  /* depth sort the quads */
  n=0;
  for(a=0;a<mc-1;a++)
    for(b=0;b<mr-1;b++){
      k=a*mr+b;
      q[n].d=(float)(0.25*(dv[k]+dv[k+1]+dv[k+mr]+dv[k+mr+1]));
      q[n].i=a;q[n].j=b;
      n++;
    }
  qsort(q,n,sizeof(QUAD),qcmp);

  /* paint from the back to the front */
  for(k=0;k<n;k++){
    int idx=q[k].i*mr+q[k].j;
    int colr;
    zn=0.25*(zv[idx]+zv[idx+1]+zv[idx+mr]+zv[idx+mr+1]);
    p=new_prim(sc,P_POLY);
    if(p==NULL)break;
    p->n=4;
    p->x[0]=px[idx];     p->y[0]=py[idx];
    p->x[1]=px[idx+mr];  p->y[1]=py[idx+mr];
    p->x[2]=px[idx+mr+1];p->y[2]=py[idx+mr+1];
    p->x[3]=px[idx+1];   p->y[3]=py[idx+1];
    if(COLOR){
      colr=FIRSTCOLOR+(int)(zn*color_total);
      if(colr<FIRSTCOLOR)colr=FIRSTCOLOR;
      if(colr>FIRSTCOLOR+color_total)colr=FIRSTCOLOR+color_total;
      p->colr=colr;
      p->stroke=s_mesh;
    }
    else{
      p->colr=C_BACK;
      p->stroke=1;
    }
  }

  if(!axis_back)add_zaxis(sc,cxs[left],cys[left],zlo,zhi);

  tlo=TRANS+aplot.nstart*NJMP*DELTA_T;
  thi=tlo+(nr-1)*aplot.nskip*NJMP*DELTA_T;

  /* axis names: on the lower (front) edge of each floor direction, just outside
     the floor outline, drawn last so the surface cannot cover them */
  {
    double mx=0.25*(fsx[0]+fsx[1]+fsx[2]+fsx[3]),my=0.25*(fsy[0]+fsy[1]+fsy[2]+fsy[3]);
    int e1=(fsy[0]+fsy[1]>=fsy[2]+fsy[3])?0:2;     /* edge y=-1 or y=+1: runs along index */
    int e2=(fsy[1]+fsy[2]>=fsy[3]+fsy[0])?1:3;     /* edge x=+1 or x=-1: runs along time  */
    int ee[2],t,c;
    ee[0]=e1;ee[1]=e2;
    for(t=0;t<2;t++){
      int k0=ee[t],k1=(ee[t]+1)%4;
      double ex=0.5*(fsx[k0]+fsx[k1]),ey=0.5*(fsy[k0]+fsy[k1]);
      double vx=ex-mx,vy=ey-my,len=sqrt(vx*vx+vy*vy);
      if(len<1e-6){vx=0;vy=1;len=1;}
      add_text(sc,ex+16.0*vx/len,ey+16.0*vy/len+4.0,A_CENTER,t==0?"index":"time");
      /* the two ends of the edge: first/last column, or start/end time.  Each goes
         just outside its own edge, so the two labels at a shared corner fan apart. */
      for(c=0;c<2;c++){
	int kc=c?k1:k0,an;
	double nx=vx/len,ny=vy/len,lx,ly;
	char lb[48];
	if(t==0)sprintf(lb,"%d",cxs[kc]<0?1:nc);
	else sprintf(lb,"%.5g",cys[kc]>0?tlo:thi);
	an=(nx<-0.3)?A_RIGHT:((nx>0.3)?A_LEFT:A_CENTER);
	if(kc==left){          /* the z axis is here: go below it */
	  lx=fsx[kc]+(an==A_RIGHT?-2.0:(an==A_LEFT?2.0:0.0));
	  ly=fsy[kc]+DCURYs+11.0;
	}
	else{lx=fsx[kc]+12.0*nx;ly=fsy[kc]+12.0*ny+4.0;}
	add_text(sc,lx,ly,an,lb);
      }
    }
  }

  /* annotations */
  if(export)sprintf(buf,"t = %g .. %g",tlo,thi);
  else sprintf(buf,"t = %g .. %g    az=%d el=%d",tlo,thi,(int)s_az,(int)s_el);
  add_text(sc,6,DCURYs,A_LEFT,buf);
  if(!export)
    add_text(sc,6,s_ph-4,A_LEFT,"drag: rotate   wheel/+/-: zoom   [ ]: height   m: mesh   r: reset");

  free(ci);free(rj);free(px);free(py);free(zv);free(dv);free(q);
  return 1;
}

/* ------------------------------------------------------------------ */
/*  Screen back end                                                    */
/* ------------------------------------------------------------------ */
static short clip16(double v)
{
  if(v>30000.0)return 30000;
  if(v<-30000.0)return -30000;
  return (short)v;
}

static int s_exporting=0;   /* 1 while drawing a picture for a GIF: white paper, black ink */

static unsigned long pixel_of(int colr)
{
  if(colr==C_BACK)return s_exporting?WhitePixel(display,screen):GrBack;
  if(colr==C_FORE)return s_exporting?BlackPixel(display,screen):GrFore;
  return ColorMap(colr);
}

static void paint_scene(SCENE *sc)
{
  int i,k;
  XPoint pt[6];
  for(i=0;i<sc->n;i++){
    PRIM *p=&sc->p[i];
    if(p->kind==P_TEXT){
      int len=strlen(p->t),w=len*DCURXs,x=(int)p->x[0];
      if(p->anchor==A_RIGHT)x-=w;
      else if(p->anchor==A_CENTER)x-=w/2;
      XSetForeground(display,s_gc,pixel_of(C_FORE));
      XDrawString(display,s_pix,s_gc,x,(int)p->y[0],p->t,len);
      continue;
    }
    for(k=0;k<p->n;k++){pt[k].x=clip16(p->x[k]);pt[k].y=clip16(p->y[k]);}
    if(p->kind==P_LINES){
      XSetForeground(display,s_gc,pixel_of(C_FORE));
      XDrawLines(display,s_pix,s_gc,pt,p->n,CoordModeOrigin);
    }
    else{
      XSetForeground(display,s_gc,pixel_of(p->colr));
      XFillPolygon(display,s_pix,s_gc,pt,p->n,Convex,CoordModeOrigin);
      if(p->stroke){
	pt[p->n]=pt[0];
	XSetForeground(display,s_gc,pixel_of(C_FORE));
	XDrawLines(display,s_pix,s_gc,pt,p->n+1,CoordModeOrigin);
      }
    }
  }
}

static void copy_to_window(void)
{
  if(!s_alive||s_pix==0)return;
  XCopyArea(display,s_pix,s_win,s_gc,0,0,s_pw,s_ph,0,0);
  XFlush(display);
}

static void make_pixmap(void)
{
  if(s_pix)XFreePixmap(display,s_pix);
  s_pix=XCreatePixmap(display,s_win,s_pw>0?s_pw:1,s_ph>0?s_ph:1,
		      DefaultDepth(display,screen));
}

static void render_ex(int fast,int export)
{
  SCENE sc;
  if(!s_alive||s_pix==0)return;
  s_exporting=export;
  XSetForeground(display,s_gc,pixel_of(C_BACK));
  XFillRectangle(display,s_pix,s_gc,0,0,s_pw,s_ph);
  build_scene(fast,export,&sc);
  paint_scene(&sc);
  free(sc.p);
  s_exporting=0;
  copy_to_window();
}

static void render(int fast)
{
  render_ex(fast,0);
}

/* ------------------------------------------------------------------ */
/*  Saved files                                                        */
/* ------------------------------------------------------------------ */
static double rgbc[1024][3];
static char   rgbok[1024];

static void rgb_of(int colr,double *r,double *g,double *b)
{
  XColor xc;
  if(colr==C_BACK){*r=*g=*b=1.0;return;}
  if(colr==C_FORE){*r=*g=*b=0.0;return;}
  if(colr>=0&&colr<1024&&rgbok[colr]){*r=rgbc[colr][0];*g=rgbc[colr][1];*b=rgbc[colr][2];return;}
  xc.pixel=ColorMap(colr);
  XQueryColor(display,DefaultColormap(display,screen),&xc);
  *r=xc.red/65535.0;*g=xc.green/65535.0;*b=xc.blue/65535.0;
  if(colr>=0&&colr<1024){rgbc[colr][0]=*r;rgbc[colr][1]=*g;rgbc[colr][2]=*b;rgbok[colr]=1;}
}

/* the file name and parameters, in XPP's usual footer, as plain text */
static char *pars_text(void)
{
  FILE *tf=tmpfile();
  char *txt;
  long len;
  if(tf==NULL)return NULL;
  ps_write_pars(tf);
  len=ftell(tf);
  rewind(tf);
  txt=(char *)malloc(len+1);
  if(txt){
    len=(long)fread(txt,1,len,tf);
    txt[len]=0;
  }
  fclose(tf);
  return txt;
}

static void ps_text_out(FILE *fp,const char *s)
{
  fputc('(',fp);
  for(;*s;s++){
    if(*s=='('||*s==')'||*s=='\\')fputc('\\',fp);
    fputc(*s,fp);
  }
  fputc(')',fp);
}

static int write_ps(const char *fn)
{
  SCENE sc;
  FILE *fp=fopen(fn,"w");
  int i,k,W=s_pw,H=s_ph;
  double r,g,b;
  char *pars;
  if(fp==NULL){err_msg("Cannot open file ");return 0;}
  memset(rgbok,0,sizeof(rgbok));
  build_scene(0,1,&sc);
  fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n%%%%BoundingBox: 0 0 %d %d\n",W,H);
  fprintf(fp,"%%%%Title: XPP surface plot\n%%%%EndComments\n");
  fprintf(fp,"/m {moveto} def /l {lineto} def\n");
  fprintf(fp,"/Q {closepath gsave fill grestore 0.4 setlinewidth stroke} def\n");
  fprintf(fp,"/B {closepath gsave fill grestore 0 setgray 0.4 setlinewidth stroke} def\n");
  fprintf(fp,"/Helvetica findfont 11 scalefont setfont\n1 setlinejoin\n");
  for(i=0;i<sc.n;i++){
    PRIM *p=&sc.p[i];
    if(p->kind==P_TEXT){
      fprintf(fp,"0 setgray %.2f %.2f m ",p->x[0],H-p->y[0]);
      if(p->anchor==A_RIGHT){ps_text_out(fp,p->t);fprintf(fp," dup stringwidth pop neg 0 rmoveto show\n");}
      else if(p->anchor==A_CENTER){ps_text_out(fp,p->t);fprintf(fp," dup stringwidth pop 2 div neg 0 rmoveto show\n");}
      else{ps_text_out(fp,p->t);fprintf(fp," show\n");}
    }
    else if(p->kind==P_LINES){
      fprintf(fp,"0 setgray 0.8 setlinewidth newpath %.2f %.2f m",p->x[0],H-p->y[0]);
      for(k=1;k<p->n;k++)fprintf(fp," %.2f %.2f l",p->x[k],H-p->y[k]);
      fprintf(fp," stroke\n");
    }
    else{
      rgb_of(p->colr,&r,&g,&b);
      fprintf(fp,"%.3f %.3f %.3f setrgbcolor newpath %.2f %.2f m",r,g,b,p->x[0],H-p->y[0]);
      for(k=1;k<p->n;k++)fprintf(fp," %.2f %.2f l",p->x[k],H-p->y[k]);
      fprintf(fp,p->stroke?" B\n":" Q\n");
    }
  }
  fprintf(fp,"showpage\n%%%%EOF\n");
  pars=pars_text();
  if(pars){fputs(pars,fp);free(pars);}
  fclose(fp);
  free(sc.p);
  return 1;
}

static void svg_text_out(FILE *fp,const char *s)
{
  for(;*s;s++){
    if(*s=='&')fputs("&amp;",fp);
    else if(*s=='<')fputs("&lt;",fp);
    else if(*s=='>')fputs("&gt;",fp);
    else fputc(*s,fp);
  }
}

static void svg_hex(int colr,char *out)
{
  double r,g,b;
  rgb_of(colr,&r,&g,&b);
  sprintf(out,"#%02x%02x%02x",(int)(r*255+0.5),(int)(g*255+0.5),(int)(b*255+0.5));
}

static int write_svg(const char *fn)
{
  SCENE sc;
  FILE *fp=fopen(fn,"w");
  int i,k,W=s_pw,H=s_ph;
  char hex[16],*pars,*c;
  if(fp==NULL){err_msg("Cannot open file ");return 0;}
  memset(rgbok,0,sizeof(rgbok));
  build_scene(0,1,&sc);
  fprintf(fp,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(fp,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n",W,H,W,H);
  fprintf(fp,"<title>XPP surface plot</title>\n");
  fprintf(fp,"<rect width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n",W,H);
  fprintf(fp,"<g stroke-linejoin=\"round\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"11\">\n");
  for(i=0;i<sc.n;i++){
    PRIM *p=&sc.p[i];
    if(p->kind==P_TEXT){
      fprintf(fp,"<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"%s\">",p->x[0],p->y[0],
	      p->anchor==A_RIGHT?"end":p->anchor==A_CENTER?"middle":"start");
      svg_text_out(fp,p->t);
      fprintf(fp,"</text>\n");
    }
    else if(p->kind==P_LINES){
      fprintf(fp,"<polyline fill=\"none\" stroke=\"#000000\" stroke-width=\"0.8\" points=\"");
      for(k=0;k<p->n;k++)fprintf(fp,"%s%.2f,%.2f",k?" ":"",p->x[k],p->y[k]);
      fprintf(fp,"\"/>\n");
    }
    else{
      svg_hex(p->colr,hex);
      /* a thin stroke in the fill color hides the hairline seams between quads */
      fprintf(fp,"<polygon fill=\"%s\" stroke=\"%s\" stroke-width=\"0.4\" points=\"",hex,p->stroke?"#000000":hex);
      for(k=0;k<p->n;k++)fprintf(fp,"%s%.2f,%.2f",k?" ":"",p->x[k],p->y[k]);
      fprintf(fp,"\"/>\n");
    }
  }
  fprintf(fp,"</g>\n</svg>\n");
  pars=pars_text();
  if(pars){
    for(c=pars;*c;c++)if(c[0]=='-'&&c[1]=='-')c[1]=' ';   /* "--" is not allowed in an XML comment */
    fprintf(fp,"<!--\n%s-->\n",pars);
    free(pars);
  }
  fclose(fp);
  free(sc.p);
  return 1;
}

static int write_gif(const char *fn)
{
  FILE *fp=fopen(fn,"w");
  char note[300];
  if(fp==NULL){err_msg("Cannot open file ");return 0;}
  render_ex(0,1);
  XSync(display,False);
  sprintf(note," %%%% Surface plot of the array plot: az=%g el=%g zoom=%g height=%g\n",s_az,s_el,s_zoom,s_zs);
  gif_set_note(note);
  screen_to_gif(s_pix,fp);
  fclose(fp);
  render(0);                /* back to the normal screen picture */
  return 1;
}

/* ask for a file name, then write */
static void save_as(int kind)
{
  char filename[XPP_MAX_NAME];
  static const char *ext[3]={"gif","ps","svg"};
  static const char *wild[3]={"*.gif","*.ps","*.svg"};
  static const char *title[3]={"Surface GIF","Surface PostScript","Surface SVG"};
  snprintf(filename,sizeof(filename),"%.*s.surface.%s",XPP_MAX_NAME-24,this_file,ext[kind]);
  if(!file_selector((char *)title[kind],filename,(char *)wild[kind]))return;
  if(kind==0)write_gif(filename);
  else if(kind==1)write_ps(filename);
  else write_svg(filename);
}

/* ------------------------------------------------------------------ */
/*  Window, buttons, events                                            */
/* ------------------------------------------------------------------ */
static void reset_view(void)
{
  s_az=30.0;s_el=35.0;s_zoom=1.0;s_zs=1.0;
}

static void draw_label(Window w)
{
  char *s=NULL;
  if(w==s_bredraw)s="Redraw";
  else if(w==s_bmesh)s=s_mesh?"Mesh:on":"Mesh:off";
  else if(w==s_breset)s="Reset";
  else if(w==s_bclose)s="Close";
  else if(w==s_bgif)s="GIF";
  else if(w==s_bps)s="PS";
  else if(w==s_bsvg)s="SVG";
  if(s)XDrawString(display,w,small_gc,0,CURY_OFFs,s,strlen(s));
}

void surface_expose(Window w)
{
  if(!s_alive)return;
  if(w==s_win){copy_to_window();return;}
  draw_label(w);
}

void surface_create(void)
{
  XGCValues vals;
  XGCValues fv;
  if(s_alive){
    XRaiseWindow(display,s_base);
    render(0);
    return;
  }
  s_top=2*DCURYs+16;
  s_base=make_plain_window(RootWindow(display,screen),0,0,s_w,s_h,1);
  XStoreName(display,s_base,"Surface");
  s_pw=s_w-4;
  s_ph=s_h-s_top-4;
  s_bredraw=br_button(s_base,0,0,"Redraw",0);
  s_bmesh=br_button(s_base,0,1,"Mesh",0);
  s_breset=br_button(s_base,0,2,"Reset",0);
  s_bclose=br_button(s_base,0,3,"Close",0);
  s_bgif=br_button(s_base,1,0,"GIF",0);
  s_bps=br_button(s_base,1,1,"PS",0);
  s_bsvg=br_button(s_base,1,2,"SVG",0);
  wm_register_close(s_base,s_bclose);
  s_win=make_plain_window(s_base,0,s_top,s_pw,s_ph,2);
  s_gc=XCreateGC(display,s_win,0,&vals);
  if(XGetGCValues(display,small_gc,GCFont,&fv))
    XSetFont(display,s_gc,fv.font);
  s_pix=0;
  make_pixmap();
  {
    int ncs=aplot.ncskip<1?1:aplot.ncskip;
    int ncol=aplot.nacross/ncs,nrow=aplot.ndown;
    s_mesh=((ncol-1)*(nrow-1)<=MESH_AUTO_MAX);
  }
  s_alive=1;
  render(0);
}

void surface_close(void)
{
  if(!s_alive)return;
  s_alive=0;
  if(s_pix){XFreePixmap(display,s_pix);s_pix=0;}
  XFreeGC(display,s_gc);
  XDestroySubwindows(display,s_base);
  XDestroyWindow(display,s_base);
}

void surface_redraw(void)
{
  if(s_alive)render(0);
}

static void zoom_by(double f)
{
  s_zoom*=f;
  if(s_zoom<0.2)s_zoom=0.2;
  if(s_zoom>8.0)s_zoom=8.0;
  render(0);
}

static int is_button(Window w)
{
  return w==s_bredraw||w==s_bmesh||w==s_breset||w==s_bclose||
         w==s_bgif||w==s_bps||w==s_bsvg;
}

static void button_hit(Window w)
{
  if(w==s_bredraw)render(0);
  else if(w==s_bmesh){
    s_mesh=!s_mesh;
    XClearArea(display,s_bmesh,0,0,0,0,True);
    render(0);
  }
  else if(w==s_breset){reset_view();render(0);}
  else if(w==s_bgif)save_as(0);
  else if(w==s_bps)save_as(1);
  else if(w==s_bsvg)save_as(2);
  else if(w==s_bclose)surface_close();
}

void do_surface_events(XEvent ev)
{
  Window w;
  int dx,dy,nw,nh;
  KeySym ks;
  char buf[16];
  if(!s_alive)return;
  w=ev.xany.window;
  switch(ev.type){
  case ConfigureNotify:
    if(w!=s_base)return;
    nw=ev.xconfigure.width;
    nh=ev.xconfigure.height;
    if(nw==s_w&&nh==s_h)return;
    s_w=nw;s_h=nh;
    s_pw=s_w-4;
    s_ph=s_h-s_top-4;
    if(s_pw<50||s_ph<50)return;
    XResizeWindow(display,s_win,s_pw,s_ph);
    make_pixmap();
    render(0);
    break;
  case EnterNotify:
  case LeaveNotify:
    if(is_button(w))
      XSetWindowBorderWidth(display,w,ev.type==EnterNotify?2:1);
    break;
  case ButtonPress:
    if(w==s_win){
      if(ev.xbutton.button==Button4){zoom_by(1.1);break;}
      if(ev.xbutton.button==Button5){zoom_by(1.0/1.1);break;}
      s_drag=1;
      s_lx=ev.xbutton.x;
      s_ly=ev.xbutton.y;
    }
    else button_hit(w);
    break;
  case MotionNotify:
    if(w==s_win&&s_drag){
      dx=ev.xmotion.x-s_lx;
      dy=ev.xmotion.y-s_ly;
      s_lx=ev.xmotion.x;
      s_ly=ev.xmotion.y;
      s_az+=0.6*dx;
      s_el+=0.6*dy;
      if(s_el>90.0)s_el=90.0;
      if(s_el<-90.0)s_el=-90.0;
      while(s_az>180.0)s_az-=360.0;
      while(s_az<-180.0)s_az+=360.0;
      render(1);
    }
    break;
  case ButtonRelease:
    if(w==s_win&&s_drag){
      s_drag=0;
      render(0);
    }
    break;
  case KeyPress:
    if(w!=s_win&&w!=s_base)return;
    if(XLookupString(&ev.xkey,buf,sizeof(buf)-1,&ks,NULL)<1)return;
    switch(buf[0]){
    case '+':case '=':zoom_by(1.15);break;
    case '-':zoom_by(1.0/1.15);break;
    case ']':s_zs*=1.2;if(s_zs>10.0)s_zs=10.0;render(0);break;
    case '[':s_zs/=1.2;if(s_zs<0.05)s_zs=0.05;render(0);break;
    case 'm':s_mesh=!s_mesh;XClearArea(display,s_bmesh,0,0,0,0,True);render(0);break;
    case 'r':reset_view();render(0);break;
    }
    break;
  }
}
