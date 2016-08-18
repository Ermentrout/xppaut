#include "matrixalgebra.h"

#include <math.h>

#include "../markov.h"
#include "../my_rhs.h"
#include "../numerics.h"


/* --- Forward declarations --- */
static void get_band_jac(double *a, double *y, double t, double *ypnew, double *ypold, int n, double eps, double scal);
static void hqrx(int n, int low, int igh, double *h, double *ev, int *ierr);
static int isamax(int n, double *sx, int incx);
static void orthesx(int n, int low, int igh, double *a, double *ort);
static void saxpy(int n, double sa, double *sx, int incx, double *sy, int incy);
static double sdot(int n, double *sx, int incx, double *sy, int incy);
static void sscal(int n, double sa, double *sx, int incx);


/* --- Functions --- */
/* factors the matrix */
int bandfac(double *a, int ml, int mr, int n) {
	int i,j,k;
	int row,rowi,m,r0,ri0;
	double al;

	int n1=n-1;
	int mt=ml+mr+1;

	for(row=0;row<n;row++) {
		r0=row*mt+ml;
		if((al=a[r0])==0.0) {
			return(-1-row);
		}
		al=1.0/al;
		m=MIN(mr,n1-row);
		for(j=1;j<=m;j++) {
			a[r0+j]=a[r0+j]*al;
		}
		a[r0]=al;
		for(i=1;i<=ml;i++) {
			rowi=row+i;
			if(rowi>n1) {
				break;
			}
			ri0=rowi*mt+ml;
			al=a[ri0-i];
			if(al==0.0) {
				continue;
			}
			for(k=1;k<=m;k++) {
				a[ri0-i+k]=a[ri0-i+k]-(al*a[r0+k]);
			}
			a[ri0-i]=-al;
		}
	}
	return(0);
}


/* requires that the matrix be factored */
void bandsol(double *a, double *b, int ml, int mr, int n) {
	int i,j,k,r0;

	int mt=ml+mr+1;
	int m,n1=n-1,row;

	for(i=0;i<n;i++) {
		r0=i*mt+ml;
		m=MAX(-ml,-i);
		for(j=m;j<0;j++) {
			b[i] += a[r0+j]*b[i+j];
		}
		b[i] *= a[r0];
	}
	for(row=n1-1;row>=0;row--) {
		m=MIN(mr,n1-row);
		r0=row*mt+ml;
		for(k=1;k<=m;k++) {
			b[row]=b[row]-a[r0+k]*b[row+k];
		}
	}
}


void eigen(int n, double *a, double *ev, double *work, int *ierr) {
	orthesx(n,1,n,a,work);
	hqrx(n,1,n,a,ev,ierr);
}


void get_evec(double *a, double *anew, double *b, double *bp, int n, int maxit,
			  double err, int *ipivot, double eval, int *ierr) {
	int j,iter,jmax;
	double temp;
	double zz=fabs(eval);
	if(zz<err) {
		zz=err;
	}
	*ierr=0;
	for(j=0;j<n*n;j++) {
		anew[j]=a[j];
	}
	for(j=0;j<n;j++) {
		anew[j*(1+n)]=anew[j*(1+n)]-eval-err*err*zz;
	}
	sgefa(anew,n,n,ipivot,ierr);
	if(*ierr!=-1) {
		printf(" Pivot failed\n");
		return;
	}
	for(j=0;j<n;j++) {
		b[j]=1+.1*ndrand48();
		bp[j]=b[j];
	}
	iter=0;
	*ierr=0;
	while(1) {
		sgesl(anew,n,n,ipivot,b,0);
		temp=fabs(b[0]);
		jmax=0;

		for(j=0;j<n;j++) {
			if(fabs(b[j])>temp) {
				temp=fabs(b[j]);
				jmax=j;
			}
		}
		temp=b[jmax];
		for(j=0;j<n;j++) {
			b[j]=b[j]/temp;
		}
		temp=0.0;
		for(j=0;j<n;j++) {
			temp=temp+fabs(b[j]-bp[j]);
			bp[j]=b[j];
		}
		if(temp<err) {
			break;
		}
		iter++;
		if(iter>maxit) {
			printf(" max iterates exceeded\n");
			*ierr=1;
			break;
		}
	}
	if(*ierr==0) {
		temp=fabs(b[0]);
		jmax=0;
		for(j=0;j<n;j++) {
			if(fabs(b[j])>temp) {
				temp=fabs(b[j]);
				jmax=j;
			}
		}
		temp=b[jmax];
		for(j=0;j<n;j++) {
			b[j]=b[j]/temp;
		}
	}
	return;
}


void getjac(double *x, double *y, double *yp, double *xp, double eps,
			double *dermat, int n) {
	int i,j,k;
	double r;
	my_rhs(0.0,x,y,n);
	if(METHOD==0) {
		for(i=0;i<n;i++) {y[i]=y[i]-x[i];
		}
	}
	for(i=0;i<n;i++) {
		for(k=0;k<n;k++) {
			xp[k]=x[k];
		}
		r=eps*MAX(eps,fabs(x[i]));
		xp[i]=xp[i]+r;
		my_rhs(0.0,xp,yp,n);

		if(METHOD==0) {
			for(j=0;j<n;j++) {
				yp[j]=yp[j]-xp[j];
			}
		}
		for(j=0;j<n;j++) {
			dermat[j*n+i]=(yp[j]-y[j])/r;
		}
	}
}


void getjactrans(double *x,double *y,double *yp,double *xp, double eps, double *dermat, int n) {
	int i,j,k;
	double r;
	my_rhs(0.0,x,y,n);
	for(i=0;i<n;i++) {
		for(k=0;k<n;k++) {
			xp[k]=x[k];
		}
		r=eps*MAX(eps,fabs(x[i]));
		xp[i]=xp[i]+r;
		my_rhs(0.0,xp,yp,n);

		for(j=0;j<n;j++) {
			dermat[j+n*i]=(yp[j]-y[j])/r;
		}
	}
}


/* this assumes that yp is already computed */
void get_the_jac(double t, double *y, double *yp, double *ypnew, double *dfdy, int neq, double eps, double scal) {
	int i,j;
	double yold,del,dsy;
	if(cv_bandflag) {
		get_band_jac(dfdy,y,t,ypnew,yp,neq,eps,scal);
	} else {
		for(i=0;i<neq;i++) {
			del=eps*MAX(eps,fabs(y[i]));
			dsy=scal/del;
			yold=y[i];
			y[i]=y[i]+del;
			my_rhs(t,y,ypnew,neq);
			for(j=0;j<neq;j++) {
				dfdy[j*neq+i]=dsy*(ypnew[j]-yp[j]);
			}
			y[i]=yold;
		}
	}
}


/**
 * Factors a real matrix by Gaussian elimination.
 *
 * @param a (in/out) the matrix to be factored.
 * @param lda the leading dimension of the array A.
 * @param n the order of the matrix A.
 * @param ipvt (out) an integer vector of pivot indices.
 * @param info (out) = 0  normal value.
 *             = k  if  u(k,k) .eq. 0.0 .  This is not an error
 *                  condition for this subroutine, but it does
 *                  indicate that SGESL or SGEDI will divide by zero
 *                  if called.  Use  RCOND  in SGECO for a reliable
 *                  indication of singularity.
 */
void sgefa(double *a, int lda, int n, int *ipvt, int *info) {
	int j,k,kp1,l,nm1;
	double t;
	*info=-1;
	nm1=n-1;
	if(nm1>0) {
		for(k=1;k<=nm1;k++) {
			kp1=k+1;
			l=isamax(n-k+1,&a[(k-1)*lda+k-1],lda)+k-1;
			ipvt[k-1]=l;
			if(a[l*lda+k-1]!=0.0) {
				if(l!=(k-1)) {
					t=a[l*lda+k-1];
					a[l*lda+k-1]=a[(k-1)*lda+k-1];
					a[(k-1)*lda+k-1]=t;
				}
				t=-1.0/a[(k-1)*lda+k-1];
				sscal(n-k,t,(a+k*lda+k-1),lda);
				for(j=kp1;j<=n;j++)	{
					t=a[l*lda+j-1];
					if(l!=(k-1)) {
						a[l*lda+j-1]=a[(k-1)*lda+j-1];
						a[(k-1)*lda+j-1]=t;
					}
					saxpy(n-k,t,(a+k*lda+k-1),lda,(a+k*lda+j-1),lda);
				}
			} else {
				*info=k-1;
			}
		}
	}
	ipvt[n-1]=n-1;
	if(a[(n-1)*lda+n-1]==0.0) {
		*info=n-1;
	}
}


/**
 * Solves the real system
 *   A * X = B  or  transpose(A) * X = B
 * using factors computed by sgeco or sgefa.
 *
 * @param a the output from sgeco or sgefa.
 * @param lda the leading dimension of the array A.
 * @param n the order of the matrix A.
 * @param ipvt the pivot vector from sgeco or sgefa.
 * @param b (in/out) the right hand side vector.
 * @param job whether to solve the transpose(A) problem.
 */
void sgesl(double *a, int lda, int n, int *ipvt, double *b, int job) {
	int k,kb,l,nm1;
	double t;
	nm1=n-1;
	if(job==0) {
		if(nm1>=1) {
			for(k=1;k<=nm1;k++) {
				l=ipvt[k-1];
				t=b[l];
				if(l!=(k-1)) {
					b[l]=b[k-1];
					b[k-1]=t;
				}
				saxpy(n-k,t,(a+lda*k+k-1),lda,(b+k),1);
			}
		}
		for(kb=1;kb<=n;kb++) {
			k=n+1-kb;
			b[k-1]=b[k-1]/a[(k-1)*lda+k-1];
			t=-b[k-1];
			saxpy(k-1,t,(a+k-1),lda,b,1);
		}
		return;
	}
	for(k=1;k<=n;k++) {
		t=sdot(k-1,(a+k-1),lda,b,1);
		b[k-1]=(b[k-1]-t)/a[(k-1)*lda+k-1];
	}
	if(nm1>0) {
		for(kb=1;kb<=nm1;kb++) {
			k=n-kb;
			b[k-1]=b[k-1]+sdot(n-k,(a+k*lda+k-1),lda,b+k,1);
			l=ipvt[k-1];
			if(l!=(k-1)) {
				t=b[l];
				b[l]=b[k-1];
				b[k-1]=t;
			}
		}
	}
}


double sqr2(double z) {
	return(z*z);
}



/* --- Static functions --- */
static void get_band_jac(double *a, double *y, double t, double *ypnew, double *ypold, int n, double eps, double scal) {
	int ml=cv_bandlower,mr=cv_bandupper;
	int i,j,k,n1=n-1,mt=ml+mr+1;
	double yhat;
	double dy;
	double dsy;

	for(i=0;i<(n*mt);i++) {
		a[i]=0.0;
	}
	for(i=0;i<n;i++) {
		yhat=y[i];
		dy=eps*(eps+fabs(yhat));
		dsy=scal/dy;
		y[i] += dy;
		my_rhs(t,y,ypnew,n);
		for(j=-ml;j<=mr;j++) {
			k=i-j;
			if(k<0 || k>n1) {
				continue;
			}
			a[k*mt+j+ml]=dsy*(ypnew[k]-ypold[k]);
		}
		y[i]=yhat;
	}
}


static void hqrx(int n, int low, int igh, double *h, double *ev, int *ierr) {
	int i,j,k,l=0,m=0,en,ll,mm,na,its,mp2,enm2;
	double p=0.0,q=0.0,r=0.0,s,t,w,x,y,zz,norm,machep=1.e-10;
	int notlas;
	*ierr = 0;
	norm = 0.0;
	k = 1;
	for( i = 1;i<= n;i++) {
		for(j = k;j<= n;j++) {
			norm = norm + fabs(h[i-1+(j-1)*n]);
		}
		k = i;
		if((i >= low) && ( i <= igh)) {
			continue;
		}
		ev[(i-1)*2] = h[i-1+(i-1)*n];
		ev[1+(i-1)*2] = 0.0;
	}
	en = igh;
	t = 0.0;
l60:
	if(en < low) {
		return;
	}
	its = 0;
	na = en - 1;
	enm2 = na - 1;
l70:
	for(ll=low; ll<=en; ll++) {
		l = en + low - ll;
		if(l==low) {
			break;
		}
		s = fabs(h[l-2+(l-2)*n]) + fabs(h[l-1+(l-1)*n]);
		if(s==0.0) {
			s = norm;
		}
		if(fabs(h[l-1+(l-2)*n]) <= machep*s) {
			break;
		}
	}
	x = h[en-1+(en-1)*n];
	if(l==en) {
		goto l270;
	}
	y = h[na-1+(na-1)*n];
	w = h[en-1+(na-1)*n] * h[na-1+(en-1)*n];
	if(l==na) {
		goto l280;
	}
	if(its==30) {
		goto l1000;
	}
	if((its != 10) && (its != 20)) {
		goto l130;
	}
	t = t + x;
	for(i = low;i<= en;i++) {
		h[i-1+(i-1)*n] = h[i-1+(i-1)*n] - x;
	}
	s = fabs(h[en-1+(na-1)*n]) + fabs(h[na-1+(enm2-1)*n]);
	x = 0.75 * s;
	y = x;
	w = -0.4375 * s * s;
l130:
	its++; /*its = its++; This may be undefined. Use its++ instead.*/
	for(mm = l;mm <= enm2;mm++) {
		m = enm2 + l - mm;
		zz = h[m-1+(m-1)*n];
		r = x - zz;
		s = y - zz;
		p = (r * s - w) / h[m+(m-1)*n] + h[m-1+m*n];
		q = h[m+m*n] - zz - r - s;
		r = h[m+1+m*n];
		s = fabs(p) + fabs(q) + fabs(r);
		p = p / s;
		q = q / s;
		r = r / s;
		if(m==l) {
			break;
		}
		if((fabs(h[m-1+(m-2)*n])*(fabs(q)+fabs(r)))<=
		   (machep*fabs(p)*(fabs(h[m-2+(m-2)*n])+ fabs(zz) + fabs(h[m+m*n])))) {
			break;
		}
	}
	mp2 = m + 2;
	for( i = mp2;i<= en;i++) {
		h[i-1+(i-3)*n] = 0.0;
		if(i==mp2) {
			continue;
		}
		h[i-1+(i-4)*n] = 0.0;
	}
	/*260 */
	for( k = m;k<= na;k++) {
		notlas=0;
		if(k != na) {
			notlas=1;
		}
		if(k==m) {
			goto l170;
		}
		p = h[k-1+(k-2)*n];
		q = h[k+(k-2)*n];
		r = 0.0;
		if(notlas) {
			r = h[k+1+(k-2)*n];
		}
		x=fabs(p) + fabs(q) + fabs(r);
		if(x==0.0) {
			continue;
		}
		p = p / x;
		q = q / x;
		r = r / x;
l170:
		s = SIGN(sqrt(p*p+q*q+r*r),p);
		if(k != m) {
			h[k-1+(k-2)*n] = -s * x;
		} else if(l != m) {
			h[k-1+(k-2)*n] = -h[k-1+(k-2)*n];
		}
		p = p + s;
		x = p / s;
		y = q / s;
		zz = r / s;
		q = q / p;
		r = r / p;
		for(j = k;j<= en;j++) {
			p = h[k-1+(j-1)*n] + q * h[k+(j-1)*n];
			if(notlas) {
				p = p + r * h[k+1+(j-1)*n];
				h[k+1+(j-1)*n] = h[k+1+(j-1)*n] - p * zz;
			}
			h[k+(j-1)*n] = h[k+(j-1)*n] - p * y;
			h[k-1+(j-1)*n] = h[k-1+(j-1)*n] - p * x;
		}
		j = MIN(en,k+3);
		for(i = l;i<= j ;i++) {
			p = x * h[i-1+(k-1)*n] + y * h[i-1+k*n];
			if(notlas) {
				p = p + zz * h[i-1+(k+1)*n];
				h[i-1+(k+1)*n] = h[i-1+(k+1)*n] - p * r;
			}
			h[i-1+k*n] = h[i-1+k*n] - p * q;
			h[i-1+(k-1)*n] = h[i-1+(k-1)*n] - p;
		}
	}
	goto l70;
l270:
	ev[(en-1)*2]=x+t;
	ev[1+(en-1)*2]=0.0;
	en = na;
	goto l60;
l280:
	p = (y - x) / 2.0;
	q = p * p + w;
	zz = sqrt(fabs(q));
	x = x + t;
	if(q < 0.0) {
		goto l320;
	}
	zz = p + SIGN(zz,p);
	ev[(na-1)*2] = x + zz;
	ev[(en-1)*2] = ev[(na-1)*2];
	if(zz != 0.0) {
		ev[(en-1)*2] = x-w/zz;
	}
	ev[1+(na-1)*2] = 0.0;
	ev[1+(en-1)*2] = 0.0;
	goto l330;
l320:
	ev[(na-1)*2] = x+p;
	ev[(en-1)*2] = x+p;
	ev[1+(na-1)*2] = zz;
	ev[1+(en-1)*2] = -zz;
l330:
	en = enm2;
	goto l60;
l1000:
	*ierr = en;
}


static int isamax(int n, double *sx, int incx) {
	int i,ix,imax;
	double smax;
	if(n<1)
		return(-1);
	if(n==1)
		return(0);
	if(incx!=1) {
		ix=0;
		imax=0;
		smax=fabs(sx[0]);
		ix+=incx;
		for(i=1;i<n;i++,ix+=incx) {
			if(fabs(sx[ix])>smax) {
				imax=i;
				smax=fabs(sx[ix]);
			}
		}
		return(imax);
	}
	imax=0;
	smax=fabs(sx[0]);
	for(i=1;i<n;i++) {
		if(fabs(sx[i])>smax) {
			imax=i;
			smax=fabs(sx[i]);
		}
	}
	return(imax);
}


static void orthesx(int n, int low, int igh, double *a, double *ort) {
	int i,j,m,ii,jj,la,mp,kp1;
	double f,g,h,scale;
	la = igh - 1;
	kp1 = low + 1;
	if(la < kp1) {
		return;
	}
	/*180*/
	for(m=kp1;m<=la;m++) {
		h = 0.0;
		ort[m-1] = 0.0;
		scale = 0.0;
		for(i = m;i<= igh;i++) {
			scale = scale + fabs(a[i-1+(m-2)*n]);
		}
		if(scale==0.0) {
			continue;
		}
		mp = m + igh;
		/*100*/
		for( ii = m;ii<= igh;ii++) {
			i = mp - ii;
			ort[i-1] = a[i-1+(m-2)*n] / scale;
			h = h + ort[i-1] * ort[i-1];
		}
		g = -SIGN(sqrt(h),ort[m-1]);
		h = h - ort[m-1] * g;
		ort[m-1] = ort[m-1] - g;
		/*130 */
		for(j = m;j<= n;j++) {
			f = 0.0;
			for( ii = m;ii<= igh;ii++) {
				i = mp - ii;
				f = f + ort[i-1] * a[i-1+(j-1)*n];
			}
			f = f / h;
			for(i = m;i<= igh;i++) {
				a[i-1+(j-1)*n] = a[i-1+(j-1)*n] - f * ort[i-1];
			}
		}
		/*160*/
		for(i = 1;i<= igh;i++) {
			f = 0.0;
			/*140 */
			for( jj = m;jj<= igh;jj++) {
				j = mp - jj;
				f = f + ort[j-1] * a[i-1+(j-1)*n];
			}
			f = f / h;
			for(j = m;j<= igh;j++) {
				a[i-1+(j-1)*n] = a[i-1+(j-1)*n] - f * ort[j-1];
			}
		}
		ort[m-1] = scale * ort[m-1];
		a[m-1+(m-2)*n] = scale * g;
	}
}


static void saxpy(int n, double sa, double *sx, int incx, double *sy, int incy) {
	int i,ix,iy;
	if(n<=0)
		return;
	if(sa==0.0)
		return;
	ix=0;
	iy=0;
	if(incx<0) {
		ix=-n*incx;
	}
	if(incy<0) {
		iy=-n*incy;
	}
	for(i=0;i<n;i++,ix+=incx,iy+=incy) {
		sy[iy]=sy[iy]+sa*sx[ix];
	}
}


static double sdot(int n, double *sx, int incx, double *sy, int incy) {
	int i,ix,iy;
	double stemp=0.0;
	if(n<=0)
		return(0.0);
	ix=0;
	iy=0;
	if(incx<0) {
		ix=-n*incx;
	}
	if(incy<0) {
		iy=-n*incy;
	}
	for(i=0;i<n;i++,ix+=incx,iy+=incy) {
		stemp+=sx[ix]*sy[iy];
	}
	return(stemp);
}


static void sscal(int n, double sa, double *sx, int incx) {
	int i,nincx;
	if(n<=0)
		return;
	nincx=n*incx;
	for(i=0;i<nincx;i+=incx) {
		sx[i]*=sa;
	}
}
