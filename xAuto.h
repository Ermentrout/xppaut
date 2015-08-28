 typedef struct {
    int ndim;
    int ips;
    int irs;
    int ilp;
    int nicp;
    int icp[100];
    int ntst;
    int ncol;
    int iad;
    int isp;
    int isw;
    int iplt;
    int nbc;
    int nint;
    int nmx;
    double rl0;
    double rl1;
    double a0;
    double a1;
    int npr;
    int mxbf;
    int iid;
    int itmx;
    int itnw;
    int nwtn;
    int jac;
    double epsl;
    double epss;
    double epsu;
    double ds;
    double dsmax;
    double dsmin;
    int iads;
    int nthl;
    int ithl[100];
    double thl[100];
    int nuzr;
    int iuz[20];
    double vuz[20];
   /*  these are for the homoclinic stuff */
    int nunstab;
   int nstab;
   int iequib; /* +1 homoclinic   -2 heteroclinic */



  } XAUTO;

/*  for homcont  - itwist=0, istart=2, nrev=0,nfixed=0,npsi=0 */

 
