// Routine contains the structures for reading the user input through
// input_fform.cc. The default values for all these variables are defined in
// input.config.
//
// Input variables that have been previously declared in common.h have to be
// re-declared here, in a consistant structure.   
/* Arsen
#ifdef intel
#define workfc WORKFC 
#define fronts FRONTS 
#define newdim NEWDIM
#define timer4 TIMER4 
#define extrat EXTRAT 
#define spongevar SPONGEVAR 
#define turbvar TURBVAR 
#define turbvari TURBVARI 
#define spebcvr SPEBCVR 
#define aerfrc AERFRC 
#define astore ASTORE 
#define conpar CONPAR 
#define shpdat SHPDAT 
#define datpnt DATPNT 
#define errpar ERRPAR
#define elmpar ELMPAR 
#define genpar GENPAR 
#define inpdat INPDAT 
#define intdat INTDAT 
#define mio MIO
#define mioname MIONAME 
#define itrpar ITRPAR 
#define itrpnt ITRPNT 
#define matdat MATDAT
#define matramp MATRAMP 
#define mmatpar MMATPAR 
#define outpar OUTPAR 
#define point POINT
#define precis PRECIS 
#define propar PROPAR 
#define resdat RESDAT 
#define solpar SOLPAR 
#define timdat TIMDAT 
#define timpar TIMPAR
#define incomp INCOMP
#define mtimer1 MTIMER1 
#define mtimer2 MTIMER2
#define timer3 TIMER3 
#define title TITLE
#define sclrs SCLRS
#define levlset LEVLSET
#define bubstudy BUBSTUDY
#define nomodule NOMODULE
#define sequence SEQUENCE
#define amgvarr AMGVARR
#define amgvari AMGVARI
#define pcboiling PCBOILING
#define contactangle CONTACTANGLE
#elif ( !defined ibm)*/

#include <FCMangle.h>

#define workfc workfc_
#define fronts fronts_ 
#define newdim newdim_ 
#define timer4 timer4_ 
#define extrat extrat_ 
#define spongevar spongevar_ 
#define turbvar turbvar_ 
#define turbvari turbvari_ 
#define spebcvr spebcvr_ 
#define sclrs sclrs_
#define aerfrc aerfrc_ 
#define astore astore_ 
#define conpar conpar_ 
#define levlset levlset_
#define bubstudy bubstudy_
#define shpdat shpdat_ 
#define datpnt datpnt_
#define elmpar elmpar_ 
#define errpar errpar_
#define genpar genpar_ 
#define inpdat inpdat_
#define mio mio_ 
#define mioname mioname_ 
#define itrpar itrpar_ 
#define itrpnt itrpnt_ 
#define matdat matdat_ 
#define matramp matramp_
#define mmatpar mmatpar_ 
#define outpar outpar_ 
#define point point_ 
#define precis precis_ 
#define propar propar_ 
#define resdat resdat_ 
#define solpar solpar_ 
#define timdat timdat_ 
#define timpar timpar_ 
#define incomp incomp_ 
#define mtimer1 mtimer1_ 
#define mtimer2 mtimer2_ 
#define timer3 timer3_ 
#define title title_ 
#define intdat intdat_ 
#define nomodule nomodule_
#define pcboiling pcboiling_
#define contactangle contactangle_
#define sequence sequence_
#define amgvarr amgvarr_
#define amgvari amgvari_
//#endif

#define MAXBLK   5000
#define MAXSURF  20  
#define MAXTS    100
#define MAXTOP   5
#define MAXQPT   125
#define MAXSH    125
#define NSD      3
#define machin   'RS/6000'
#define machfl   4
#define zero   0.0000000000000000000000000000000d0
#define pt125  0.1250000000000000000000000000000d0
#define pt25   0.2500000000000000000000000000000d0
#define pt33   0.3333333333333333333333333333333d0
#define pt39   0.3968502629920498686879264098181d0
#define pt5    0.5000000000000000000000000000000d0
#define pt57   0.5773502691896257645091487805020d0
#define pt66   0.6666666666666666666666666666667d0
#define pt75   0.7500000000000000000000000000000d0
#define one    1.0000000000000000000000000000000d0
#define sqt2   1.4142135623730950488016887242097d0
#define onept5 1.5000000000000000000000000000000d0
#define two    2.0000000000000000000000000000000d0
#define three  3.0000000000000000000000000000000d0
#define four   4.0000000000000000000000000000000d0
#define five   5.0000000000000000000000000000000d0
#define pi     3.1415926535897932384626433832795d0

#ifdef __cplusplus
extern "C" {
#endif
  extern struct { 
    int master;
    int numpe;
    int myrank;
  } workfc;

  extern struct { 
    int maxfront;
    int nlwork;
    int idirstep;
    int idirtrigger;
  } fronts;

  extern struct { 
    int numper;
    int nshgt;
    int nshg0;
  } newdim;

  extern struct { 
    double birth;
    double death;
    double comtim;
  } timer4;

  extern struct { 
    double ttim[100];
  } extrat;

  extern struct {
    double zoutsponge;
    double radsponge;
    double zinsponge;
    double grthosponge;
    double grthisponge;
    double betamax;
    int spongecontinuity;
    int spongemomentum1;
    int spongemomentum2;
    int spongeenergy;
    int spongemomentum3;
  } spongevar;

  extern struct {
    double eles;
    double ylimit[9][3]; /* 9 = 5 + 4 = puvwT + 4Scalars */
    double rmutarget;
    double pzero;
    double wtavei;
    double dtavei;
    double dke;
    double fwr1;
    double flump;
    int ierrcalc;
    int ihessian;
    int itwmod;
    int ngaussf;
    int idim;
    int nlist;
    int nintf[MAXTOP];
  } turbvar;

  extern struct {
    int irans;
    int iles;
    int idns;
    int isubmod;
    int ifproj;
    int i2filt;
    int modlstats;
    int idis;
    int nohomog;
    int ierrsmooth;
  } turbvari;

  extern struct { 
    int irscale;
    int intpres;
    double plandist;
    double thetag;
    double ds;
    double tolerence;
    double radcyl;
    double rbltin;
    double rvscal;
  } spebcvr;

  extern struct {
    double scdiff[5];
    double tdecay;
    int nsclr;
    int isclr;
    int nsolt;
    int nosource;
    int consrv_sclr_conv_vel;
  } sclrs;

  extern struct { 
    double flxID[MAXSURF+1][10];
    double Force[3];
    double HFlux;
    int nsrflist[MAXSURF+1];
    int isrfIM;
    double flxIDsclr[MAXSURF][4];
  } aerfrc;

  extern struct { 
    double a[100000];
  } astore;

  extern struct { 
    int numnp;
    int numel;
    int numelb;
    int numpbc;
    int nen;
    int nfaces;
    int numflx;
    int ndof;
    int iALE;
    int icoord;
    int navier;
    int iblk;
    int irs;
    int iexec;
    int necho;
    int ichem;
    int iRK;
    int nedof;
    int nshg;
    int nnz;
    int istop;
    int nflow;
    int nnz_tot;
    int idtn;
    int ncorpsize;    // Arsen
    int iownnodes;    // Arsen
  } conpar;

  extern struct { 
    double epsilon_ls;
    double epsilon_lsd;
    double dtlset;
    double dtlset_cfl;
    double redist_toler;
    double redist_toler_curr;
    double r_int_buffer;
    double r_int_elem_size;
    double phvol[2];
    double AdjRedistVelCFL;
    double BubRad;
    double vf_target;
    double C_int_adjust;
    double vf_now;
    double vfcontrcoeff;
    double C_int_cap;
    double epsilonBT;
    double xdistancesum;
    double ydistancesum;
    double zdistancesum;
    double totalxdist;
    double totalydist;
    double totalzdist;
    double avgxdistance;
    double avgydistance;
    double avgzdistance;
    double avgxdistold;
    double avgydistold;
    double avgzdistold;
    double xdistideal;
    double ydistideal;
    double zdistideal;
    double dx_new;
    double dy_new;
    double dz_new;
    double ddxvel;
    double ddyvel;
    double ddzvel;
    double xvelsum;
    double yvelsum;
    double zvelsum;
    double totalxvel;
    double totalyvel;
    double totalzvel;
    double avgxvel;
    double avgyvel;
    double avgzvel;
    double avgxvelold;
    double avgyvelold;
    double avgzvelold;
    double velwghtsum;
    double totalvelwght;
    double bubvolsum;
    double totbubvol;
    double denssum;
    double totbubdens;
    double xcforcesum;
    double ycforcesum;
    double zcforcesum;
    double totalxcforce;
    double totalycforce;
    double totalzcforce;
    double totalxcforceold;
    double totalycforceold;
    double totalzcforceold;
    double avgxcforce;
    double avgycforce;
    double avgzcforce;
    double avgxcforceold;
    double avgycforceold;
    double avgzcforceold;
    double avgxcf;
    double avgycf;
    double avgzcf;
    double x_c_f;
    double y_c_f;
    double z_c_f;
    double xforcenewtsum;
    double yforcenewtsum;
    double zforcenewtsum;
    double totxfnewtsum;
    double totyfnewtsum;
    double totzfnewtsum;
    double xcfnewtons;
    double ycfnewtons;
    double zcfnewtons;
    double rholiq;
    double rhogas;
    double coalbubrad;
    double avgxcoordold[100];
    double avgycoordold[100];
    double avgzcoordold[100];
    int i_res_cf;
    int nzinBsum;
    int ntotnzinB;
    int iLSet;
    int iuse_vfcont_cap;
    int i_num_bubbles;
    int ivconstraint;
    int iSolvLSSclr1;
    int iSolvLSSclr2;
    int i_redist_loop_flag;
    int i_redist_max_iter;
    int i_spat_var_eps_flag;
    int i_dtlset_cfl;
    int i_check_prox;
    int i_gradphi;
    int i_focusredist;
    int i_AdjRedistVel;
    int buintcfl;
    int iBT;
    int iFT;
    int icoalCtrl;
    int coalcon;
    int update_coalcon;
    int coaltimtrak;
    int coalest;
    int coalcon_rem[100];
  } levlset;

  extern struct {
    double xcfcoeff[10];
    double ycfcoeff[9];
    double zcfcoeff[9];
    double DomainSize[6];
    double shear_rate;
    double vel_centre;
    double y_drag_flip;
    int icforz;
    int icforz_where;
    int numts_histyavg;
    int iClrLiq;
    int iBK;
    int Nbubtot;
    int Nghost;
  } bubstudy;

  extern struct { 
    int nshape;
    int nshapeb;
    int maxshb;
    int nshl;
    int nshlb;
    int nfath;
    int ntopsh;
    int nsonmax;
  } shpdat;

  extern struct { 
    int mshp;
    int mshgl;
    int mwght;
    int mshpb;
    int mshglb;
    int mwghtb;
    int mmut;
    int mrhot;
    int mxst;
  } datpnt;

  extern struct { 
    int lelCat;
    int lcsyst;
    int iorder;
    int nenb;
    int nelblk;
    int nelblb;
    int ndofl;
    int nsymdl;
    int nenl;
    int nfacel;
    int nenbl;
    int intind;
    int mattyp;
  } elmpar;

  extern struct {
    int numerr;
  } errpar;

  extern struct { 
    double E3nsd;
    int I3nsd;
    int nsymdf;
    int ndofBC;
    int ndiBCB;
    int ndBCB;
    int Jactyp;
    int jump;
    int ires;
    int iprec;
    int iprev;
    int ibound;
    int idiff;
    int lhs;
    int itau;
    int ipord;
    int ipred;
    int lstres;
    int iepstm;
    double dtsfct;
    double dtsfctsclr;
    double taucfct;
    int ibksiz;
    int iabc;
    int isurf;
    int idflx;
    double Bo;
    double CoalInvSigma;
    double presavg;
    int EntropyPressure;
  } genpar;

  extern struct { 
    double epstol[8];  /* 1+ max number of scalars  (beginning of the end of time sequences) */
    double Delt[MAXTS];
    double CFLfl[MAXTS];
    double CFLsl[MAXTS];
    int nstep[MAXTS];
    int niter[MAXTS];
    int impl[MAXTS];
    double rhoinf[MAXTS];
    int LHSupd[6];
    int loctim[MAXTS];
    double deltol[2][MAXTS];
    double CFLfl_max;
    int iCFLfl_maxelem;
    int iflag_cfl_dt;
    double CFL_limit[2];
    double timestart; 
    double CFLls_max;
    int iCFLls_maxelem;
    int svLSFlag;
    double CFLbuint_max;
    double factor_buint;
    double factor_CFLfl;
  } inpdat;

  extern struct { 
    int iin;
    int igeom;
    int ipar;
    int ibndc;
    int imat;
    int iecho;
    int iout;
    int ichmou;
    int irstin;
    int irstou;
    int ihist;
    int iflux;
    int ierror;
    int itable;
    int iforce;
    int igraph;
    int itime;
    int ivol;
    int istat;
    int ivhist;
  } mio;

  extern struct { 
    double fin;
    double fgeom;
    double fpar;
    double fbndc;
    double fmat;
    double fecho;
    double frstin;
    double frstou;
    double fhist;
    double ferror;
    double ftable;
    double fforce;
    double fgraph;
    double ftime;
    double fvol;
    double fstat;
    double fvhist;
  } mioname;

  extern struct { 
    double eGMRES;
    int lGMRES;
    int iKs;
    int ntotGM;
  } itrpar;

  extern struct { 
    int mHBrg;
    int meBrg;
    int myBrg;
    int mRcos;
    int mRsin;
  } itrpnt;

  extern struct { 
    double datmat[MAXTS][7][3];
    int matflg[MAXTS][6];
    int nummat;
    int mexist;
  } matdat;

  extern struct {
    double tmu;
    double trho;
    double omu;
    double orho;
    double qrts0;
    double qrts1;
    int iramp;
    int nrts0;
    int nrts1;
  } matramp;

  extern struct { 
    double pr, Planck, Stephan, Nh, Rh, Rgas;
    double gamma, gamma1, s0, xO2, xN2, Msh[5], h0sh[5], Trot[5], sigs[5];
	double Tvib[5], g0s[5], dofs[5], Rs[5], h0s[5], cpsh[5], cps[5], cvs[5];
	double s0sh[5];
    //, const;
    //double yN2,    yO2,
    //double   cps[5], ;
    //,ithm;
  } mmatpar;

  extern struct { 
    double ro;
    double vel;
    double temper;
    double press;
    double entrop;
    int ntout;
    int ioform;
    int iowflux;
    int iofieldv;
    char iotype[80];
    int ioybar;
  } outpar;

  extern struct { 
    int mbeg;
    int mend;
    int mprec;
  } point;

  extern struct { 
    double epsM;
    int iabres;
  } precis;

  extern struct { 
    int npro;
  } propar;

  extern struct { 
    double resfrt;
  } resdat;

  extern struct { 
    int imap;
    int ivart;
    int iDC;
    int iPcond;
    int Kspace;
    int nGMRES;
    int iconvflow;
    int iconvsclr;
    int idcsclr[2];
  } solpar;

  extern struct { 
    double time;
    double CFLfld;
    double CFLsld;
    double Dtgl;
    double Dtmax;
    double alpha;
    double etol;
    int lstep;
    int ifunc;
    int itseq;
    int istep;
    int iter;
    int nitr;
    double almi;
    double alfi;
    double gami;
    double flmpl;
    double flmpr;
    double dtol[2];
    int iCFLworst;
  } timdat;

  extern struct { 
    int LCtime;
    int ntseq;
  } timpar;

  extern struct { 
    int numeqns[100];
    int minIters;
    int maxIters;
    int iprjFlag;
    int nPrjs;
    int ipresPrjFlag;
    int nPresPrjs;
    double prestol;
    double statsflow[6];
    double statssclr[6];
    int iverbose;
  } incomp;

  extern struct { 
    double ccode[13];
  } mtimer1;

  extern struct { 
    double flops;
    double gbytes;
    double sbytes;
    int iclock;
    int icd;
    int icode;
    int icode2;
    int icode3;
  } mtimer2;

  extern struct { 
    double cpu[11];
    double cpu0[11];
    int nacess[11];
  } timer3;

  extern struct { 
    double title;
    int ititle;
  } title;

  extern struct {
    int intg[MAXTS][2];
	int intpt[3];
	int intptb[3];
  } intdat;

  extern struct {
    double shvebct;
    double bcttimescale;
    double ValueListResist[MAXSURF+1];
    double rhovw;
    double thicknessvw;
    double evw;
    double rnuvw;
    double rshearconstantvw;
    double betai;
    int icardio;
    int itvn;
    int ipvsq;
    int numResistSrfs;
    int nsrflistResist[MAXSURF+1];
    int numImpSrfs;
    int nsrflistImp[MAXSURF+1];
    int impfile;
    int ideformwall;  
    int iwallmassfactor;
    int iwallstiffactor;      
    int tvbcswitch;
    int ibcb_conv_p;
    int ibcb_conv_p_norm;
    int itvbc;
    int iRBCT;
 } nomodule;

  extern struct {
    double bubboil;
    double solheat;
    double bubgrow;
    double h_fg;
    double T_sat;
    double epsilon_lst;
    double breakup2w;
    double breaksclr;
    double delt_T[100];
    double numshell_out[100];
    double numshell_in[100];
    double bubble_vol[100];
    double bubble_tempG[100];
    double bubbleID[100]; 
    double numshell_out2[100];  
  } pcboiling;  

  extern struct {
    double CA_flag;
    double theta_adv;
    double theta_rec;
    double Forcecont;
    double stretch;
    double Fapp_thick;
    double Fapp_heigh;
    double CAramp;
    double CAF_upper;
    double CAF_interval;
    double CA_startT;
    double bubble_avgCA;
  } contactangle;

  extern struct {
    int seqsize;
    int stepseq[100];
  } sequence;

  extern struct {
    double strong_eps;      /* strong criterion Stuben factor    */
    double ramg_eps;        /* AMG convergence eps               */
    double ramg_relax;       /* relaxation factor Gauss-Seidel/Jac*/
    double ramg_trunc;      /* truncation select */
 } amgvarr;
  
  extern struct {
    int irun_amg;           /* Employ AMG feature solfar.f      */
    int irun_amg_sa;        /* Run AMG stand alone ,tamg or SAMG */
    int irun_amg_prec;      /* Run AMG as preconditioner to CG */
    int iamg_verb;          /* amg verbosity flag                */
    int iamg_neg_sten;      /* neg only stencil or neg and pos   */
    int iamg_nlevel;        /* number of levels 2-V etc.         */
    int iamg_c_solver;     /* solve fine level iter. method     */
    int iamg_prescale;       /* diagonal scale AMG             */
    int iamg_init;           /* setup flag */
    int iamg_ppe_frez;       /* how many solfars to re extract ppe */
    int iamg_setup_frez;    /* how many solfars to re setup amg */
    int iamg_iai_frez;      /* how many solfars to re iai amg */
    int iamg_interp;        /* interpolation select */
    int maxnev;             /* total eigenvectors used for ggb*/
    int maxncv;             /* total iterative vectors for ggb*/
    int mlsdeg;             /* Polynomial Smoothing (MLS) degree */
    int iamg_scale;          /* control the scaling of PPE */
    int iamg_reduce;        /* Run a reduced case */
 } amgvari;

#ifdef __cplusplus
}
#endif
