//This version computes mean time and fluctuating values, budget of t' and v't'//
//and a contour file of the vorticity//
#include <iostream>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "phastaIO.h"
#include <fstream>

/* #include <map>
#define for if (0) {} else for
#ifndef WIN32
#include <unistd.h>
#include <strings.h>
#else
#include <direct.h>
void  bzero(void* ptr, size_t sz) {
  int i;
  char *cptr;
  cptr = (char*) ptr;
  for (i=0; i < sz; i++) {
    cptr[i]=0;
  }
  return;
}
#define ijktran_ IJKTRAN
#endif */


using namespace std;
//Routine for ijk transform

/* extern "C" {
void ijktran_(int, double*, double*, double*, double*, double*);
}*/

void read_real_array_(int*, double*, int*);
void read_int_array_(int*, int*, int*);
void read_header_(int*,const char*,int*,int*,int*);
void write_restart_(int* , int* , int* , int* , double* , double*);
void write_timeavg_(int* , int* , int* , int*, int* , int* , double*);
extern FILE* frest;
extern FILE* fgeombc;

int geometry_and_connectivity(int* array, double* &xglobal, int* &ien, 
                               int** &ncorp2d, bool RequestedVolCheck);
/* This routine assembles the global array from the local (processor */
/* specific) array                                                   */
void reduce(int numvar, int lsize, int gsize, int proc,
	    double* local, double* global, int** ncorp){
    int i,j;
/*    double qmax[20], qmin[20];
    for(i=0; i< 20; i++){
        qmax[i]=-100000000;
        qmin[i]= 100000000;
    } */
    for(i=0; i< numvar; i++){
        for(j=0; j< lsize ; j++){ 
            global[i*gsize+ncorp[proc][j]-1] = local[i*lsize+j];
//            if(qmax[i] < local[i*lsize+j]) qmax[i]=local[i*lsize+j];
//            if(qmin[i] > local[i*lsize+j]) qmin[i]=local[i*lsize+j];
        }
    }
/*    for(i=0; i< numvar; i++)
        printf("var.=%d, min=%f, max=%f \n", i,qmin[i],qmax[i]); */
    return;
}


//   nv is the number of degrees of freedom per node
//   numel is the total number of elements
//   ien is the TRUE global connectivity  
//       (index local element node number and element
//        number to TRUE global node number)
//   numnp is the total number of nodes                          
//   x is the coordinates (3*numnp)
//   q is the solution (nv*numnp)
//   nen is the number of nodes in an element


int main(int argc, char* argv[])
{
  int *nshg, *nnod, numvar, nsd, nshgtot;
  int numprocs,i,j,k,l, lstep;
  double **xyz;
  double **qstep;
  double **qavrg;
  double *qglobal, *xglobal;
  int **ncorp2d, *ien;
  int iarray[20];
  char rfname[40];
  char gfname[40];
  char fname1[255];
  int intfromfile[50];
  int ione=1, itwo=2, ithree=3,  iseven=7;
  int igeom, irstin;
  int ierr;
  int ixsiz, iqsiz;
  int clread=0, step1=0;
  int size, nitems;
  char* iotype;
  char rfile[255];
  int magic_number = 362436;
  int* mptr = &magic_number;
  // variable used for budget calculation
  double pp,up,vp,wp,tp,*dt2,*pdtdy,*tdpdy,*mdx1,*mdx2,*mdy1,*mdy2,*qglobal1,**qinst,*qiglo;
  double *xc, *qs,x1,x2,x3,dtdx,dtdy,dtdz,dvdx,dvdy,dvdz,*disvtx,*disvty,*disvtz,*vor,dpdy;
  int *invmap, iorig,ii;
  int Nx,Ny,Nz,nnp,numvar1,nendx,neltot;
  double Pr=0.71; //molecular Prandtl number
  FILE* frest1;

  // variables used in processing the commandline
  int iarg, arglength;
  string tmpstr;

  int startstep, skipstep, stopstep;
  bool StepNumberAvailable;
  int indxu1, indxu2, indxx1, indxx2;
  double theta, cost, sint, rat, ur, ut;
  int nsamples, numvarTimeAvg;
  double nsamplesInverse;

  igeom=1; irstin=2;
  iotype = "binary";

  /* BEGIN PROCESSING COMMAND-LINE ARGUMENTS */
  /* Assume the command line is okay */
  bool BogusCmdLine = false;
  /* Assume no options specified at command line */
  bool RequestedHelp = false;
  bool RequestedTimeAvg = false;
  bool RequestedCylindrical = false;
  bool VolCheck = false;
  
  /* argc is the number of strings on the command-line */
  /*  starting with the program name */
  for(iarg=1; iarg<argc; iarg++){
    arglength = strlen(argv[iarg]);
    /* replace 0..arglength-1 with argv[iarg] */
    tmpstr.replace(0,arglength,argv[iarg],0,arglength);
    if(tmpstr=="-h"){
      RequestedHelp = true;
      cout << endl;
      cout << "usage:" <<endl;
      cout << "  Reduce -np <nParts> -ta <start> <skip> <stop> [-cyl]" << endl;
      cout << endl;
    }
    else if(tmpstr=="-np"){
      iarg++;
      if(iarg>=argc){
        cout << "Numparts not given, so exiting"<<endl;
        return(0);
      }
      numprocs = atoi(argv[iarg]);
    }
    else if(tmpstr=="-ta"){
      RequestedTimeAvg = true;
      iarg++;
      if(iarg>=argc){
        cout << "Startstep not given, so exiting"<<endl;
        return(0);
      }
      startstep = atoi(argv[iarg]);
      iarg++;
      if(iarg>=argc){
        cout << "Skipstep not given, so exiting"<<endl;
        return(0);
      }
      skipstep = atoi(argv[iarg]);
      iarg++;
      if(iarg>=argc){
        cout << "Stopstep not given, so exiting"<<endl;
        return(0);
      }
      stopstep = atoi(argv[iarg]);
      StepNumberAvailable=true;
    } 
    else if(tmpstr=="-cyl"){
      RequestedCylindrical = true;
    }
    else {
      BogusCmdLine = true;
    }
    /* reset tmpstr for next argument */
    tmpstr.erase(0,arglength);
  }
  
  /*In the case of a bogus command line, print proper usage and exit */
  if(BogusCmdLine){
    cout << endl;
    cout << "usage:" <<endl;
    cout << "  Reduce -np <nParts> -start <startstep> -skip <skipstep> -stop <stopstep> [-cyl]" << endl;
    cout << endl;
    return(0);
  }

  cout << "Will gather time-averaged fields, sampling every " << skipstep;
  cout <<  " steps, starting at step " << startstep << " and ending at step ";
  cout << stopstep << endl;
  if(RequestedCylindrical){
    cout << "Will transform cartesian velocity components to cylindrical components." << endl;
  }
  // Keep these last
  if(!StepNumberAvailable){
    cout << "No step number or range of steps given, so exiting." << endl;
    return(0);
  }
  if(RequestedHelp){
    cout << endl;
    cout << "Exiting before performing any of the above operations due to -h flag";
    cout << endl;
    return(0);
  }
  /* END PROCESSING COMMAND-LINE ARGUMENTS */
  
  
  nshg = new int[numprocs];
  nnod = new int[numprocs];
  
  if(RequestedCylindrical){
    // Read coordinates
    xyz = new double*[numprocs];
    for(i=0;i<numprocs;i++){
      bzero((void*)gfname,40);
      sprintf(gfname,"geombc.dat.%d",i+1);
      cout << "Opening " << gfname << " to read coordinate array..."<<endl;
    //
      fgeombc = fopen(gfname,"r");
      if(!fgeombc){
        cout << "Could not successfully open " << gfname << ", so exiting"<< endl;
        return(0);
      }
      fclose(fgeombc);
    //
      fgeombc = fopen(gfname,"r");
      read_header_(&igeom,"co-ordinates",intfromfile,&itwo,&ierr);
      nnod[i]=intfromfile[0];
      nsd=intfromfile[1];
      ixsiz=nnod[i]*nsd;
      xyz[i]=new double[ixsiz];
      read_real_array_(&igeom,xyz[i],&ixsiz);
      fclose(frest);
    }
  }  

  // Scan for solution field dimensions
  for(i=0;i<numprocs;i++){
    bzero((void*)rfname,40);
    sprintf(rfname,"restart.%d.%d",startstep,i+1);
    cout << "Opening " << rfname << " to scan sizes..."<<endl;
    //
    frest = fopen(rfname,"r");
    if(!frest){
      cout << "Could not successfully open " << rfname << ", so exiting"<< endl;
      return(0);
    }
    fclose(frest);
    //
    frest = fopen(rfname,"r");
    read_header_(&irstin,"solution",intfromfile,&ithree,&ierr);
    nshg[i]=intfromfile[0];
    numvar=intfromfile[1];
    fclose(frest);
  }  

  // We get the time-average of 5 flow variables, 5 squares of
  // these variables, and 6 cross-terms from the reynolds stress
  // tensor. Also, t'2 u',t'2 v', t'v'u', t'v'2, p't'. */
  numvarTimeAvg = 21;
  //number of fluctuating components passed:t', p', v',u',w' and instantaneous values V,U,W
  int nfc=8;
  qavrg = new double*[numprocs];
  qstep = new double*[numprocs];
  //For instantaneous parameters
  qinst = new double*[numprocs];

  for(i=0;i<numprocs;i++){
    qavrg[i]=new double[nshg[i]*numvarTimeAvg];
    qinst[i]=new double[nshg[i]*nfc];
    for(j=0;j<nshg[i]*numvarTimeAvg;j++) qavrg[i][j]=0.0;
    for(j=0;j<nshg[i]*nfc;j++) qinst[i][j]=0.0;
    qstep[i]=new double[nshg[i]*numvar];
  }

  nsamples=0;

  // Loop over time steps
  int stepnumber;
  for(stepnumber=startstep; stepnumber<=stopstep; stepnumber+=skipstep){
    cout<<endl<<"Reading step "<<stepnumber<<":";
    for(i=0; i<numprocs; i++){
      // Open this partition's solution file
      sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
      cerr<<" "<<i+1;
      frest = fopen(rfname,"r");
      if(!frest){
        cout << "Could not open "<<rfname<<" so exiting"<<endl;
        return(0);
      }
      // Read solution
      read_header_(&irstin,"solution",intfromfile,&ithree,&ierr);
      iqsiz=nshg[i]*numvar;
      read_real_array_(&irstin,qstep[i],&iqsiz);
      // Close this partition's solution file
      fclose(frest);
      
      if(RequestedCylindrical){
        for(k=0;k<nnod[i]; k++){
          indxx1=k;
          indxx2=k+nnod[i];
          indxu1=k+nnod[i];
          indxu2=k+2*nnod[i];
          theta=0;
          if(xyz[i][indxx2]>0) {
            rat=-xyz[i][indxx1]/xyz[i][indxx2];
            theta=atan(rat);
          }
          sint=sin(theta);
          cost=cos(theta);
          
          ur=-sint*qstep[i][indxu1]+cost*qstep[i][indxu2];
          ut=-cost*qstep[i][indxu1]-sint*qstep[i][indxu2];
          qstep[i][indxu1]=ur;
          qstep[i][indxu2]=ut;
        }
      }
      //MEAN TIME VALUES
      for(k=0; k< nshg[i]; k++){
       //P
        qavrg[i][0*nshg[i]+k] = qavrg[i][0*nshg[i]+k] + qstep[i][0*nshg[i]+k];
       //U
        qavrg[i][1*nshg[i]+k] = qavrg[i][1*nshg[i]+k] + qstep[i][1*nshg[i]+k];
       //V
        qavrg[i][2*nshg[i]+k] = qavrg[i][2*nshg[i]+k] + qstep[i][2*nshg[i]+k];
       //W
        qavrg[i][3*nshg[i]+k] = qavrg[i][3*nshg[i]+k] + qstep[i][3*nshg[i]+k];
       //T
        qavrg[i][4*nshg[i]+k] = qavrg[i][4*nshg[i]+k] + qstep[i][4*nshg[i]+k];
      }
    } //End loop over partitions
    nsamples = nsamples + 1;
  } // End loop over time steps

  // Complete the time-averaging process
  nsamplesInverse = 1.0/nsamples;
  for(i=0;i<numprocs;i++){
    for(k=0;k<5*nshg[i];k++){
      qavrg[i][k]=qavrg[i][k]*nsamplesInverse;
    }
  }

//    write_timeavg_(&startstep, &skipstep, &stopstep, &i,&nshg[i], &numvarTimeAvg, qavrg[i]);
  
  cout <<endl;
  cout << "END CALCULATION MEAN VALUES"<<endl;
  cout <<endl;
  cout << "START OF  RMS CALCULATION" <<endl;
  // Scan for solution field dimensions
  for(i=0;i<numprocs;i++){
    bzero((void*)rfname,40);
    sprintf(rfname,"restart.%d.%d",startstep,i+1);
    cout << "Opening " << rfname << " to scan sizes..."<<endl;
    //
    frest = fopen(rfname,"r");
    if(!frest){
      cout << "Could not successfully open " << rfname << ", so exiting"<< endl;
      return(0);
    }
    fclose(frest);
    //
    frest = fopen(rfname,"r");
    read_header_(&irstin,"solution",intfromfile,&ithree,&ierr);
    nshg[i]=intfromfile[0];
    numvar=intfromfile[1];
    fclose(frest);
  }  

  nsamples=0;

  // Loop over time steps
  for(stepnumber=startstep; stepnumber<=stopstep; stepnumber+=skipstep){
    cout<<endl<<"Reading step "<<stepnumber<<":";
    for(i=0; i<numprocs; i++){
      // Open this partition's solution file
      sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
      cerr<<" "<<i+1;
      frest = fopen(rfname,"r");
      if(!frest){
        cout << "Could not open "<<rfname<<" so exiting"<<endl;
        return(0);
      }
      // Read solution
      read_header_(&irstin,"solution",intfromfile,&ithree,&ierr);
      iqsiz=nshg[i]*numvar;
      read_real_array_(&irstin,qstep[i],&iqsiz);
      // Close this partition's solution file
      fclose(frest);
      
      if(RequestedCylindrical){
        for(k=0;k<nnod[i]; k++){
          indxx1=k;
          indxx2=k+nnod[i];
          indxu1=k+nnod[i];
          indxu2=k+2*nnod[i];
          theta=0;
          if(xyz[i][indxx2]>0) {
            rat=-xyz[i][indxx1]/xyz[i][indxx2];
            theta=atan(rat);
          }
          sint=sin(theta);
          cost=cos(theta);
          
          ur=-sint*qstep[i][indxu1]+cost*qstep[i][indxu2];
          ut=-cost*qstep[i][indxu1]-sint*qstep[i][indxu2];
          qstep[i][indxu1]=ur;
          qstep[i][indxu2]=ut;
        }
      }

      for(k=0; k< nshg[i]; k++){
        //p'2
        pp=qstep[i][0*nshg[i]+k]-qavrg[i][0*nshg[i]+k];        
        qavrg[i][5*nshg[i]+k] += pow(pp,2);
        //u'2
        up=qstep[i][1*nshg[i]+k]-qavrg[i][1*nshg[i]+k];
        qavrg[i][6*nshg[i]+k] += pow(up,2);
        //v'2
        vp=qstep[i][2*nshg[i]+k]-qavrg[i][2*nshg[i]+k];
        qavrg[i][7*nshg[i]+k] += pow(vp,2);
        //w'2
        wp=qstep[i][3*nshg[i]+k]-qavrg[i][3*nshg[i]+k];
        qavrg[i][8*nshg[i]+k] += pow(wp,2);
        //t'2
        tp=qstep[i][4*nshg[i]+k]-qavrg[i][4*nshg[i]+k];
        qavrg[i][9*nshg[i]+k] += pow(tp,2);
        
        //u'v'
        qavrg[i][10*nshg[i]+k] += up*vp;
        //u'w'
        qavrg[i][11*nshg[i]+k] += up*wp;
        //v'w'
        qavrg[i][12*nshg[i]+k] += vp*wp;
        //u't'
        qavrg[i][13*nshg[i]+k] += up*tp;
        //v't'
        qavrg[i][14*nshg[i]+k] += vp*tp;
        //w't'
        qavrg[i][15*nshg[i]+k] += wp*tp;

	//Budget of temperature variance, t'2/2
	//t'2 u'
        qavrg[i][16*nshg[i]+k] += pow(tp,2)*up;
        //t'2 v'
        qavrg[i][17*nshg[i]+k] += pow(tp,2)*vp;
	//Budget of wall-normal turbulent heat fluxes
	//t'v'u'
        qavrg[i][18*nshg[i]+k] += tp*vp*up;
	//t'v'2
        qavrg[i][19*nshg[i]+k] += tp*pow(vp,2);
	//p't'
        qavrg[i][20*nshg[i]+k] += pp*tp;
	//instantaneous t'
        qinst[i][k] = tp;
	//instantaneous p'
        qinst[i][1*nshg[i]+k] = pp;
	//instantaneous v'
        qinst[i][2*nshg[i]+k] = vp;
        //instantaneous u'
        qinst[i][3*nshg[i]+k] = up;
        //instantaneous w'
        qinst[i][4*nshg[i]+k] = wp;
        //instantaneous V
        qinst[i][5*nshg[i]+k] = qstep[i][2*nshg[i]+k];
        //instantaneous U
        qinst[i][6*nshg[i]+k] = qstep[i][1*nshg[i]+k];
        //instantaneous W
        qinst[i][7*nshg[i]+k] = qstep[i][3*nshg[i]+k];
      }
    } //End loop over partitions
    nsamples = nsamples + 1;
    //Starting to transform parameters into structured arrays i,j,k
    /* Calculates global nodal coordinates, ien and ncorp2d arrays*/
    if (stepnumber == startstep){
    clread = geometry_and_connectivity(iarray, xglobal, ien, ncorp2d, VolCheck);
    if (clread == 1) return(0);
    nshgtot  = iarray[1];
    qiglo = (double *) malloc( nfc*nshgtot * sizeof(double));
    }
    /* map solution to global */
    for(i=0; i<numprocs; i++){
      reduce(nfc, nshg[i], nshgtot, i, qinst[i], qiglo, ncorp2d);
    }

    //ifstream geom("geom.old");
    //geom >> nnp >> nsd;
    
    //xc = (double*) malloc(nnp*3 * sizeof(double));
    //qs = (double*) malloc(nnp*nfc * sizeof(double));
    //invmap = (int*) malloc( nnp * sizeof(int));
    
    if (stepnumber == startstep){
      //Reading Nx, Ny and Nz
      ifstream nxnynz("nxnynz.dat");
      if (nxnynz){
	cout<<"nxnynz.dat exists"<<endl;
        nxnynz >> Nx >> Ny >> Nz >> neltot;
        nxnynz.close();
      }
      else{
	cout<<"File nxnynz.dat does not exist, so exiting"<<endl;
        return(1);
      }
    
      //Reading coordinates
    ifstream geom("geom.old");
    geom >> nnp >> nsd;
    xc = (double*) malloc(nnp*3 * sizeof(double));
    qs = (double*) malloc(nnp*nfc * sizeof(double));
    invmap = (int*) malloc( nnp * sizeof(int));
    dt2 = (double*) malloc(nnp * sizeof(double));
    pdtdy = (double*) malloc(nnp * sizeof(double));
    tdpdy = (double*) malloc(nnp * sizeof(double));
    mdx1 = (double*) malloc(nnp * sizeof(double));
    mdx2 = (double*) malloc(nnp * sizeof(double));
    mdy1 = (double*) malloc(nnp * sizeof(double));
    mdy2 = (double*) malloc(nnp * sizeof(double));
    disvtx = (double*) malloc(nnp * sizeof(double));
    disvty = (double*) malloc(nnp * sizeof(double));
    disvtz = (double*) malloc(nnp * sizeof(double));
    vor = (double*) malloc(nnp*8 * sizeof(double));

    for (i=0; i<nnp;i++){
    dt2[i]=0.0;
    pdtdy[i]=0.0;
    tdpdy[i]=0.0;
    mdx1[i]=0.0;
    mdx2[i]=0.0;
    mdy1[i]=0.0;
    mdy2[i]=0.0;
    disvtx[i]=0.0;
    disvty[i]=0.0;
    disvtz[i]=0.0;
      }
    cout<<"vectors allocated"<<endl;
    //}

    for (i=0; i<nnp;i++){
      geom>>j>>j>>x1>>x2>>x3;
      xc[i]=x1;
      xc[i+nnp]=x2;
      xc[i+2*nnp]=x3;
    }
    geom.close();
    ifstream geomk("geom.kmap");
    for (i=0; i<nnp;i++){
      geomk>>iorig;
      //invmap[iorig+1]=i;
      invmap[iorig]=i;
    }
    geomk.close();
    }//done in the first stepnumber
    //if (stepnumber != startstep){
    //qs = (double*) malloc(nnp*nfc * sizeof(double));
    //}

    //Transforming to structured grid
    for (j=0; j<nfc;j++){
      for (i=0; i<nnp;i++){
	ii=invmap[i];
       qs[j*nnp+i] = qiglo[j*nnp+ii];
      }
    }
    
    // Dissipation for t'2 and v't'
    for (i=0; i<(nnp-(Nz*Ny));i++){
      dtdx=(qs[i+(Nz*Ny)]-qs[i])/(xc[i+(Nz*Ny)]-xc[i]);
      dtdy=(qs[i+(Nz)]-qs[i])/(xc[i+Nz+nnp]-xc[i+nnp]);
      dtdz=(qs[i+1]-qs[i])/(xc[i+1+2*nnp]-xc[i+2*nnp]);
      dt2[i] += dtdx*dtdx+dtdy*dtdy+dtdz*dtdz;
      dvdx=(qs[i+(Nz*Ny)+2*nnp]-qs[i+2*nnp])/(xc[i+(Nz*Ny)]-xc[i]);
      dvdy=(qs[i+(Nz)+2*nnp]-qs[i+2*nnp])/(xc[i+Nz+nnp]-xc[i+nnp]);
      dvdz=(qs[i+1+2*nnp]-qs[i+2*nnp])/(xc[i+1+2*nnp]-xc[i+2*nnp]);
      disvtx[i] += dtdx*dvdx;
      disvty[i] += dtdy*dvdy;
      disvtz[i] += dtdz*dvdz;
	      }
    //product p'dt/dy (pressure-temp. gradient correlation)
    for (i=0; i<(nnp-(Nz*Ny));i++){
      dtdy=(qs[i+(Nz)]-qs[i])/(xc[i+Nz+nnp]-xc[i+nnp]);
      pdtdy[i] += qs[i+1*nnp]*dtdy;
    }
    //product t'dp/dy (temp. pressure gradient correlation)
    for (i=0; i<(nnp-(Nz*Ny));i++){
      dpdy=(qs[i+(Nz)+1*nnp]-qs[i+1*nnp])/(xc[i+Nz+nnp]-xc[i+nnp]);
      tdpdy[i] += qs[i]*dpdy;
    }
    //product mdx=(1/Pr v'dt/dx + t'dv/dx) and mdy=(1/Pr v'dt/dy+t'dv/dy) for molecular diff. of v't'
    for (i=0; i<(nnp-(Nz*Ny));i++){
      dtdx=(qs[i+(Nz*Ny)]-qs[i])/(xc[i+(Nz*Ny)]-xc[i]);
      dvdx=(qs[i+(Nz*Ny)+2*nnp]-qs[i+2*nnp])/(xc[i+(Nz*Ny)]-xc[i]);
      mdx1[i] += qs[i+2*nnp]*dtdx;
      mdx2[i] += qs[i]*dvdx;
      dtdy=(qs[i+(Nz)]-qs[i])/(xc[i+Nz+nnp]-xc[i+nnp]);
      dvdy=(qs[i+(Nz)+2*nnp]-qs[i+2*nnp])/(xc[i+Nz+nnp]-xc[i+nnp]);
      mdy1[i] += qs[i+2*nnp]*dtdy;
      mdy2[i] += qs[i]*dvdy;
	}

    //cout<<"dt2[100] "<<dt2[100]<<" "<<xc[100]<<" "<<xc[100+nnp]<<endl;
    //cout<<"test "<<qs[1000+(Nz)]<<" "<<qs[1000]<<endl;
    //cout<<"test "<<xc[1000+2*(Nz)]<<" "<<xc[1000+(Nz)]<<" "<<xc[1000]<<endl;
    //cout<<"test "<<xc[1000+2*(Nz)+nnp]<<" "<<xc[1000+(Nz)+nnp]<<" "<<xc[1000+nnp]<<endl;
    //cout<<"test "<<xc[102+2*nnp]<<" "<<xc[101+2*nnp]<<" "<<xc[100+2*nnp]<<endl;

    //Iso-countour of instantaneous vorticity
    if (stepnumber==stopstep){
      cout<<"computing vorticity "<<stepnumber<<endl;
      //fluglobal = (double *) malloc( 5*nshgtot * sizeof(double));
      for (i=0; i<(nnp-(Nz*Ny));i++){
	vor[i]=(qs[i+(Nz)+7*nnp]-qs[i+7*nnp])/(xc[i+Nz+nnp]-xc[i+nnp])-(qs[i+1+5*nnp]-qs[i+5*nnp])/(xc[i+1+2*nnp]-xc[i+2*nnp]);
        vor[i+1*nnp]=(qs[i+1+6*nnp]-qs[i+6*nnp])/(xc[i+1+2*nnp]-xc[i+2*nnp])-(qs[i+(Nz*Ny)+7*nnp]-qs[i+7*nnp])/(xc[i+(Nz*Ny)]-xc[i]);
        vor[i+2*nnp]=(qs[i+(Nz*Ny)+5*nnp]-qs[i+5*nnp])/(xc[i+(Nz*Ny)]-xc[i])-(qs[i+(Nz)+6*nnp]-qs[i+6*nnp])/(xc[i+Nz+nnp]-xc[i+nnp]);
        vor[i+3*nnp]=(qs[i+(Nz)+4*nnp]-qs[i+4*nnp])/(xc[i+Nz+nnp]-xc[i+nnp])-(qs[i+1+2*nnp]-qs[i+2*nnp])/(xc[i+1+2*nnp]-xc[i+2*nnp]);
        vor[i+4*nnp]=(qs[i+1+3*nnp]-qs[i+3*nnp])/(xc[i+1+2*nnp]-xc[i+2*nnp])-(qs[i+(Nz*Ny)+4*nnp]-qs[i+4*nnp])/(xc[i+(Nz*Ny)]-xc[i]);
        vor[i+5*nnp]=(qs[i+(Nz*Ny)+2*nnp]-qs[i+2*nnp])/(xc[i+(Nz*Ny)]-xc[i])-(qs[i+(Nz)+3*nnp]-qs[i+3*nnp])/(xc[i+Nz+nnp]-xc[i+nnp]);
        vor[i+6*nnp]=pow((vor[i]*vor[i]+vor[i+1*nnp]*vor[i+1*nnp]+vor[i+2*nnp]*vor[i+2*nnp]),0.5);
        vor[i+7*nnp]=pow((vor[i+3*nnp]*vor[i+3*nnp]+vor[i+4*nnp]*vor[i+4*nnp]+vor[i+5*nnp]*vor[i+5*nnp]),0.5);
      }
      //Printing iso-countour of vorticity
      frest1 = fopen("vorticity.dat","w");
      numvar1=8;
      nendx=8;
      //neltot=170079; //lowRe
      //neltot=169884;
      fprintf(frest1,"Title = \"FE-VOLUME BRICK DATA SET\"\n");
      fprintf(frest1,"VARIABLES = \"X\", \"Y\", \"Z\", \"WX\", \"WY\", \"WZ\", \"wpx\", \"wpy\", \"wpz\", \"MW\",\"mwp\"\n");
      fprintf(frest1,"ZONE N= %d E= %d  F=FEPOINT  ET=BRICK\n", nshgtot, neltot);
      for(i=0; i< nshgtot; i++){
	for(j=0; j < 3; j++)
	  fprintf(frest1,"%e  ",xglobal[j*nshgtot+i]);
	for(k=0; k < numvar1; k++)
	  fprintf(frest1,"%22.16e  ",vor[k*nshgtot+i]);
	fprintf(frest1,"\n");
      }

      fprintf(frest1,"\n");
      for(i=0; i< neltot; i++){
	for(j=0; j < nendx; j++)
	  fprintf(frest1,"%d  ",ien[j*neltot+i]);
	fprintf(frest1,"\n");
      }
      fclose(frest1);
      cout<<"iso-vorticity file printed"<<endl;

    }//end loop when last step

   
    //free(qs);
    //free(qiglo);

    //Ending the transformation for this timestep

  } // End loop over time steps

  // Complete the time-averaging process
  nsamplesInverse = 1.0/nsamples;
  for(i=0;i<numprocs;i++){
    for(k=5*nshg[i];k<numvarTimeAvg*nshg[i];k++){
      qavrg[i][k]=qavrg[i][k]*nsamplesInverse;
    }
  }

  //END OF RMS CALCULATION
  

    /* Calculates global nodal coordinates, ien and ncorp2d arrays*/
    clread = geometry_and_connectivity(iarray, xglobal, ien, ncorp2d, VolCheck);
    
    if (clread == 1) return(0);

    nshgtot  = iarray[1];

    //Modification: new parameters are added to qglobal
    
    qglobal1 = (double *) malloc( numvarTimeAvg*nshgtot * sizeof(double));
  
  /* map solution to global */    
  for(i=0; i<numprocs; i++){
    reduce(numvarTimeAvg, nshg[i], nshgtot, i, qavrg[i], qglobal1, ncorp2d);  
  }

  //Modification: numvarTimeAvg is updated
  int newvar;
  newvar=6;

  numvarTimeAvg=numvarTimeAvg+newvar;
  qglobal = (double *) malloc( numvarTimeAvg*nshgtot * sizeof(double));
  for(i=0; i<(numvarTimeAvg-newvar)*nnp; i++){
    qglobal[i]=qglobal1[i];
  }
  //free(qglobal1);
  for(i=0; i<nnp; i++){
   qglobal[(numvarTimeAvg-newvar)*nnp+i]=dt2[i]*nsamplesInverse;
   qglobal[(numvarTimeAvg-newvar+1)*nnp+i]=pdtdy[i]*nsamplesInverse;
   qglobal[(numvarTimeAvg-newvar+2)*nnp+i]=(1./Pr*mdx1[i]+mdx2[i])*nsamplesInverse;
   qglobal[(numvarTimeAvg-newvar+3)*nnp+i]=(1./Pr*mdy1[i]+mdy2[i])*nsamplesInverse;
   qglobal[(numvarTimeAvg-newvar+4)*nnp+i]=(disvtx[i]+disvty[i]+disvtz[i])*nsamplesInverse;
   qglobal[(numvarTimeAvg-newvar+5)*nnp+i]=tdpdy[i]*nsamplesInverse;
  }
  cout<<"ATTENTION! qglobal has " << numvarTimeAvg <<" variables" << endl ;
  //cout<<"qglobal "<<qglobal[(numvarTimeAvg-newvar)*nnp+100]<<" "<<dt2[100]<<" "<<nsamplesInverse<<endl;
  sprintf(rfile,"%s.%d-%d.%s","stats",startstep,stopstep,"out");
  openfile_(rfile, "write" , &irstin);

  writestring_( &irstin,"# PHASTA Input File Version 2.0\n");
  writestring_( &irstin, "# Byte Order Magic Number : 362436 \n");

  bzero( (void*)rfile, 255 );
  sprintf(rfile,"# Output generated by phPost version 2.7:  \n");
  writestring_( &irstin, rfile );
	  
  size = 1;
  nitems = 1;
  iarray[0] = 1;
  writeheader_( &irstin, "byteorder magic number ",
                  (void*)iarray, &nitems, &size, "integer", iotype );
      
  writedatablock_( &irstin, "byteorder magic number ",
                     (void*)mptr, &nitems, "integer", iotype );

  bzero( (void*)rfile, 255 );
  sprintf(rfile,"number of modes : < 0 > %d\n", nshgtot);
  writestring_( &irstin, rfile );
    
  bzero( (void*)rfile, 255 );
  sprintf(rfile,"number of variables : < 0 > %d\n", numvar);
  writestring_( &irstin, rfile );
    
  size =  numvarTimeAvg*nshgtot;
  nitems = 5;
  iarray[0] = nshgtot;
  iarray[1] = numvarTimeAvg;
  iarray[2] = startstep;
  iarray[3] = stopstep;
  iarray[4] = 1;

  writeheader_( &irstin, "statistics ",
                  (void*)iarray, &nitems, &size,"double", iotype);
      
  nitems = numvarTimeAvg*nshgtot;
  writedatablock_( &irstin, "statistics ",
                     (void*)(qglobal), &nitems, "double", iotype );

  closefile_( &irstin, "write" );

  for(i=0; i< numprocs; i++) free(ncorp2d[i]);
  free(ncorp2d);
  free(qglobal);
  free(qglobal1);
  free(xglobal);
  free(ien);
  free(xc);
  free(qs);
  free(qiglo);
  free(invmap);
  free(dt2);
  free(pdtdy);
  free(tdpdy);
  free(mdx1);
  free(mdx2);
  free(mdy1);
  free(mdy2);
  free(disvtx);
  free(disvty);
  free(disvtz);
  free(vor);
  return 0;
}
