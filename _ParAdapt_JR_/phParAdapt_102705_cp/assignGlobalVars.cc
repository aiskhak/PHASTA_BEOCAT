#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <iostream.h> 
#include <fstream.h>
#include <strings.h>

#include "Input.h"
#include "phParAdapt.h"

/* global variables */
extern int globalP;
extern int fType;
extern int rStart;
extern int per;
extern int intBC;
extern int timing;
extern int wGraph;
extern int prCd;
extern int  zScale;
int WRITEASC=0;
int ensa_dof=5;
int CUBES=0;
int rRead=5;
int NSFTAG=-1; 
int lStart=0; 
int FortFormFlag=0; 
int phastaIO=1; 
int old_format=0;

char mname[256];
char fname[256];
char gname[256];

// for adaptation 
int strategy;
int nErrorVars;
double factor;
double hmax;
double hmin;
int adaptOption;
double refthreshold;
char outputFormat[100];
int adaptFlag=0;//use adaptor first
double* wght;

int preLBforAdaptivity;
double masterProcWgt;

int timeStepNumber=0;

char* iformat = "binary";
char* oformat;

extern int justpart;
extern int parted;
extern int errorName;
// should be removed in later versions
// now still used in attachIBC.cc
// the issue is whether the option should be kept
// to be able to read in a SINGLE restart file
// or if that should be moved to the partitioner
int multipleRestarts = 0;
int numRestartsToRead = 1;

char version[256];
char meshoutname[256];
char phastaVersion[256];

int numelX,Idirection;
int BYPASS=0;
int SONFATH = 0;

void 
assignGlobalVars()
{
    std::ifstream inf("adapt.inp");
    if (!inf) {
        cerr << "Error: adapt.inp not found" << endl;
        exit(0);
    }
    Input inputObject(inf);

//******************my changes, I am reading files adapt.inp Azat 09.16.04 12.58PM
//******************Help Part***************************************************

    refthreshold  = inputObject.threshold(); // reference threshold value(fraction of GlobalMaxError) above which the edges would be marked
    nErrorVars   = inputObject.inp_nErrorVars();//Number of error indicators to be read
    double *wght1 = inputObject.weights();   // Weights of the error indicators
    wght = new double[nErrorVars];
    for(int iEvar=0; iEvar<nErrorVars; iEvar++)
      wght[iEvar] = wght1[iEvar];
    //sets the global polynomial order
    globalP = inputObject.inp_globalP();

    // time step number
    timeStepNumber = inputObject.inp_timeStepNumber();
    numelX = inputObject.inp_numelX();
    NSFTAG = inputObject.inp_NSFTAG();
    ensa_dof =inputObject.inp_ensa_dof();
    strcpy(fname,inputObject.inp_fname());
    strcpy(mname,inputObject.inp_mname());
    strcpy(gname,inputObject.inp_gname() ) ;
    Idirection = inputObject.inp_Idirection();
    BYPASS  =inputObject.inp_BYPASS();
    fType  = inputObject.inp_fType();
    zScale =  inputObject.inp_zScale();
    adaptFlag  = inputObject.inp_adaptFlag();
    errorName  = inputObject.inp_errorName();
    SONFATH  = inputObject.inp_SONFATH();
    lStart  = inputObject.inp_lStart();
    rRead   = inputObject.inp_rRead();
    rStart   = inputObject.inp_rStart();
    strategy   = inputObject.inp_strategy();
    factor   = inputObject.inp_factor();

    nErrorVars   = inputObject.inp_nErrorVars();
    hmax   = inputObject.inp_hmax();
    hmin   = inputObject.inp_hmin();
    adaptOption = inputObject.inp_adaptOption ();
    multipleRestarts   = inputObject.inp_multipleRestarts();
    per   = inputObject.inp_per();
    prCd   = inputObject.inp_prCD();
    timing   = inputObject.inp_timing();
    wGraph   = inputObject.inp_wGraph();
    strcpy(phastaVersion , inputObject.inp_phastaVersion());
    old_format   = inputObject.inp_old_format();
    FortFormFlag   = inputObject.inp_FortFormFlag();
    strcpy(outputFormat ,inputObject.inp_outputFormat());

    oformat = outputFormat; 

    CUBES   = inputObject.inp_CUBES();
    intBC   = inputObject.inp_intBC();
    strcpy(version , inputObject.inp_version());
    WRITEASC   = inputObject.inp_WRITEASC();
    phastaIO   = inputObject.inp_phastaIO();
    parted     =  inputObject. inp_parted();

    preLBforAdaptivity =  inputObject. inp_preLBforAdaptivity();
    masterProcWgt         =  inputObject. inp_masterProcWgt();
}
