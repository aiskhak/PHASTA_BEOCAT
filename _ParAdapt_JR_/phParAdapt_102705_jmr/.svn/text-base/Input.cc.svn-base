#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <iostream>

#include "Input.h"

#include "SimPMesh.h"

using  std::ifstream;
using std::cerr;
using std::cout;
using std::endl;

const int lineSize = 1024;

Input::Input(ifstream &inf)
{
  char *inBuf = new char[lineSize];

  // set default values
  dnErrorVars  = 0;
  //now ref_threshhold is defining the fraction of the GlobalMaxError
  //default value is the 0.5 and all nodes, where error value >=0.5*GlobalMaxError will
  //be marked
  d_threshold = 0.5;

  dpreLBforAdaptivity = 0;
  dmasterProcWgt = 0.0;

  // parse the file
  while (inf.getline(inBuf, lineSize,' ')){
    
    // comment line
    if (inBuf[0] == '#') {
        inf.getline(inBuf,lineSize,'\n');
        continue;
    }
    if ( !strcmp(inBuf,"numberErrorVars") ){
        inf >> dnErrorVars;
    } 
    else if ( !strcmp(inBuf,"globalP")) {
        inf >> dglobalP;
    } 
    else if ( !strcmp(inBuf,"timeStepNumber")) {
        inf >> dtimeStepNumber;
    } 
    else if ( !strcmp(inBuf,"numelX")) {
        inf >> dnumelX;
    } 
    else if ( !strcmp(inBuf,"NSFTAG")) {        
        inf >> dNSFTAG;
    } 
    else if ( !strcmp(inBuf,"ensa_dof")) {
        inf >> densa_dof;
    }
    else if ( !strcmp(inBuf,"attributeFileName")) {
        inf >> dfname;
    }
    else if ( !strcmp(inBuf,"meshFileName")) {
        inf >> dmname;
    }
    else if ( !strcmp(inBuf,"modelFileName")) {     
        inf >> dgname;
    }
    else if ( !strcmp(inBuf,"Idirection")) {
        inf >> dIdirection;
    }
    else if ( !strcmp(inBuf,"BYPASS")) {   
        inf >> dBYPASS;
    }
    else if ( !strcmp(inBuf,"oldPhastaStyle")) {    
        inf >> dfType;
    }
    else if ( !strcmp(inBuf,"zScale")) {
        inf >> dzScale;
    }
    else if ( !strcmp(inBuf,"adaptFlag")) {
        inf >> dadaptFlag;
    }
    else if ( !strcmp(inBuf,"errorName")) {     
        inf >> derrorName;
    }
    else if ( !strcmp(inBuf,"SONFATH")) {
        inf >> dSONFATH;
    }
    else if ( !strcmp(inBuf,"lStart")) {  
        inf >> dlStart;
    }
    else if ( !strcmp(inBuf,"rRead")) { 
        inf >> drRead; 
    }
    else if ( !strcmp(inBuf,"rStart")) { 
        inf >> drStart;
    }
    else if ( !strcmp(inBuf,"AdaptStrategy")) {
        inf >> dstrategy;
    }
    else if ( !strcmp(inBuf,"AdaptFactor")) {
        inf >> dfactor;
    }
    else if ( !strcmp(inBuf,"AdaptOption")) {
        inf >> dadaptOption;
    }
    else if ( !strcmp(inBuf,"hmax")) {
        inf >> dhmax;
    }
    else if ( !strcmp(inBuf,"hmin")) {        
        inf >> dhmin;
    }
    else if ( !strcmp(inBuf,"multipleRestarts")) {
        inf >> dmultipleRestarts;
    }
    else if ( !strcmp(inBuf,"Periodic")) {
        inf >> dper;
    }
    else if ( !strcmp(inBuf,"prCD")) {        
        inf >> dprCD;
    }
    else if ( !strcmp(inBuf,"timing")) {    
        inf >> dtiming;
    } 
    else if ( !strcmp(inBuf,"wGraph")) {
        inf >> dwGraph;
    }
    else if ( !strcmp(inBuf,"phastaVersion")) {
        inf >> dphastaVersion;
    }
    else if ( !strcmp(inBuf,"old_format")) {
        inf >> dold_format;
    } 
    else if ( !strcmp(inBuf,"FortFormFlag")) {   
        inf >> dFortFormFlag;   
    }
    else if ( !strcmp(inBuf,"outputFormat")) {
        inf >> doutputFormat;
    } 
    else if ( !strcmp(inBuf,"CUBES")) {
        inf >> dCUBES;
    }
    else if ( !strcmp(inBuf,"internalBCNodes")) {       
        inf >> dintBC;  
    }
    else if ( !strcmp(inBuf,"version")) {
        inf >> dversion;
    }
    else if ( !strcmp(inBuf,"WRITEASC")) {
        inf >> dWRITEASC;
    }
    else if ( !strcmp(inBuf,"phastaIO")) {  
        inf >> dphastaIO;
    }
    else if ( !strcmp(inBuf,"Epsilon_value")) {
        inf >> depsilon;
cerr << "depsilon = "<<depsilon <<endl;
    }
    else if ( !strcmp(inBuf,"Maximum_size")) {
        inf >> dmax_size;
cerr << "dmax_size = "<<dmax_size <<endl;
    }
    else if ( !strcmp(inBuf,"Minimum_size")) {
        inf >> dmin_size;
cerr << "dmin_size = "<<dmin_size <<endl;
    }
    else if ( !strcmp(inBuf,"Size_flag")) {
        inf >> dsize_flag;
cerr << "dsize_flag = "<<dsize_flag <<endl;
    }
    else if ( !strcmp(inBuf,"refThreshold")) {
        inf >> d_threshold;
    } 
    else if ( !strcmp(inBuf,"parted")) {
        inf >> dparted;
    }
    else if ( !strcmp(inBuf,"preLBforAdaptivity")) {
        inf >> dpreLBforAdaptivity;
    }
    else if ( !strcmp(inBuf,"masterProcWgt")) {
        inf >> dmasterProcWgt;
    }
    else if ( !strcmp(inBuf,"refWeights")) {
      if (dnErrorVars == 0) {
	printf("Error. Check numberErrorVars must be set before weights in adapt.inp\n");
	exit(0);
      }
        d_wght = new double[dnErrorVars];
        for (int i=0; i < dnErrorVars; i++) {
            inf >> d_wght[i];
        }
    }
    else {
      if(PMU_rank()==0)
	printf("\n Warning [in Input()] : <%s> phrase in adapt.inp unrecognized\n\n",inBuf);
    }

    // read rest of the line
    inf.getline(inBuf,lineSize,'\n');
  }//while( 
  
  delete [] inBuf;
}

