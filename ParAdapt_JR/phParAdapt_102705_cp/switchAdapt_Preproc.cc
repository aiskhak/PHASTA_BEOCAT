///////////////////////////////////////////////////////////////////////////
// switchAdapt_Preproc.cc
// switch adaptor or preprocessor 
//
// in case the adaptor is called ahead mesh construction
// is to be handled differently
//
// also a datastructure has to be provided for the solution transfer
///////////////////////////////////////////////////////////////////////////

#include <stdio.h>

#include "phParAdapt.h"
#include "func.h"
#include "ccfunc.h"

#include "MeshSimInternal.h"
#include "SimPMesh.h"
#include "SimMeshTools.h"
#include "MeshSimAdapt.h"

#if ( defined MODELER_PARASOLID )
  #include "SimParasolidKrnl.h"
#elif ( defined MODELER_DISCRETE )
  #include "SimDiscrete.h"
#endif

extern "C" int procSize();

extern int adaptFlag;
extern int timeStepNumber;
extern int rStart;
extern int multipleRestarts;
extern int strategy;
extern int ensa_dof;// number of field(=solution) variables
extern int nErrorVars;
extern double* wght;
extern double factor;
extern double hmax;
extern double hmin;
extern int adaptOption;

extern int preLBforAdaptivity;
extern double masterProcWgt;

int lstep;
extern time_t wtimePoints[32];

int
switchAdapt_Preproc(int argc, char *argv[]){

    // process command line arguments
    // procArgs overwrites geom.sms by the mesh directory name specified
    // fname:attribute file name ,mname mesh directory name 
    procArgs(argc, argv);

    assignGlobalVars();

    lstep = timeStepNumber;

    // model and mesh declaration
    pGModel model; 
    pParMesh pmesh;
    
    pACase acase;
    char    fname[100];
    char    mname[100];
  
  // default attribute/mesh filename
    strcpy(fname, "geom.spj");

    if((adaptFlag && (strategy!=7 || strategy!=8))|| multipleRestarts ){
        strcpy(fname, "geomNOIC.spj");
    }

    strcpy(mname, "geom_p.sms");


    pAManager attmngr = AMAN_retrieve(fname);
    if (attmngr == 0){
      if (PMU_rank() == 0){
	fprintf(stderr, "could not open attribute file %s\n", fname);
	exit(-1);
      }           
    }
    else{
      if (PMU_rank() == 0){
	printf("\n AttMan_retrieve(fname) success\n");
      }
    }
    acase = AMAN_findCase(attmngr, "geom");
    if (acase == NULL){
      char* casename="geom";
      if ( ( acase = AMAN_findCase(attmngr, casename)) == NULL) {
	printf("[%d] Error: could not find attribute case %s\n",PMU_rank(),casename);
	exit (-1);
      }
      else{
	if(PMU_rank() == 0)
	  printf("\n AMAN_findCase success\n");
      }
    }
    // associate the attribute case with the model


    // associate the case with the model (case: e.g. set of BCs)
    AttCase_associate(acase);
    if(PMU_rank()==0)
      printf("\n acase->associate(model) SIM success\n");
    
    model = (pGModel)AttCase_model( acase );
    if(PMU_rank()==0)
      printf("\n AttCase_model( acase ) success\n"); 

#ifdef PARALLEL

    if(PMU_rank()==0) {
      printf("\n");
      printf("Meshing   Library build ID : %s\n",SimMeshing_buildID());
      printf("PMesh     Library build ID : %s\n",SimPMesh_buildID());
      printf("MeshTools Library build ID : %s\n",SimMeshTools_buildID());
    }
    
    wtimePoints[0] = time(0);

    pmesh = PM_new(0,model);
    PM_read(pmesh, mname);
    
    wtimePoints[1] = time(0);

    if(PM_verify  (  pmesh ,0 ) == 0){
      if (PMU_rank() == 0){
	printf("\nerror in adapt.cc: invalid parallel mesh read in\n");
      }
      SimPMesh_stop();
      exit(1);
    }
    
    // cout<<"\n["<<PMU_rank()<<"] pMesh read success...\n";

    // make sure there is only one part per proc
    if(PM_numParts(pmesh) > 1){
      printf("\n[%d] error in adapt:",PMU_rank());
      printf("only one part per proc at the moment allowed \n");
      SimPMesh_stop();
      exit(1);
    }

#endif
    
    // use the adaptor first
    // created new mesh, new restart files
    // total number of DOFS defined here
    int nshgTot=0;
    if(adaptFlag){
        
        printf("\n[%2d] memory usage before mesh adaptation: %d (KB)\n",PMU_rank(),phParAdaptProcSize());

        // flags needed by the preprocessor
        // to directlty take the ICs from the restart files
        // reads files: ./restart.%d.%d
        rStart=1;
        nshgTot = adapt (pmesh,
                         model,
                         timeStepNumber,
                         strategy,
                         factor,
                         // number of field(=solution) variables
                         ensa_dof,
                         // number of variables for error indicators (EI)
                         // (e.g., 5 for ybar & 10 for residual-based)
                         nErrorVars,
                         hmax,
                         hmin,
                         adaptOption);

        printf("\n[%2d] memory usage after mesh adaptation: %d (KB)\n",PMU_rank(),phParAdaptProcSize());
        
    }

    if(nErrorVars)
      delete [] wght; 

    // continue with the usual preprocessing
    mdb2phasta (fname,mname,model, pmesh,nshgTot );

    printf("\n[%2d] memory usage after running preprocessor: %d (KB)\n",PMU_rank(),phParAdaptProcSize());

    
    AttCase_unassociate(acase);

    AMAN_delete(attmngr);

    /* Delete mesh */
    PM_delete(pmesh);

    /* Delete the model. */
    GM_delete(model);
    
    return 1;
    
}
