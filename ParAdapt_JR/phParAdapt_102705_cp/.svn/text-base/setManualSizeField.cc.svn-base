#include "phParAdapt.h"
#include "Eigen.h"
#include <iostream.h>
#include <fstream.h>
#include "attachData.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>


//  using std::cout;
//  using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId phasta_solution;
extern pMeshDataId errorIndicatorID;
extern pMeshDataId nodalGradientID;
extern pMeshDataId nodalHessianID;


void 
setManualSizeField(pMesh mesh,
		   pMSAdapt simAdapter, 
        	   int strategy) {
  double dir[3][3];
  char option[28];
  
  switch(strategy) {
  case -1:
    {
      sprintf(option,"constant");

      dir[0][0]=0.1;//x-direction
      dir[1][1]=0.1;//y-direction
      dir[2][2]=6.0;//z-direction
    
      dir[0][1]=dir[0][2]=0.;
      dir[1][0]=dir[1][2]=0.;
      dir[2][0]=dir[2][1]=0.;

      pVertex vertex;
      VIter vit=M_vertexIter(mesh);
      while(vertex=VIter_next(vit)) {   
	MSA_setAnisoVertexSize(simAdapter, 
			       vertex,
			       dir);
      }
      VIter_delete(vit);
    }
    break;
  case -2:
    {
      sprintf(option,"cylindrical");

      double norm,tol=1.e-8;
    
      double sizeR=0.1;
      double sizeTheta=0.1;
      double sizeZ=1.0;

      pVertex vertex;
      VIter vit=M_vertexIter(mesh);
      while(vertex=VIter_next(vit)) {  
	double xyz[3];
	V_coord(vertex,xyz);
	norm=sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	if( norm>tol ) 
	  {
	    dir[0][0]=sizeR*(xyz[0]/norm);
	    dir[0][1]=sizeR*(xyz[1]/norm);
	    dir[0][2]=0.;
	    dir[1][0]=-sizeTheta*(xyz[1]/norm);
	    dir[1][1]=sizeTheta*(xyz[0]/norm);
	    dir[1][2]=0.;
	    dir[2][0]=0.;
	    dir[2][1]=0.;
	    dir[2][2]=sizeZ;
	  }
	else
	  {
	    dir[0][0]=sizeR;
	    dir[0][1]=0.;
	    dir[0][2]=0.;
	    dir[1][0]=0.;
	    dir[1][1]=sizeTheta;
	    dir[1][2]=0.;
	    dir[2][0]=0.;
	    dir[2][1]=0.;
	    dir[2][2]=sizeZ;
	  }

	MSA_setAnisoVertexSize(simAdapter, 
			       vertex,
			       dir);
      }
      VIter_delete(vit);
    }
    break;
  default :
    cout<<"check strategy [in setManualSizefield(...)]"<<endl;
    exit(-1);
    break; 
  }

  ofstream adaptSimLog("phAdapt.log");
  adaptSimLog<<"Strategy chosen for adaptation is size-field driven"<<endl;
  adaptSimLog<<"Mesh size-field is set manually"<<endl;
  adaptSimLog<<"Size-field option : "<<option<<endl;
  adaptSimLog.close();
}

#ifdef __cplusplus
}
#endif
