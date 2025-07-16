#include "phParAdapt.h"
#include <iostream.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId nodalGradientID;

// reconstruct the element hessian : 6-component (symmetric)
void
elementHessian(pRegion region, double* elemHessian)
{
  pVertex v;
  double* nodalGradient;
  
  // build the linear system
  double matrix[16];
  buildSystem(region,matrix);
  
  double  fieldVectorX[4];
  double fieldVectorY[4];
  double fieldVectorZ[4];
  int four = 4;
  int indx[4];
  double fnumber;
  
  pPList regionVerts;
  regionVerts = R_vertices(region, 0);
  
  // get the field vals
  for (int i=0; i<PList_size(regionVerts); i++) {
    v = (pVertex)PList_item(regionVerts,i);
    if(!EN_getDataPtr((pEntity)v, nodalGradientID,
		      (void**)&nodalGradient)){
      cout<<"\nerror in elementHessian: no data attached to vertex\n";
      exit(0);
    }
    fieldVectorX[i]= nodalGradient[0];
    fieldVectorY[i]= nodalGradient[1];
    fieldVectorZ[i]= nodalGradient[2];

#if  ( defined  DEBUG )
//       printf("[%i]:\n",PMU_rank());
//       printf("\nELEM HESS: NODAL GRAD :\n %f\n %f\n %f\n",nodalGradient[0],nodalGradient[1],nodalGradient[2]);
//       printf("for vertex:");
//       double c[3];
//       V_coord(v,c);
//       printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif 


  }
  PList_delete(regionVerts);    

  ludcmp_(matrix , &four, &four, indx, &fnumber);
  
  // fieldVector is going to be overridden and will
  // contain the solution
  lubksb_(matrix , &four, &four, indx, fieldVectorX );
  lubksb_(matrix , &four, &four, indx, fieldVectorY );
  lubksb_(matrix , &four, &four, indx, fieldVectorZ );

  // each poly is 
  // a0 + a1*x + a2*y + a3*z  
  elemHessian[0] = fieldVectorX[1];// xx
  elemHessian[1] = fieldVectorX[2];// xy
  elemHessian[2] = fieldVectorX[3];// xz
  elemHessian[3] = fieldVectorY[2];// yy
  elemHessian[4] = fieldVectorY[3];// yz
  elemHessian[5] = fieldVectorZ[3];// zz
  
#if  ( defined  DEBUG )
//    double c[3];
//    for(int i=0; i<6; i++) {
     
//    }
//    printf("for vertex:");
//    V_coord((pVertex)en,c);
//    printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif
 
}
   
#ifdef __cplusplus
}
#endif
