#include "phParAdapt.h"
#include <stdio.h>
#include <unistd.h>
#include <iostream.h>

#ifdef __cplusplus
extern "C" {
#endif

extern pMeshDataId errorIndicatorID;


// reconstruct the element gradient
void
elementGradient(pRegion region, double* elementGradient)
{
  pVertex v;
  double* nodalData;
  
  // build the linear system
  double matrix[16];
  buildSystem(region,matrix);
  double fieldVector[4];
  
  int four = 4;
  int indx[4];
  double fnumber;
  
  pPList regionVerts;
  regionVerts = R_vertices(region, 0);
  
  // get the field vals
  for (int i=0; i<PList_size(regionVerts); i++) {
    v = (pVertex)PList_item(regionVerts,i);
    if(!EN_getDataPtr((pEntity)v, errorIndicatorID,
		      (void**)&nodalData)){
      cout<<"\nerror in elementGradient: no data attached to vertex\n";
      exit(0);
    }
//  #if  ( defined  DEBUG )
//    printf("\nretrieved EI data: %f \n",nodalData[4]);
//  #endif 


    fieldVector[i]=nodalData[4];// speed component/ temperature
    
#if  ( defined  DEBUG )
//     printf("\n[%i]: ELEM GRAD: NODAL DATA :\n %f\n",PMU_rank(),nodalData[4]);
//     printf("for vertex:");
//     double c[3];
//     V_coord(v,c);
//     printf("coords %f %f %f\n",c[0],c[1],c[2]);
#endif 


  }
  PList_delete(regionVerts);    
  
  ludcmp_(matrix , &four, &four, indx, &fnumber);
  
  // fieldVector is going to be overridden and will
  // contain the solution
  lubksb_(matrix , &four, &four, indx, fieldVector );
  
  for(int i=0; i<3; i++) {
    // the poly's structure is
    // a0 + a1*x + a2*y + a3*z
    elementGradient[i] = fieldVector[i+1];
  }
//  #if  ( defined  DEBUG )
//    printf("\nretrieved elem gradient: %f %f %f\n",elementGradient[0],elementGradient[1],elementGradient[2]);
//  #endif 
}
   
#ifdef __cplusplus
}
#endif
