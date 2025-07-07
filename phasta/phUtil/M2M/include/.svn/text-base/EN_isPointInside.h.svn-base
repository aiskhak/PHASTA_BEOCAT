#ifndef _EN_ISPOINTINSIDE_H_
#define _EN_ISPOINTINSIDE_H_

#include "MeshSim.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef double dArray[3];

  double M_getTolerance();
  // internal functions
  void diffVt(double [],double [],double []);
  double dotProd(double [],double []);
  void crossProd(double [],double [],double []);
  double XYZ_distance2(double [3], double [3]);
  double XYZ_volume(dArray *xyz);
  void F_coord(pFace, dArray*);
  // internal functions

  // interpolate solution
  double* InterpolateSolution( pRegion region,
                               pMeshDataId nodalSolutionID,
                               double xi[3],
                               int ndof );

  int inverseMap( pRegion region,
                  double* qpt,
                  double* pt );
  // interpolate solution

  int EN_isPointInside(pEntity entity, double *xyz);

  int F_isPointInside(pFace face, double *xyz);
  int XYZ_isPointInsideTriangle(double tri_xyz[3][3], double *xyz);

  double XYZ_volume2(dArray *X_Y_Z, int rTopoType);
  double R_volume2(pRegion region);

  int R_isPointInside(pRegion region, double *xyz);
  int XYZ_isPointInside(dArray *X_Y_Z, double *xyz);

#ifdef __cplusplus
}
#endif

#endif
