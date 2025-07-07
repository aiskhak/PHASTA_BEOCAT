#include "EN_isPointInside.h"

#ifdef __cplusplus
extern "C" {
#endif

double XYZ_volume (dArray *xyz)
{
  double v01[3],v02[3],v03[3];
  double normal[3];


  diffVt(xyz[1],xyz[0],v01);
  diffVt(xyz[2],xyz[0],v02);
  crossProd(v01,v02,normal);
  diffVt(xyz[3],xyz[0],v03);

  return(dotProd(normal,v03)/6);
}

#ifdef __cplusplus
}
#endif

