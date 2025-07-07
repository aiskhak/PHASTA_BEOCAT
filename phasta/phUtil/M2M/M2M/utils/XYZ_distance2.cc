#include "EN_isPointInside.h"

#ifdef __cplusplus
extern "C" {
#endif

double XYZ_distance2(double P1[3], double P2[3])
{
 double vec[3];

 diffVt(P1,P2,vec);

 return (dotProd(vec,vec));
}

#ifdef __cplusplus
}
#endif
