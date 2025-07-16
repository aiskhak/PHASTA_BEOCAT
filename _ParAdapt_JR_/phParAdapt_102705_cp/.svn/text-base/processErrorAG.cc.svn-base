#include <iostream.h>
#include <fstream>
#include <mpi.h>
#include "phParAdapt.h"

//  using std::cout;
//  using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

extern double* wght;


// processing error indicator values
// provided by phasta according to your needs
double
processErrorAG(double* nodalErrorSet, int nvar)
{
    double scalarVal = 0;

    // using linear weights provided by input file adapt.inp
    for(int i=0; i<nvar;i++){
        
        scalarVal +=  nodalErrorSet[i]* wght[i];
    }
    return scalarVal;
}
        
        



#ifdef __cplusplus
}
#endif
	
