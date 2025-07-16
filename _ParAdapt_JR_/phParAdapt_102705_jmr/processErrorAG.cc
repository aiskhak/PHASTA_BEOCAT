#include <iostream.h>
#include <fstream>
#include <mpi.h>
#include "phParAdapt.h"
#include <math.h>
//  using std::cout;
//  using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

extern double* wght;
extern double epsilon;


// processing error indicator values
// provided by phasta according to your needs
double
processErrorAG(double* nodalErrorSet, double* nodalSolutionSet, int nvar, int option)
{
    double scalarVal = 0;
    double phi;

// option =2 - base error on proximity to interface
    if (option == 2) {
      phi = sqrt(nodalSolutionSet[6]*nodalSolutionSet[6]);
//    cout << "phi = "<<phi<<"\n";
      if (phi < epsilon) {
        scalarVal = 1.0 - 0.5*(1.0-phi/epsilon+1.0/3.141593*sin(3.141593*phi/epsilon)); 
      } else {
        scalarVal = 0.001;
      }
    }
    else {
    // using linear weights provided by input file adapt.inp
    for(int i=0; i<nvar;i++){
        
        scalarVal +=  nodalErrorSet[i]* wght[i];
      }
    }
//
    return scalarVal;
}
        
        



#ifdef __cplusplus
}
#endif
	
