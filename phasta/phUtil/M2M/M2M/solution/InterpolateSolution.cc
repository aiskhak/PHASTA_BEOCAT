#include "MeshSim.h"

#include <iostream>
#include <stdlib.h>

using std::cout;
using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

double* 
InterpolateSolution( pRegion region, 
                     pMeshDataId nodalSolutionID,
                     double xi[3], 
                     int ndof ) 
{

/*    This routine takes the region and the location of the vertex in */
/*    parametric coordinates, polynomial order, and the attached solution*/
/*    and returns the solution interpolated at the specified location. Higher*/
/*    order mode coefficients are also evaluated    */
    int exitCode=0;
    double* vcc;
    double* ecc;
    double* fcc;
    double sfval;
    double* q = new double[ndof];
    double u=xi[0];
    double v=xi[1];
    double w=xi[2];
    // max. index is 8 (for hex)
    double L[8];

    int rTopoType = R_topoType(region);
    switch(rTopoType) {
    case Rtet : {
      L[0] =  1.-u-v-w;
      L[1] =  u;
      L[2] =  v;
      L[3] =  w;
    }
    break;
    case Rpyramid : {
      double ui, vi;
      ui = 2.*u/(1.-w); vi = 2.*v/(1.-w);

      L[0] = 0.125*(1.-ui)*(1.-vi)*(1.-w);
      L[1] = 0.125*(1.+ui)*(1.-vi)*(1.-w);
      L[2] = 0.125*(1.+ui)*(1.+vi)*(1.-w);
      L[3] = 0.125*(1.-ui)*(1.+vi)*(1.-w);
      L[4] = 0.5*(1.+w);
    }
    break;
    case Rwedge : {
      double k = 1.-u-v;
      L[0] =  (k) * (0.5*(1.-w));
      L[1] =  (u) * (0.5*(1.-w));
      L[2] =  (v) * (0.5*(1.-w));
      L[3] =  (k) * (0.5*(1.+w));
      L[4] =  (u) * (0.5*(1.+w));
      L[5] =  (v) * (0.5*(1.+w));
    }
    break;
    default:
      cout<<"\nError in InterpolateSolution()..."<<endl;
      cout<<"element topology ["<<rTopoType<<"] NOT supported\n"<<endl;
      exit(0);
    }

    for (int i=0; i < ndof; i++) q[i] = 0.0;
  
//     pVertex vrts[4];
//     pEdge edgs[6];
//     pFace facs[4];

//     R_entities(region, vrts, edgs, facs);

    pPList verts = R_vertices(region,0);
    int nshl = PList_size(verts);

    pVertex *vrts = new pVertex[nshl];
    for(int k=0; k<nshl; k++)
      vrts[k] = (pVertex)PList_item(verts,k);
    PList_delete(verts);
  
    for (int k=0; k<nshl; k++) { // loop over the vertices of (parent) region
        
      if(! EN_getDataPtr( ( pEntity )vrts[k], nodalSolutionID, (void
								**)&vcc )){
	
	cout<<"\nerror in InterpolateSolution: wanted to retrieve a solution on a vertex "
	    <<"that is empty\n";
	
	exit(-1);
      }

      /*evaluate the shape function at the given location xi*/
//       sfval = Entity_shapeFunction(( pEntity )vrts[k],
// 				   ( pEntity )region, 1, 0, L);
      sfval = L[k];

      for (int i=0; i < ndof; i++){
	/*sigma(N_a(xi)*y_a)*/
	q[i] += vcc[i]*sfval;
      }
    }

    delete [] vrts;

    // get the edge and face coefficients
    // not implemented at the moment
//      for (int k=0; k < 6; k++) {    // loop over edges
//          int nem =0;
//          solution(( pEntity )edgs[k], phSol, &ecc ); 
//          numberofmodes((pEntity)edgs[k], modes, &nem );
//          if ( ecc ) {
//              for( int e=0; e < nem ; e++ ) {
//                  sfval = Entity_shapeFunction(( pEntity )edgs[k],
//                                               ( pEntity )region, e+2, 0, L);
//                  for (int j=0; j < ndof; j++) 
//                      q[j] += ecc[e*ndof+j]*sfval;
//              }
//          }
//      }
//      for (int k=0; k < 4; k++) {    // loop over faces
//          int nem=0;
//          solution(( pEntity )facs[k], phSol, &ecc ); 
//          numberofmodes((pEntity)facs[k], modes, &nem );
//          if ( ecc ) {
//              for( int e=0; e < nem ; e++ ) {
//                  sfval = Entity_shapeFunction(( pEntity )facs[k],
//                                               ( pEntity )region, e+3, 0, L);
//                  for (int j=0; j < ndof; j++) 
//                      q[j] += ecc[e*ndof+j]*sfval;
//              }
//          }
//      }
    return q;
}

#ifdef __cplusplus
}
#endif
