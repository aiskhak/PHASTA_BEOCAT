#include "phParAdapt.h"
#include "shapeFunction.h"
#include "MeshSimInternal.h"
#include "SimPMesh.h"
#include "MeshSim.h"
#include "MeshSimAdapt.h"

#include <stdio.h>
#include <stdlib.h>

extern pMeshDataId phasta_solution;


double* 
InterpolateSolution( pRegion region, 
                     double xi[3], 
                     int ndof,
                     pMeshDataId modes ) {

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
    double r=xi[0];
    double s=xi[1];
    double t=xi[2];
    double L[4] = { r, s, t, 1-r-s-t};
    
    for (int i=0; i < ndof; i++) q[i] = 0.0;
  
    //pVertex vrts[4];
    pPList vrts;
    pEdge edgs[6];
    pFace facs[4];


// may mal function depending on entities that come together with "region"
//    R_entitiesAdapt(region, vrts, edgs, facs);
    vrts = R_vertices(region,0);


  
    for (int k=0; k < 4; k++) {        // loop over the vertices of (parent) region
        

        if(! EN_getDataPtr( (pEntity)PList_item(vrts,k), phasta_solution, (void
                                                               **)&vcc )){
           
            V_info((pVertex)PList_item(vrts,k));
            printf("\n[%d]error in InterpolateSolution: wanted to retrieve a solution on a vertex "
                "that is empty (phasta_solution)\n",PMU_rank());
            
            exit(-1);
            
        }

        /*evaluate the shape function at the given location xi*/
        sfval = Entity_shapeFunction(( pEntity )(pVertex)PList_item(vrts,k),
                                     ( pEntity )region, 1, 0, L);
        for (int i=0; i < ndof; i++){
            /*sigma(N_a(xi)*y_a)*/
            q[i] += vcc[i]*sfval;
        }
    }
    PList_delete(vrts);
  
    // get the edge and face coefficients
    // not implemented at the moment
    // may be in conflict with callback function paradigm:
    // a prefefined set of adjacent entities would have to be passed
    // along with a given entity
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
