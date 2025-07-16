#include <stdio.h>
#include <iostream.h>
#include <fstream>
#include <mpi.h>
#include "phParAdapt.h"

//  using std::cout;
//  using std::endl;

#ifdef __cplusplus
extern "C" {
#endif


extern pMeshDataId errorIndicatorID;


// nvar is number of error indicators
// as provided by phasta
void
transformToScalarErrorVal(pMesh mesh, int nvar)
{

    // loop over vertices
    pVertex v;
    VIter vIter=M_vertexIter(mesh);

    while(v = VIter_next(vIter)) {

        double *nodalErrorSet;
        double* scalarValue = new double;
        if(!EN_getDataPtr((pEntity)v, errorIndicatorID,(void**)&nodalErrorSet)){
            
            cout<<"\nerror in transformToScalarErrorVal: no data attached to  vertex\n";
            V_info(v);
            exit(0);
        }
        *scalarValue = processErrorAG(nodalErrorSet,nvar);

	delete [] nodalErrorSet;
        
        EN_deleteData((pEntity)v,errorIndicatorID );


        EN_attachDataPtr( (pEntity)v, errorIndicatorID, (void *)
		       scalarValue);


#ifdef DEBUG
        double *testScalarValue;
        if(!EN_getDataPtr((pEntity)v, errorIndicatorID,(void**)&testScalarValue)){
            
            cout<<"\nerror in transformToScalarErrorVal: no data attached to  vertex\n";
            V_info(v);
            exit(0);
        }
//         else{
//             cout<<"\nin transformToScalarErrorVal: attached a scalar"
//                 <<"errorIndicatorID value : "<<testScalarValue[0]<<"\n";
//         }
#endif//DEBUG        


    }
    VIter_delete(vIter);

}

#ifdef __cplusplus
}
#endif
	
