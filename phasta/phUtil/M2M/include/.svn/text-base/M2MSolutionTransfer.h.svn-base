#ifndef _M2MSOLUTIONTRANSFER_H_
#define _M2MSOLUTIONTRANSFER_H_

#include "MeshSim.h"

#ifdef __cplusplus
extern "C" {
#endif

void M2MSolutionTransfer(pGModel source_model, pGModel dest_model, pMesh source_mesh, pMesh dest_mesh, pMeshDataId solutionID, int numVars, int rfOption, int farOption, int halfToFullOption, int entTagsOption);

// memory is allocated for solution on each vertex
void M2MSolutionTransfer0(pGModel source_model, pGModel dest_model, pMesh source_mesh, pMesh dest_mesh, pMeshDataId solutionID, int numVars, int farOption, int halfToFullOption, int entTagsOption);

// memory is allocated for solution on each vertex
void M2MSolutionTransfer1(pGModel source_model, pGModel dest_model, pMesh source_mesh, pMesh dest_mesh, pMeshDataId solutionID, int numVars, int farOption, int halfToFullOption, int entTagsOption);

#ifdef __cplusplus
}
#endif

#endif
