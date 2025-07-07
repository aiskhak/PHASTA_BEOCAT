#ifndef H_MeshSimInternal1
#define H_MeshSimInternal1

#include "MeshSim.h"
#include "SimPMesh.h"
#include "MeshSimAdapt.h"


typedef double dArray[3];

#ifdef __cplusplus
extern "C" {
#endif

void diffVt(double [],double [],double []);
int normVt(double [], double []);
double dotProd(double [],double []);
void crossProd(double [],double [],double []);

void M_setTolerance(pMesh);

int GEN_dataI(pGEntity , const char *tag);
void F_coord(pFace, dArray*);
int inList(pPList, void*);
void * GEN_dataP(pGEntity , const char *tag);
pEdge E_exists(pVertex v1, pVertex v2);
pEdge R_gtOppEdg(pRegion region, pEdge edge);
pFace   R_vtOpFc(pRegion ,pVertex);
pVertex R_fcOpVt(pRegion rgn, pFace face);
double R_volume(pRegion);
int F_conToFace(pFace,pFace);
void EN_attachDataI(pEntity , const char *tag, int data);
int EN_dataI(pEntity , const char *tag);
void * EN_dataP(pEntity , const char *tag);

void R_info(pRegion);
void F_info(pFace);
void E_info(pEdge);

void F_setWhatIn(pFace  , pGEntity what);
void E_setWhatIn(pEdge  , pGEntity what);
double E_length(pEdge edge);

pVertex F_oppositeVertex(pFace face, pEdge edge);

pVertex M_nextVertex(pMesh, void **restart);

void F_normalVector(pFace face, int dir, double* normal);
void F_chDir(pFace);

void GEN_attachDataP(pGEntity , const char *tag, void *data);
void GEN_attachDataI(pGEntity , const char *tag, int data);
int GEN_modifyDataP(pGEntity, const char *tag, void * data);
void EN_attachDataP(pEntity , const char *tag, void *data);


// left overs
// JM 8/9/04
pPList R_verticesLeft(pRegion);
void V_info(pVertex);
pEdge M_nextEdge(pMesh, void **restart);
pFace M_nextFaceCancel(pMesh, void **restart);
pFace M_nextFace(pMesh, void **restart);
pRegion M_nextRegionCancel(pMesh, void **restart);
pRegion M_nextRegion(pMesh, void **restart);


//from PMops.h
#ifdef B64
  typedef long P_int;
#else
  typedef int P_int;
#endif





#ifdef __cplusplus
}
#endif

typedef struct BdryProcIter *PBEntProcIter;
typedef struct BdryProcIter *PBFProcIter;
typedef struct BdryProcIter *PBEProcIter;
typedef struct BdryProcIter *PBVProcIter;

typedef struct BdryPartIter *PBEntPartIter;
typedef struct BdryPartIter *PBFPartIter;
typedef struct BdryPartIter *PBEPartIter;
typedef struct BdryPartIter *PBVPartIter;

typedef struct BdryNPartsIter *PBEntNPartsIter;
typedef struct BdryNPartsIter *PBFNPartsIter;
typedef struct BdryNPartsIter *PBENPartsIter;
typedef struct BdryNPartsIter *PBVNPartsIter;


#endif /* not H_MeshSimInternal */
