#ifndef H____func
#define H____func

#include <time.h>
#include <stdio.h>
#include <stdlib.h>


#define MAXNT 100
//#define MAX(a, b) ((a) < (b) ? (b) : (a) )
#define MIN(a, b) ((a) > (b) ? (b) : (a) )

#include "parallel.h"
#include "SimPMesh.h"



#ifdef __cplusplus
extern "C" {
#endif
void UserWeight(pMesh mesh);
void generate_keyphrase(char* target, char* prefix, blockKey* tpblock);
void generate_keyphraseNew(char* target, char* prefix, blockKey* tpblock);
pPList return_vertices(pRegion region);

int topology(pRegion);
void sureToContinue();// just prompt for a y/n sattement
void procArgs(int argc, char *argv[]);


#ifdef PARALLEL
void readBC(pGModel, pParMesh);
void writeNCorp(pParMesh pmesh, globalInfo *info, int **ifath, int **nsons);
void initLocalInfo(pParMesh pmesh, globalInfo* info);
void setupGlobalTasks(pParMesh pmesh, globalInfo *info);
void printInfo(pParMesh pmesh, globalInfo *info);
void genILWork(globalInfo *info, int *ilwork);
#else
void readBC(pGModel model, pMesh mesh);
void writeNCorp(pMesh mesh, globalInfo *info, int **ifath, int **nsons);
void initLocalInfo(pMesh mesh, globalInfo* info);
void printInfo(pMesh mesh, globalInfo *info);

#endif

/* Calcs info->nlwork */
void nlworkCalc(globalInfo *info);

/* returns if ent is periodic, if it is then master is returned thru ment */
int isPeriodic(pEntity ent, pMatch *ment);
/* is master from above on current proc? */
int isOnThisProc(pMatch ment);

void getX(pMesh mesh, double **x);
void getConnectivity(pMesh mesh, pPList bdry, int ***ien, int ***ienb,
                     globalInfo *info);

void R_entities(pRegion region, pVertex *vrts, pEdge *edgs, pFace *fcs,
                int nen);
void R_entitiesBdry(pRegion region, pFace face, pVertex *vrts, pEdge *edgs,
                    pFace *fcs, int nen);
void attachPeriodicBC(pMesh mesh, globalInfo *info, int *iper);

void initGlobalInfo(globalInfo *info);
void freeGlobalInfo(globalInfo *info);

void addSendSegment(globalInfo *info, pEntity ent, pEntity ment, int pid);
void addRecvSegment(globalInfo *info, pEntity ent, int pid);


pPList nsR_edges(pRegion region);


#ifdef __cplusplus
}

#endif

#endif
