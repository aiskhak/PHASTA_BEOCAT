#include "phParAdapt.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int BCInflowFaceTag;

extern pMeshDataId vertID;
extern pMeshDataId inflowOwnedVertID;

extern 
void PMU_commuInt(int* ns, int* nr);
extern
void PMU_commuArr(void **s, int *ns, void **r, int *nr, MPI_Datatype type,
                  int mult);

// communicates contributions for verticies on inflow face and
// on partition bdry.

// owner vertex on inflow face and on part. bdry.
// (and have been tagged with vertID)
// take the id of owner vertex in connectivity

// data are gathered via
// vertID
void
commuBCInflowFaceInfo(pParMesh pmesh, pMesh mesh)
{
    int mult = sizeof(P_int)/sizeof(int);
    int *ns = (int*) malloc(sizeof(int)*PMU_size());
    int *nr = (int*) malloc(sizeof(int)*PMU_size());
    int *nsNodes = (int*) malloc(sizeof(int)*PMU_size());
    int *nrNodes = (int*) malloc(sizeof(int)*PMU_size());
    void ***se = (void***) malloc(sizeof(void**)*PMU_size());
    void ***re = (void***) malloc(sizeof(void**)*PMU_size());
    pEntity ent;
    int** dsend = new int*[PMU_size()];
    int** drecv = new int*[PMU_size()];

    for(int i=0; i<PMU_size(); i++) {
        ns[i] = 0;
	nsNodes[i] = 0;
    }

    int inflowOwnedVertTag;
    pVertex vertex;
    VIter vit = M_vertexIter(mesh);
    while(vertex = VIter_next(vit)) {
      if(EN_getDataInt((pEntity)vertex,inflowOwnedVertID,
		       &inflowOwnedVertTag)) {
	// need this for accumulatedID
	// ids of owned vertex by higher rank procs would
	// be added by accumulatedID
	for(int i=0; i<PMU_size(); i++)
	  nsNodes[i]++;
      }
    }
    VIter_reset(vit);

    pEntity en;
    pBdryProcIter myBdryIter = PM_bdryProcIter(pmesh,0);
    // determine the size of *ns[i]
    while(en = BdryProcIter_next(myBdryIter)) {
      if(EN_getDataInt((pEntity)en,inflowOwnedVertID,
		       &inflowOwnedVertTag)) {
	pEntCopies vCopies = EN_copies(en); 

	// how many of them on OTHER parts/=procs ?
	int vCopiesSize = EntCopies_size(vCopies);
        
	// fill the send array with ents/vertices
	for(int i=0; i<vCopiesSize; i++) {

	  // global id of partition/proc on which en's n'th copy lies.
	  int vCopyGid = EntCopies_gid(vCopies,i);

	  if (vCopyGid != PMU_rank()) {
	    // increment number of copies shared between this proc and proc
	    // with rank vCopyGid
	    ns[vCopyGid]++;
	  }
	}
      }
    }
    BdryProcIter_reset(myBdryIter);
    
    for(int i=0; i < PMU_size(); i++) {
      // ns[i] is   number of ent-copies shared by this proc with
      // proc i

      // each se[i] is of length ns[i] and
      // therefore each se[i] will be filled with (off-proc) copies of all ents
      // that are shared with proc i

      // along with that we communicate an int array of length 2*ns[i]: dsend, drecv
      // -- id. of owned vertex
      // -- owner gid
        se[i] = (void**) malloc(sizeof(void*)*ns[i]);
        dsend[i] = new int[2*ns[i]];
        ns[i] = nr[i] = 0;
	nrNodes[i] = 0;
    }

    BdryProcIter_reset(myBdryIter);
    while(en = BdryProcIter_next(myBdryIter)) {
      if(EN_getDataInt((pEntity)en,inflowOwnedVertID,
		       &inflowOwnedVertTag)) {
        // copies of the vertex on other parts/=procs
        // does not include the entity en that vCopies was created with
        pEntCopies vCopies = EN_copies(en); 

        // how many of them on OTHER parts/=procs ?
        int vCopiesSize = EntCopies_size(vCopies );
        
        // fill the send array with ents/vertices
        for(int i=0; i<vCopiesSize; i++) {

	  // global id of partition/proc on which en's n'th copy lies.
	  int vCopyGid = EntCopies_gid(vCopies,i);

	  if (vCopyGid != PMU_rank()) {
	    // fill send array to proc of rank vCopyGid with the
	    // (ns[vCopyGid])-th copy
	    // among all copies of verts shared
	    // between this proc and proc vCopyGid
	    se[vCopyGid][ ns[vCopyGid] ] = EntCopies_ent(vCopies,i);

	    int ID;
	    // get the local data
	    if(!EN_getDataInt(en,vertID,&ID)) {
	      cout<<"\nerror in commuBCInflowFaceInfo: no vertex id. attached to vertex (during send)\n";
	      exit(0);
	    }
		
	    dsend[vCopyGid][ 2*ns[vCopyGid] ] = ID;
	    dsend[vCopyGid][ 2*ns[vCopyGid]+1 ] = PMU_rank();

	    // increment number of copies shared between this proc and proc
	    // with rank vCopyGid
	    ns[vCopyGid]++;
	  }
        }
      }
    }
    BdryProcIter_reset(myBdryIter);
    BdryProcIter_delete(myBdryIter);

    // communicate the NUMBER of copies that are shared with other procs
    PMU_commuInt(ns, nr);
    // communicate the NUMBER of verts owned by each proc
    PMU_commuInt(nsNodes, nrNodes);

    // allocate length of recv array
    // the size (number) of verts copies that are shared with this proc
    // on each proc i
    for(int i=0; i<PMU_size(); i++) {
        re[i] = (void**) malloc(sizeof(void*)*nr[i]);
        drecv[i] = new int[2*nr[i]];
    }
    
    // communicate the actual copies and data
    PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);
    PMU_commuArr((void**)dsend, ns, (void**)drecv, nr, MPI_INT, 2);

    int accumulatedID = 0;
    // add for the procs with lower ranks (excluding this proc)
    for(int i=0; i<PMU_rank(); i++)
      accumulatedID += nrNodes[i];

    // retrieve communicated data and globalize
    for(int i=0; i<PMU_size(); i++) {
      // nr[i] is number of ent-copies shared by this proc with
      // proc i
      for(int j=0; j<nr[i]; j++) {            
	// the vertex itself
	ent = (pEntity) re[i][j];

	int ID;
	if(EN_getDataInt(ent,vertID,&ID)) {
	  cout<<"\nerror in commuBCInflowFaceInfo globalize: already have vertex id. data attached to vertex (during receive)\n";
	  exit(0);
	}
	
	int accuID = 0;
	for(int k=0; k<drecv[i][2*j+1]; k++) {
	  if(k==PMU_rank())
	    accuID += nsNodes[0];
	  else
	    accuID += nrNodes[k];
	}
	
	ID = drecv[i][2*j] + accuID;
	EN_attachDataInt(ent,vertID,ID);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int id;
    VIter_reset(vit);
    while(vertex = VIter_next(vit)) {
      if(EN_getDataInt((pEntity)vertex,inflowOwnedVertID,
		       &inflowOwnedVertTag)) {
	if(!EN_getDataInt(vertex,vertID,&id)) {
	  cout<<"\nerror in commuBCInflowFaceInfo globalize: no vertex id. data attached to vertex\n";
	  exit(0);
	}
	
	id = id + accumulatedID;
	EN_modifyDataInt((pEntity)vertex,vertID,id);
      }
    }
    VIter_delete(vit);

    for(int i=0; i<PMU_size(); i++)  {
        free(se[i]);
        free(re[i]);
        delete [] dsend[i];
        delete [] drecv[i];
    }

    delete [] dsend;
    delete [] drecv;

    free(se);
    free(re);
    free(ns);
    free(nr);
    free(nsNodes);
    free(nrNodes);
}

#ifdef __cplusplus
}
#endif
