#ifdef PARALLEL

/************************************************************************/
/* This function sets up the interprocessor communication tasks.      */
/**********************************************************************/
#include <stdio.h>
#include "func.h"
#include "parallel.h"
#include "SimPMesh.h"

// called in mdb2phasta
void setupGlobalTasks(pParMesh pmesh, globalInfo *info)
{
  pEntity ent;
  pEntity em;
  pMatch ment;
  pEntOrig eo;
  pEntCopies ec;

  PBEntProcIter ei;

  void *tmps, *tmpm;
  int i, j, etype, pid, mpid, mult = sizeof(P_int)/sizeof(int);
  int *ns = (int*) malloc(sizeof(int)*PMU_size());
  int *nr = (int*) malloc(sizeof(int)*PMU_size());
  void ***se = (void***) malloc(sizeof(void**)*PMU_size());
  void ***re = (void***) malloc(sizeof(void**)*PMU_size());

  /* first take care of periodic pairs in info that we made in localInfo */
  for (tmps=0, tmpm=0; (ent = (pEntity) PList_next(info->perSlv, &tmps))
                    && (ment = (pMatch) PList_next(info->perMst, &tmpm));)  {
    /* the master is offproc, we made sure of it in localInfo */

      // PMU_proc(), gives rank, given global partition id (gid)
      // Match_gid(), Returns global id of partition on which entity in ment lies.
      pid = PMU_proc(Match_gid(ment));

      addSendSegment(info, ent, Match_ent(ment), pid);

  }

  for (i=0; i < PMU_size(); i++)  {
    ns[i] = PM_numBdryEntOwner(pmesh,0,0,i) + PM_numBdryEntOwner(pmesh,0,1,i)
            + PM_numBdryEntOwner(pmesh,0,2,i);
    se[i] = (void**) malloc(sizeof(void*)*3*ns[i]);
    ns[i] = nr[i] = 0;
  }

  for (etype = 0; etype < 3; etype++)  {
      if ((etype==1 && !info->edgeson) || (etype==2 && !info->faceson)){

          continue;
      }
      // neiter (etype==1 && !info->edgeson) || (etype==2 && !info->faceson)
      ei = PM_bdryProcIter(pmesh, etype);

      while (ent = BdryProcIter_next(ei))  {

          if (isPeriodic(ent, &ment))  {

              mpid = PMU_proc(Match_gid(ment));  
              ec = EN_copies(ent);
              for (i=0; i < EntCopies_size(ec); i++)  {
                  
                  // if the owner proc of the current copy is NOT PMU_rank(),
                  // i.e. it is off-proc
                  // fill up the send array with the copie's (off-proc) pointer
                  if ((pid = PMU_proc(EntCopies_gid(ec,i))) != PMU_rank())  {
                      se[pid][3*ns[pid]] = EntCopies_ent(ec, i);

                      
                      se[pid][3*ns[pid]+1] = Match_ent(ment);
                      se[pid][3*ns[pid]+2] = (void*) mpid;
                      ns[pid]++;
                  }
              }
          }
      }
      BdryProcIter_delete(ei);
      
  }// for (etype ...
  // Communicates integer arrays across processes, only in parallel
  // collective call
  PMU_commuInt(ns, nr);
  for (i=0; i < PMU_size(); i++){
      
      re[i] = (void**) malloc(sizeof(void*)*3*nr[i]);
  }

  // Communicates arrays across processes, only in parallel
  PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, 3*mult);

  for (i=0; i < PMU_size(); i++)  {
    for (j=0; j < nr[i]; j++)  {
      ent = (pEntity) re[i][3*j];
      EN_attachDataI(ent, "REI+", i+1);  /* +1 is to avoid the 0 ambiguity */
      EN_attachDataI(ent, "REJ+", j);
    }
  }

  for (etype = 0; etype < 3; etype++)  {
      if ((etype==1 && !info->edgeson) || (etype==2 && !info->faceson)){
          
          continue;
      }
      // neiter (etype==1 && !info->edgeson) || (etype==2 && !info->faceson)
      ei = PM_bdryProcIter(pmesh, etype);

      while (ent = BdryProcIter_next(ei))  {
          
          if ((pid = EN_ownerProc(ent)) != PMU_rank())  {  /* found a slave */

              if (i = EN_dataI(ent, "REI+"))  {

                  j = EN_dataI(ent, "REJ+");
                  em = (pEntity) re[i-1][3*j+1];
                  mpid = (P_int) re[i-1][3*j+2];
                  if (mpid != PMU_rank())
                      addSendSegment(info, ent, em, mpid);
                  else
                      EN_attachDataP(ent, "IPER", em);
              }
              else  {
                  eo = EN_original(ent);
                  addSendSegment(info, ent, EntOrig_ent(eo), pid);
                  EntOrig_delete(eo);
              }
          }
      }
      BdryProcIter_delete(ei); 
  }//for

  /* setup send stuff */
  for (i=0; i < PMU_size(); i++)  {
    free(se[i]);
    free(re[i]);

    nr[i] = 0;
    ns[i] = info->stask[i] ? info->stask[i]->numSeg : 0;
    se[i] = (void**) malloc(sizeof(void*)*ns[i]);
    for (j=0; j < ns[i]; j++)
      se[i][j] = (pEntity) PList_item(info->stask[i]->ments, j);
  }

  // Communicates integer arrays across processes, only in parallel
  // collective call
  PMU_commuInt(ns, nr);

  for (i=0; i < PMU_size(); i++)
    re[i] = (void**) malloc(sizeof(void*)*nr[i]);

  PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);

  for (i=0; i < PMU_size(); i++)  {
    for (j=0; j < nr[i]; j++)
      addRecvSegment(info, (pEntity) re[i][j], i);
  }

  for (i=0; i < PMU_size(); i++)  {
    free(se[i]);
    free(re[i]);
  }
  free(se);
  free(re);
  free(ns);
  free(nr);
}

#endif
