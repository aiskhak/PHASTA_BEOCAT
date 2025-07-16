/* collect the periodic boundary condition array by finding
the degrees of freedom associated with periodic mesh entities */

#include <stdio.h>
#ifndef SIM
#include "MSops.h"
#else
#include "MeshSim.h"
#endif
#include "func.h"

/**********************************************************************/
/* compute the periodic boundary condition array. Periodic boundary   */
/* conditions which reside on other processors are treated as a       */
/* communication, and thus are not included in iper.                  */
/**********************************************************************/
// called in writeEnsaFiles.cc
void attachPeriodicBC(pMesh mesh, globalInfo *info, int *iper)
{
  pEntity ent, em;
  pMatch ment;
  VIter viter;  EIter eiter;  FIter fiter;
  int i, nd;

  viter = M_vertexIter(mesh);

  printf("\nentering attachPeriodicBC\n"); 


  // "IPER" is set setupGlobalTasks.c
  // based on 
  while (ent = (pEntity) VIter_next(viter))  {
    if ((em = EN_dataP(ent, "IPER")) || (isPeriodic(ent, &ment)

            && isOnThisProc(ment) && (em = Match_ent(ment))))                              
        iper[EN_dataI(ent, "MYCT")] = EN_dataI(em, "MYCT");//MYCT from localinfo 
  }
  VIter_delete(viter);

  if (info->edgeson)  {
    eiter = M_edgeIter(mesh);
    while (ent = (pEntity) EIter_next(eiter))  {
      if (nd = EN_dataI(ent, "NDOF"))  {
        if ((em = EN_dataP(ent, "IPER")) || (isPeriodic(ent, &ment)

                   && isOnThisProc(ment) && (em = Match_ent(ment))))  {
          for (i=0; i < nd; i++)
            iper[EN_dataI(ent, "MYCT") + i] = EN_dataI(em, "MYCT") + i;
        }
      }
    }
    EIter_delete(eiter);
  }

  if (info->faceson)  {
    fiter = M_faceIter(mesh);
    while (ent = (pEntity) FIter_next(fiter))  {
      if (nd = EN_dataI(ent, "NDOF"))  {
        if ((em = EN_dataP(ent, "IPER")) || (isPeriodic(ent, &ment)

             && isOnThisProc(ment) && (em = Match_ent(ment))))  {                                 
          for (i=0; i < nd; i++)
            iper[EN_dataI(ent, "MYCT") + i] = EN_dataI(em, "MYCT") + i;
        }
      }
    }
    FIter_delete(fiter);
  }
}
