#ifdef PARALLEL

/**********************************************************************/
/* this function generates the local work array for parallel */
/* communication for the given processor */
/**********************************************************************/
// called in writeEnsaFiles.cc: writeEnsaFiles
#include <stdio.h>
#include "func.h"
#include "parallel.h"
#include "assert.h"


void ilworkTask(Task *t, int *ilwork, int *count, FILE *fp)
{
  int i;
  pEntity ent;
  if (t)  {
    ilwork[0]++;
    ilwork[*count+1] = t->tag;
    ilwork[*count+2] = t->type;
    ilwork[*count+3] = t->other_pid + 1;
    ilwork[*count+4] = t->numSeg;
   

    fprintf(fp, "%d %d %d %d\n", t->tag, t->type, t->other_pid, t->numSeg);
    for (i=0; i < t->numSeg; i++)  {
      ent = (pEntity) PList_item(t->ents, i);
      ilwork[*count + 5 + 2*i] = EN_dataI(ent, "MYCT") + 1;  /* Fortran idx */
      ilwork[*count + 6 + 2*i] = EN_dataI(ent, "NDOF");
    }
    *count += 4 + 2*t->numSeg;
  }
}


void genILWork(globalInfo *info, int *ilwork)
{
  char fname[20];
  FILE *ilwf;
  int i, count;

  sprintf(fname, "graph.out.%d", PMU_rank()+1);
  ilwf = fopen(fname, "w");
  

  fprintf(ilwf, "tag, type, other, numSeg\n");

  ilwork[0] = 0;
  count = 0;

  /* Exchange ordering of ilwork from all-recvs-follow-all-sends to vice-versa
   * by swapping rtask and stask below */
  for (i=0; i < PMU_size(); i++)
    ilworkTask(info->rtask[i], ilwork, &count, ilwf);
  for (i=0; i < PMU_size(); i++)
    ilworkTask(info->stask[i], ilwork, &count, ilwf);
  assert(count+1 == info->nlwork);  /* +1 for ilwork[0] which has numtasks */
  fclose(ilwf);
}

#endif
