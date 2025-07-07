#include "func.h"
#include "parallel.h"

#ifdef PARALLEL
void printInfo(pParMesh pmesh, globalInfo *info)
#else
void printInfo(pMesh mesh, globalInfo *info)
#endif
{
#ifdef DEBUG_PER
  pEntity ent;
  int count=0;
  int nem, nfm, nrm, nsh, nsb;
  void *temp = 0;
#endif
  char fn[20];
  FILE *fp;

#ifdef PARALLEL
  pMesh mesh = PM_mesh(pmesh, 0);
  sprintf(fn, "info.out.%d", PMU_rank()+1);
#else
  strcpy(fn, "info.out");
#endif

  fp = fopen(fn, "w");

  fprintf(fp, "\n==========> ENSA-5C preprocessing information <==========\n\n"
             "\t%d (times n_flow) total degrees of freedom (dofs) \n", info->nshgTot);
  fprintf(fp, "Global mesh entity information \n"
             "------------------------------ \n"
             "\t%6d Vertices\n\t%6d Edges\n\t%6d Faces\n\t%6d Regions\n\n",
   M_numVertices(mesh), M_numEdges(mesh), M_numFaces(mesh), M_numRegions(mesh));

  fprintf(fp, "\n\nIndividual partition information \n"
                 "--------------------------------\n"
             "Total number of processors: %d\n\n", PMU_size());
  fprintf(fp, "  Processor %d \n", PMU_rank());
  fprintf(fp, "  ------------ \n"
             "\t%6d elements\n\t%6d vertices\n\t%6d boundary elements\n\t%6d "
             "prescribed essential BCs\n\t%6d is the dimension of the local "
             "work array\n\n",
             info->numel, info->numnp, info->numelb, info->numpbc,
             info->nlwork);

#ifdef DEBUG_PER
  fprintf(fp, "\nLocal mesh information \n");
  fprintf(fp, "----------------------- \n\n");

  fprintf(fp, "Vertices: (id, master)\n");
  while (ent=M_nextVertex(mesh, &temp))
    fprintf(fp, "V%-2d ---------> %d\n", count++, EN_ownerProc(ent));

  temp = 0;
  count = 0;
  fprintf(fp, "\nEdges: (id, master)\n");
  while (ent=M_nextEdge(mesh, &temp))
    fprintf(fp, "V%-2d-V%-2d ---------> %-6d\n", EN_id(E_vertex(ent, 0)),
            EN_id(E_vertex(ent, 1)), EN_ownerProc(ent));

  temp = 0;
  count = 0;
  fprintf(fp, "\nFaces: (id, master)\n");
  while (ent=M_nextFace(mesh, &temp))
    fprintf(fp, "%6d \t%6d\n", count++, EN_ownerProc(ent));

  temp = 0;
  count = 0;
  fprintf(fp, "\nRegions: (id, master)\n");
  while (ent=M_nextRegion(mesh, &temp))
    fprintf(fp, "%6d \t%6d\n", count++, EN_ownerProc(ent));
#endif

  fclose(fp);
}
