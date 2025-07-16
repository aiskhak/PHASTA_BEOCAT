#include "parallel.h"
#include "MeshTypes.h"

/* The following function retuns a list of edges on a mesh region */
/* not in any particular order */

pPList nsR_edges(pRegion region)
{
  pPList eList = PList_new();
  pPList eList2;
  pFace  face;
  int i,nfaces;

  nfaces = R_numFaces(region);

  for(i=0; i< nfaces; i++){  /* iterating over faces */
    face = R_face(region,i);
    eList2 = F_edges(face,1,0);
    PList_appPListUnique(eList,eList2);
        PList_delete(eList2);
  }
  return eList;
}
