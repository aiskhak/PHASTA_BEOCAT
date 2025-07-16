
/**********************************************************************/
/* This function gathers all the mesh entities in the closure of a    */
/* tetrahedron.                                                       */
/**********************************************************************/
#include "parallel.h"
#include "func.h"
#include "SimPMesh.h"

/* This routine takes the region and its boundary face as inputs and
   changes the  element numbering so that the boundary face comes
   first */

void R_entitiesBdry(pRegion region, pFace face, pVertex *vrts, pEdge *edgs,
                    pFace *fcs, int nen)
{
  int i,k,dir;
  pPList list,flist,list2,vlist;
  pFace face2;
  pFace f0,f1,f2,f3,f4;
  switch (nen) {
    /* tets */
  case 4:
    /* collect the four vertices */
    /* determine the direction the face is being used by the region */
    dir = R_dirUsingFace(region,face);
    list = F_vertices(face,dir);
    for(i=0; i< 3 ; i++) vrts[i] = (pVertex)PList_item(list,i);
    vrts[3] = R_fcOpVt(region,face); /* other vertex */
    PList_delete(list);

    if (edgs)  {
      /* collect the six edges */
      edgs[0] = E_exists(vrts[0],vrts[1]);
      edgs[1] = E_exists(vrts[1],vrts[2]);
      edgs[2] = E_exists(vrts[2],vrts[0]);
      edgs[3] = E_exists(vrts[0],vrts[3]);
      edgs[4] = E_exists(vrts[1],vrts[3]);
      edgs[5] = E_exists(vrts[2],vrts[3]);
    }

    if (fcs)  {
      /* collect the four faces */
      fcs[0] = face;
      fcs[1] = R_vtOpFc(region,vrts[2]);
      fcs[2] = R_vtOpFc(region,vrts[0]);
      fcs[3] = R_vtOpFc(region,vrts[1]);
    }

    break;

    /* pyramid */
  case 5:

#ifdef SUN4
      vlist = R_vertices(region,1);
#endif

#ifdef SGI
      vlist = R_vertices(region,1);
#endif


    f0 = R_face(region,0);
    f1 = R_face(region,1);
    f2 = R_face(region,2);
    f3 = R_face(region,3);
    f4 = R_face(region,4);

    if (F_numEdges(face) == 4) {
      /* base face */
      vrts[0] = (pVertex)PList_item(vlist,0);
      vrts[1] = (pVertex)PList_item(vlist,1);
      vrts[2] = (pVertex)PList_item(vlist,2);
      vrts[3] = (pVertex)PList_item(vlist,3);
      vrts[4] = (pVertex)PList_item(vlist,4);

    } else if(face == f1) {
      /* first tri face */
      vrts[0] = (pVertex)PList_item(vlist,0);
      vrts[1] = (pVertex)PList_item(vlist,1);
      vrts[2] = (pVertex)PList_item(vlist,2);
      vrts[3] = (pVertex)PList_item(vlist,3);
      vrts[4] = (pVertex)PList_item(vlist,4);

    } else if(face == f2){
      /* we need to rotate the pyramid around the zeta axis for the rest
         of the three faces */
      vrts[0] = (pVertex)PList_item(vlist,1);
      vrts[1] = (pVertex)PList_item(vlist,2);
      vrts[2] = (pVertex)PList_item(vlist,3);
      vrts[3] = (pVertex)PList_item(vlist,0);
      vrts[4] = (pVertex)PList_item(vlist,4);

    } else if(face == f3){
      vrts[0] = (pVertex)PList_item(vlist,2);
      vrts[1] = (pVertex)PList_item(vlist,3);
      vrts[2] = (pVertex)PList_item(vlist,0);
      vrts[3] = (pVertex)PList_item(vlist,1);
      vrts[4] = (pVertex)PList_item(vlist,4);

    } else if(face == f4){
      vrts[0] = (pVertex)PList_item(vlist,3);
      vrts[1] = (pVertex)PList_item(vlist,0);
      vrts[2] = (pVertex)PList_item(vlist,1);
      vrts[3] = (pVertex)PList_item(vlist,2);
      vrts[4] = (pVertex)PList_item(vlist,4);
    }

    PList_delete(vlist);

    if (edgs)  {
      edgs[0] = E_exists(vrts[0],vrts[1]);
      edgs[1] = E_exists(vrts[1],vrts[2]);
      edgs[2] = E_exists(vrts[2],vrts[3]);
      edgs[3] = E_exists(vrts[3],vrts[0]);
      edgs[4] = E_exists(vrts[0],vrts[4]);
      edgs[5] = E_exists(vrts[1],vrts[4]);
      edgs[6] = E_exists(vrts[2],vrts[4]);
      edgs[7] = E_exists(vrts[3],vrts[4]);
    }

    if (fcs)  {
      fcs[0]=F_exists(Tvertex, vrts[0], vrts[1], vrts[2], vrts[3]);
      fcs[1]=F_exists(Tvertex, vrts[0],vrts[1],vrts[4],0);
      fcs[2]=F_exists(Tvertex, vrts[1],vrts[2],vrts[4],0);
      fcs[3]=F_exists(Tvertex, vrts[2],vrts[3],vrts[4],0);
      fcs[4]=F_exists(Tvertex, vrts[3],vrts[0],vrts[4],0);
    }

    break;

    /* wedges */
  case 6:

#ifdef SUN4
      flist = R_faces(region,1);
      vlist = R_vertices(region,1);
#endif

#ifdef SGI
      flist = R_faces(region,1);
      vlist = R_vertices(region,1);
#endif


    if (face == (pFace)PList_item(flist,0)) {
      vrts[0] = (pVertex)PList_item(vlist,0);
      vrts[1] = (pVertex)PList_item(vlist,1);
      vrts[2] = (pVertex)PList_item(vlist,2);
      vrts[3] = (pVertex)PList_item(vlist,3);
      vrts[4] = (pVertex)PList_item(vlist,4);
      vrts[5] = (pVertex)PList_item(vlist,5);
    } else if(F_numEdges(face)==3) {  /* the other tri face */
      vrts[0] = (pVertex)PList_item(vlist,3);
      vrts[1] = (pVertex)PList_item(vlist,5);
      vrts[2] = (pVertex)PList_item(vlist,4);
      vrts[3] = (pVertex)PList_item(vlist,0);
      vrts[4] = (pVertex)PList_item(vlist,2);
      vrts[5] = (pVertex)PList_item(vlist,1);
    } else if(EN_inClosure(face,(pVertex)PList_item(vlist,0)) &&
              EN_inClosure(face,(pVertex)PList_item(vlist,1))) {
      vrts[0] = (pVertex)PList_item(vlist,0);
      vrts[1] = (pVertex)PList_item(vlist,1);
      vrts[2] = (pVertex)PList_item(vlist,2);
      vrts[3] = (pVertex)PList_item(vlist,3);
      vrts[4] = (pVertex)PList_item(vlist,4);
      vrts[5] = (pVertex)PList_item(vlist,5);

    } else if(EN_inClosure(face,(pVertex)PList_item(vlist,2)) &&
              EN_inClosure(face,(pVertex)PList_item(vlist,1))) {
      vrts[0] = (pVertex)PList_item(vlist,1);
      vrts[1] = (pVertex)PList_item(vlist,2);
      vrts[2] = (pVertex)PList_item(vlist,0);
      vrts[3] = (pVertex)PList_item(vlist,4);
      vrts[4] = (pVertex)PList_item(vlist,5);
      vrts[5] = (pVertex)PList_item(vlist,3);

    } else { /* contention */
      vrts[0] = (pVertex)PList_item(vlist,3);
      vrts[1] = (pVertex)PList_item(vlist,5);
      vrts[2] = (pVertex)PList_item(vlist,4);
      vrts[3] = (pVertex)PList_item(vlist,0);
      vrts[4] = (pVertex)PList_item(vlist,2);
      vrts[5] = (pVertex)PList_item(vlist,1);
    }

    PList_delete(vlist);

    PList_delete(flist);

    if (edgs)  {
      /* create the stuff for hierarchic edge and face for wedge element */
      edgs[0] = E_exists(vrts[0],vrts[1]);
      edgs[1] = E_exists(vrts[1],vrts[2]);
      edgs[2] = E_exists(vrts[2],vrts[0]);
      edgs[3] = E_exists(vrts[3],vrts[4]);
      edgs[4] = E_exists(vrts[4],vrts[5]);
      edgs[5] = E_exists(vrts[5],vrts[3]);
      edgs[6] = E_exists(vrts[0],vrts[3]);
      edgs[7] = E_exists(vrts[1],vrts[4]);
      edgs[8] = E_exists(vrts[2],vrts[5]);
    }

    if (fcs)  {
      fcs[0]=F_exists(Tvertex, vrts[0], vrts[2], vrts[1], 0);
      fcs[1]=F_exists(Tvertex, vrts[0],vrts[1],vrts[4],vrts[3]);
      fcs[2]=F_exists(Tvertex, vrts[1],vrts[2],vrts[5],vrts[4]);
      fcs[3]=F_exists(Tvertex, vrts[2],vrts[0],vrts[3],vrts[5]);
      fcs[4]=F_exists(Tvertex, vrts[3],vrts[4],vrts[5],0);
    }

    break;

  case 8:
    /* the element is a brick */
    /* get vertices on boundary face, and faces on the region */

    dir = R_dirUsingFace(region,face);

    list = F_vertices(face,1-dir);

#ifdef SUN4
    flist = R_faces(region,1);
#endif

#ifdef SGI
    flist = R_faces(region,1);
#endif

  




    for (i=0; i < 4; i++) {
      vrts[i] = PList_item(list,i);
    }

    /* find the face on the opposite side of the hex */
    for (i=0; i < PList_size(flist); i++) {
      if (PList_item(flist,i) != face &&
            !F_conToFace(face,PList_item(flist,i))) {
        face2 = PList_item(flist,i);
        break;
      }
    }

    dir = R_dirUsingFace(region,face2);
    /* find the vertex on the other face matching vertex 0 */
    list2 = F_vertices(face2,dir);
    for (k=0; k < 4; k++) {
      for (i=0; i < 4; i++) {
        if (E_exists(PList_item(list2,i),vrts[k])) {
          vrts[4+k] = PList_item(list2,i);
          break;
        }
      }
    }

    PList_delete(list);
    PList_delete(list2);

    if (!(face = F_exists(Tvertex,vrts[0], vrts[1],vrts[2],vrts[3]))){
      fprintf(stderr,"Error: face not found in getConnectivity...");
      exit(-1);
    }

    if (edgs)  {
      /* 12 edges on this region */
      edgs[0] = E_exists(vrts[0],vrts[1]);
      edgs[1] = E_exists(vrts[1],vrts[2]);
      edgs[2] = E_exists(vrts[2],vrts[3]);
      edgs[3] = E_exists(vrts[3],vrts[0]);
      edgs[4] = E_exists(vrts[4],vrts[5]);
      edgs[5] = E_exists(vrts[5],vrts[6]);
      edgs[6] = E_exists(vrts[6],vrts[7]);
      edgs[7] = E_exists(vrts[7],vrts[4]);
      edgs[8] = E_exists(vrts[0],vrts[4]);
      edgs[9] = E_exists(vrts[1],vrts[5]);
      edgs[10] = E_exists(vrts[2],vrts[6]);
      edgs[11] = E_exists(vrts[3],vrts[7]);
    }

    if (fcs)  {
      /* there are 6 faces to this region */
      fcs[0]=face;
      fcs[1]=F_exists(Tvertex, vrts[0],vrts[1],vrts[5],vrts[4]);
      fcs[2]=F_exists(Tvertex, vrts[1],vrts[2],vrts[6],vrts[5]);
      fcs[3]=F_exists(Tvertex, vrts[2],vrts[3],vrts[7],vrts[6]);
      fcs[4]=F_exists(Tvertex, vrts[3],vrts[0],vrts[4],vrts[7]);
      fcs[5]=F_exists(Tvertex, vrts[4],vrts[5],vrts[6],vrts[7]);
    }
    break;

  default:
    fprintf(stderr,"\nError: R_entitiesBdry not impl for %d verts\n",nen);
    exit(-1);
  }
}
