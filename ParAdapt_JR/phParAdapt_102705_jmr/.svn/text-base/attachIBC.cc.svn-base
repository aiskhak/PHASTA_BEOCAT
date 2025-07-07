////////////////////////////////////////////////////////////////////////
//
// These functions set up the boundary condition information.
//
////////////////////////////////////////////////////////////////////////
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include <map>
#ifndef SIM
#include "modeler.h"
#else
#include "MeshSim.h"
#endif
#include "func.h"
#include "ccfunc.h"
#include "parallel.h"
#include "EssentialBC.h"
#include "NaturalBC.h"
#include "InitialCondition.h"

#define RGAS 288.294382801664369

// global variables
extern int rStart,intBC,prCd,ensa_dof,rRead;
extern int multipleRestarts;
extern forwardblock Bblock;
extern int adaptFlag;
extern "C" int topology(pRegion region);

int attachEssentialBC(pMesh mesh, globalInfo *info, double *qTot, int *nBC, int *iBC, double **BC)
{
  pGModel model = M_model(mesh);
  pGEntity gent;
  double   *intlist;
  int code,numpbc=0,sign;
  int GDOFnum; // DOF number in qTot 
  int nshg = info->nshg;
  int nshgTot = info->nshgTot;
  ofstream fout;
  if (prCd && PMU_rank()==0)
    fout.open("BC.out", ios::app);
  int done;

  // if there are internal boundary nodes (on a trip strip)
  // get their values now
  if (intBC == 1)  {
    GRIter grIter = GM_regionIter(model);
    gent    = (pGEntity)GRIter_next(grIter);
    GRIter_delete(grIter);
    intlist = (double *)GEN_dataP(gent,"bc  ");

    VIter vIter = M_vertexIter(mesh);
    pVertex vtx;
    while (vtx = VIter_next(vIter)) {
      code = EN_dataI((pEntity)vtx,"ibc ");
      if (code > 0) {
        nBC[EN_dataI((pEntity)vtx, "MYCT")] = numpbc;
        iBC[numpbc] = code;
        for (int i=0; i < 11; i++) {
          BC[numpbc][i] = intlist[i];
        }
        numpbc++;
      }
    }
    VIter_delete(vIter);
  }

  if (PMU_rank()==0 && prCd)  {
    fout << "\n\nEssential boundary condition codes\n";
    fout << "----------------------------------\n\n";
    fout << "Entity" << "\tiBC-Code\n";
    fout << "------" << "\t--------\n";
  }

  // model faces
  GFIter gfIter = GM_faceIter(model);
  pGFace gface;
  while(gface = GFIter_next(gfIter))  {

    EssentialBC bc(gface);

    double theta =  *(double*)(GEN_dataP((pGEntity)gface,"teta"));

    if (PMU_rank()==0 && prCd)  done = 0;
    if (bc.isSet() || GEN_dataI((pGEntity)gface,"iper")
                   || GEN_dataI((pGEntity)gface,"sper"))  {

      // mesh vertices on this face
      VIter vIter  = M_classifiedVertexIter(mesh, (pGEntity)gface, 0);
      pVertex vert;
      while (vert = VIter_next(vIter))  {
        if (EN_dataI((pEntity)vert,"ibc ") > 0) continue;
                                                 // nodal bc code attached

        if (bc.isSet())  {
//            printf("\n[%d]in attachIBC: attachEssentialBC :bc.isSet() \n", PMU_rank() );
            if ((code = bc.eval(vert,BC[numpbc])) != 0)  {
                // here is a function that will override the BC values
                // obtained by evaluating the expressions (in the
                // previous line) with those given in the initial
                // condition (typically read in solution from a previous
                // mesh) for entities with "take bc from ic" attribute set
	        if(!EN_getDataInt((pEntity)vert,info->incorp,&GDOFnum)) {
		  GDOFnum=EN_dataI((pEntity)vert, "MYCT");
	        }
                // bc.takeBCfromIC(BC[numpbc],qTot,nshgTot,GDOFnum);
                bc.takeBCfromIC(BC[numpbc],qTot,nshg,GDOFnum);

                nBC[EN_dataI((pEntity)vert, "MYCT")] = numpbc;
                iBC[numpbc] = code;
                
                sign=((EN_dataI((pEntity)vert,"swit"))? -1 : 1);
                BC[numpbc][11] = sign*theta;
                
                if (PMU_rank()==0 && prCd && !done)  {
                    fout << "F-" << GEN_tag((pGEntity)gface)<< "\t\t" ;
                    fout << code << "\n";
                    done = 1;
                }
                numpbc++;
//                printf("\n[%d]in attachIBC: attachEssentialBC  (bc.isSet incr. numpbc to %d: \n", PMU_rank(),numpbc );
            }//if ((code
        }//if (bc.isSet()) 
        // NOT (bc.isSet()
        else {
          code =0;
          if (GEN_dataI((pGEntity)gface,"iper")) {
            code = 1024;
          }
          sign=((EN_dataI((pEntity)vert,"swit"))? -1 : 1);
          BC[numpbc][11] = sign*theta;
          if (GEN_dataI((pGEntity)gface,"sper"))  code += 2048;
          nBC[EN_dataI((pEntity)vert, "MYCT")] = numpbc;
          iBC[numpbc] = code;

          if (PMU_rank()==0 && prCd && !done)  {
            fout << "F-" << GEN_tag((pGEntity)gface)  << "\t\t" ;
            fout << code << "\n";
            done = 1;
          }
          numpbc++;
//          printf("\n[%d]in attachIBC: attachEssentialBC NOT (bc.isSet incr. numpbc to %d: \n", PMU_rank(),numpbc );
        }//else: NOT (bc.isSet()
      }//while{vert=..
      VIter_delete(vIter);

      // mesh edges
      if (info->edgeson){
        EIter eIter = M_classifiedEdgeIter(mesh,(pGEntity)gface, 0);
        pEdge edge;
        int nem;
        while (edge = EIter_next(eIter))  {
          if(nem = EN_dataI((pEntity)edge, "NDOF")) {
            for (int i=0; i < nem; i++)  {
              if (bc.isSet())  {
                if ((code = bc.eval(edge,BC[numpbc],i+2)) != 0)  {
                  nBC[EN_dataI((pEntity)edge, "MYCT") + i] = numpbc;
                  iBC[numpbc] = code;
                  sign=((EN_dataI((pEntity)edge,"swit"))? -1 : 1);
                  BC[numpbc][11] = sign*theta;
                  numpbc++;
                }
              }
              else  {
                nBC[EN_dataI((pEntity)edge, "MYCT") + i] = numpbc;
                code=0;
                sign=((EN_dataI((pEntity)edge,"swit"))? -1 : 1);
                BC[numpbc][11] = sign*theta;
                if (GEN_dataI((pGEntity)gface,"iper")) {
                  code = 1024;
                }
                if (GEN_dataI((pGEntity)gface,"sper"))  code += 2048;

                iBC[numpbc] = code;
                numpbc++;
              }
            }
          }
        }
        EIter_delete(eIter);
      }
      // mesh faces
      if(info->faceson){
        FIter fIter  = M_classifiedFaceIter(mesh,(pGEntity)gface, 0);
        pFace face;
        int nfm;
        while (face = FIter_next(fIter))  {
          if(nfm = EN_dataI((pEntity)face, "NDOF")) {
            for (int i=0; i < nfm; i++)  {
              if (bc.isSet())  {
                if ((code = bc.eval(face,BC[numpbc],i+3)) != 0)  {
                  nBC[EN_dataI((pEntity)face, "MYCT") + i] = numpbc;
                  iBC[numpbc] = code;
                  numpbc++;
                }
              }
              else  {
                nBC[EN_dataI((pEntity)face, "MYCT") + i] = numpbc;
                code=0;
                sign=((EN_dataI((pEntity)face,"swit"))? -1 : 1);
                BC[numpbc][11] = sign*theta;
                if (GEN_dataI((pGEntity)gface,"iper")) {
                  code = 1024;
                }
                iBC[numpbc] = code;

                numpbc++;
              }
            }
          }
        }
	FIter_delete(fIter);
      }
    }//if (bc.isSet()
  }//while{ gface= ..
  GFIter_delete(gfIter);

  // model edges
  GEIter geIter = GM_edgeIter(model);
  pGEdge gedge;
  while(gedge = GEIter_next(geIter))  {

    // axisymmetric bc's
    double theta = *(double*)(GEN_dataP((pGEntity)gedge,"teta"));

    EssentialBC bc(gedge);

    if (PMU_rank()==0 && prCd)  done = 0;
    if (bc.isSet())  {
      // mesh vertices

      VIter vIter  = M_classifiedVertexIter(mesh,(pGEntity)gedge, 0);
      pVertex vert;

      while (vert = VIter_next(vIter))  {
        if (EN_dataI((pEntity)vert,"ibc ") > 0) continue;
        if ((code = bc.eval(vert,BC[numpbc])) != 0)  {
          nBC[EN_dataI((pEntity)vert, "MYCT")] = numpbc;
          iBC[numpbc] = code;
          sign=((EN_dataI((pEntity)vert,"swit"))? -1 : 1);
          BC[numpbc][11] = sign*theta;
          if (PMU_rank()==0 && prCd && !done)  {
            fout << "E-" << GEN_tag((pGEntity)gedge) << "\t\t" << code << "\n";
            done = 1;
          }
          numpbc++;
        }
      }
      VIter_delete(vIter);
      // mesh edges
      if (info->edgeson){
        EIter eIter  = M_classifiedEdgeIter(mesh, (pGEntity) gedge, 0);
        pEdge edge;
        int nem;
        while (edge = EIter_next(eIter))  {
          if(nem = EN_dataI((pEntity)edge, "NDOF"))  {
            for (int i=0; i < nem; i++)  {
              if ((code = bc.eval(edge,BC[numpbc],i+2)) != 0)  {
                nBC[EN_dataI((pEntity)edge, "MYCT") + i] = numpbc;
                iBC[numpbc] = code;
                sign=((EN_dataI((pEntity)edge,"swit"))? -1 : 1);
                BC[numpbc][11] = sign*theta;
                numpbc++;
              }
            }
          }
        }
        EIter_delete(eIter);
      }
    }
  }
  GEIter_delete(geIter);

  // model vertices
  GVIter gvIter = GM_vertexIter(model);
  pGVertex gvert;
  while(gvert = GVIter_next(gvIter)) {
    EssentialBC bc(gvert);

    // axisymmetric bc's
    double theta =  *(double*)(GEN_dataP((pGEntity)gvert,"teta"));

    if (PMU_rank()==0 && prCd)  done = 0;
    if (bc.isSet())  {

      // mesh vertex
      pVertex vert;
      vert = M_classifiedVertex(mesh, gvert);
      if (!vert || EN_dataI((pEntity)vert,"ibc ") > 0) continue;

      if ((code = bc.eval(vert,BC[numpbc])) != 0)  {
        nBC[EN_dataI((pEntity)vert, "MYCT")] = numpbc;
        iBC[numpbc] = code;
        sign=((EN_dataI((pEntity)vert,"swit"))? -1 : 1);
        BC[numpbc][11] = sign*theta;

        if (PMU_rank()==0 && prCd && !done)  {
          fout << "V-" << GEN_tag((pGEntity)gvert) << "\t\t" << code << "\n";
          done = 1;
        }
        numpbc++;
      }
    }
  }
  GVIter_delete(gvIter);

  GFIter gftr2 = GM_faceIter(model);
  while (gface = GFIter_next(gftr2)){// all model faces
    gent = (pGEntity) gface;
    if(GEN_dataP(gent,"teta")) {
      double* thetat = (double*) GEN_dataP(gent,"teta");
      delete thetat;
    }
  }
  GFIter_delete(gftr2);

  GEIter getr2 = GM_edgeIter(model);
  while (gedge = GEIter_next(getr2)){// all model edges
    gent = (pGEntity) gedge;
    if(GEN_dataP(gent,"teta")) {
      double* thetat = (double*) GEN_dataP(gent,"teta");
      delete thetat;
    }
  }
  GEIter_delete(getr2);

  GVIter gvtr2 = GM_vertexIter(model);
  while (gvert = GVIter_next(gvtr2)){// all model vertices
    gent = (pGEntity) gvert;
    if(GEN_dataP(gent,"teta")) {
      double* thetat = (double*) GEN_dataP(gent,"teta");
      delete thetat;
    }
  }
  GVIter_delete(gvtr2);

  if(prCd && PMU_rank()==0)
    fout.close();
  return numpbc;
}

/**********************************************************************/
/* Natural boundary conditions                                        */
/**********************************************************************/
void attachNaturalBC(pMesh mesh, int ***iBCB, double ***BCB, int numNBC)
{
  blockKey BLOCK;
  int blockid;
  std::map<int, int> iel;
  pPList ents;
  pRegion region;
  pGModel model = M_model(mesh);

  int code;
  int done;
  ofstream fout;

  if (PMU_rank()==0 && prCd) {
    fout.open("BC.out", ios::app);
    fout << "\n\nNatural boundary condition codes\n";
    fout <<     "--------------------------------\n\n";
    fout << "Entity" << "\tiBCB-Code\n";
    fout << "------" << "\t--------\n";
  }

  GFIter gfIter = GM_faceIter(model);
  pGFace gface;
  while(gface = GFIter_next(gfIter))  {

    NaturalBC bc(gface);

    if (PMU_rank()==0 && prCd)  done = 0;
    if (bc.isSet())  {                // if code is nonzero
      FIter fIter  = M_classifiedFaceIter(mesh,(pGEntity)gface,0);
      pFace face;
      while (face = FIter_next(fIter))  {
        ents = F_regions(face);
        region = (MRegion *)PList_item(ents,0);
        PList_delete(ents);
        BLOCK.nen = EN_dataI((pEntity)region,"RNEN");
        BLOCK.maxpoly = EN_dataI((pEntity)region,"POLY");
        BLOCK.nenbl = F_numEdges(face);
        BLOCK.lcsyst = topology(region);
        if(3 == BLOCK.lcsyst) BLOCK.lcsyst = BLOCK.nenbl;
        if(5 == BLOCK.lcsyst)
                      if (3 == BLOCK.nenbl) BLOCK.lcsyst = 6;
        blockid = Bblock[BLOCK]-1;
        bc.eval(face,BCB[blockid][iel[blockid]],iBCB[blockid][iel[blockid]],
                numNBC);
        code=iBCB[blockid][iel[blockid]][1];
        if (PMU_rank()==0 && prCd && !done)  {
          fout << "F-" << GEN_tag((pGEntity)gface) << "\t\t" << code << "\n";
          done = 1;
        }
        iel[blockid]++;
      }
      FIter_delete(fIter);
    }
  }
  GFIter_delete(gfIter);

   if(prCd && PMU_rank()==0)
     fout.close();

  iel.clear();
}


////////////////////////////////////////////////////////////////////////
// Initial conditions
// called in writeEnsaFiles()
//
// ICs are either retrieved from model attribute file 
// or are specified via solution-restart file(s)  
//                                                 
////////////////////////////////////////////////////////////////////////
void attachInitialCondition(pMesh mesh, globalInfo *info, double
                            *qTot, double **q)
{
 
#if  ( defined  DEBUG )
    if (PMU_rank()==0){
//        cout<<"\nentering attachInitialCondition\n";
    }
#endif
    
  pGModel model = M_model(mesh);
  pVertex  vertex;
  pEdge    edge;
  pFace    face;
  int      i,j,numnp=M_numVertices(mesh),nv,nem=0,nfm=0;
  nv=ensa_dof;//=5
  int count=0;
  int iOnProc,iTot;


  // user is supplying a restart file to use as the initial condition
  // only ONE restart file
  if (rStart!=0 && multipleRestarts==0 && adaptFlag == 0)  {
#if  ( defined  DEBUG )
    if (PMU_rank()==0){
        cout<<"\nproviding a single restart file for ICs\n";
    }
#endif
      
    int nshgTot=info->nshgTot;

    // vertex modes
    VIter vIter = M_vertexIter(mesh);
    while (vertex = VIter_next(vIter))  {

      // get count in Total mesh and on proc number to map restart

      EN_getDataInt((pEntity)vertex,info->incorp,&iTot);
      iOnProc=EN_dataI((pEntity)vertex, "MYCT");
      
      for (i=0; i < nv; i++)  {
        if (rStart == 2 && i==0)
          /* convert from density to pressure */
          q[iOnProc][0] = qTot[            iTot]*RGAS*
                          qTot[4*nshgTot + iTot];
        else
          q[iOnProc][i] = qTot[i*nshgTot+  iTot];
      }
    }
    VIter_delete(vIter);



    // restart file contains higher order modes
    /*  The following code was butchered by KEJ (well it was left
        mortally wounded by ST and I can't save it (KEJ).

    if (rStart == 4) { // only quadratic in restart file

      count = M_numVertices(mesh);
      if (info->edgeson){
        EIter eIter = M_edgeIter(mesh);
        while (edge = EIter_next(eIter))  {
          for (i=0; i < nv; i++)  {
            q[EN_dataI((pEntity)edge, "MYCT")][i] = qTot[count*nv+i];
          }
          count++;
        }
        EIter_delete(eIter);
      }
      else count++;
      
    } else if (rStart != 3) {
      count = M_numVertices(mesh);

      // edge modes
      if (info->edgeson){
        EIter eIter = M_edgeIter(mesh);
        while (edge = EIter_next(eIter))  {
          nem = EN_dataI(edge, "NDOF");
          for (j=0; j < nem; j++)  {
            for (i=0; i < nv; i++)  {
              q[EN_dataI(edge, "MYCT") + j][i] = qTot[count*nv+i];
            }
            count++;
          }
        }
        EIter_delete(eIter);
        // face modes
        if(info->faceson){
          FIter fIter = M_faceIter(mesh);
          while (face = FIter_next(fIter))  {
            nfm = EN_dataI(face, "NDOF");
            for (j=0; j < nfm; j++)  {
              for (i=0; i < nv; i++)
                q[EN_dataI(face, "MYCT") + j][i] = qTot[count*nv+i];
              count++;
            }
          }
          FIter_delete(fIter);
          // region modes
          if(info->regnson){
            RIter rIter = M_regionIter(mesh);
            while (region =RIter_next(rIter))  {
              nrm = EN_dataI(region, "NDOF");
              for (j=0; j < nrm; j++) {
                for (i=0; i < nv; i++)
                  q[EN_dataI(region, "MYCT") + j][i] = qTot[count*nv+i];
                count++;
              }
              else for (j=0; j < nrm; j++) count++;
            }
            RIter_delete(rIter);
          }
        }
      }
     stop the butcher*/

  }//if(rStart && multipleRestarts==0 && adaptFlag == 0)

  // multiple restarts are specified
  // the solution remains local on this calling proc's mesh
  else if((rStart!=0 && adaptFlag==1) || (rStart!=0 && multipleRestarts==1) ){
//  else if(rStart!=0 && adaptFlag==1){
      // break up the 1D array qTot in 
      for(j=0; j < info->nshg; j++){
            // loop num shapefuns
          for(i=0; i < nv; i++){

              // qtot is from restart() 
              q[j][i] = qTot[i*info->nshg+j];

#ifdef DEBUG
//                      cout<<"["<<PMU_rank()<<"]: attachInitialCondition - reading multiple restarts."
//                          <<" Local solu: q["<<j<<"]["<<i<<"]="<<qTot[i*info->nshg+j]<<"\n";
#endif

          }
      }
  }//else if((rStart!=0 && adaptFlag==1) || (rStart!=0 && multipleRestarts==1) )



  // intitial condition is generated from the attribute manager
  // find the initial condition attribute expression --> will skip the
  // attribute manager if no ICs have been set in there
  InitialCondition ic(model);

  // q's index is the first rank here 
  // q got overridden in case IC were read from restart files
  VIter vIter = M_vertexIter(mesh);
  while (vertex = VIter_next(vIter)){

      // q[l] has length nv (=5)
      ic.eval(vertex,q[EN_dataI((pEntity)vertex, "MYCT")]);
  }
  VIter_delete(vIter);

  // edge modes
  if (info->edgeson){
    EIter eIter = M_edgeIter(mesh);
    while (edge = EIter_next(eIter))  {
      nem = EN_dataI((pEntity)edge, "NDOF");
      if (nem > 0)  {
        for (j=0; j < nem; j++)  {
          for (i=0; i < nv; i++)  {
            ic.eval(edge,q[EN_dataI((pEntity)edge,"MYCT") + j],j+2);
          }
        }
      }
    }
    EIter_delete(eIter);
  }

  // face modes
  if (info->faceson){
    FIter fIter = M_faceIter(mesh);
    while (face = FIter_next(fIter))  {
      nfm = EN_dataI((pEntity)face, "NDOF");
      if (nfm > 0)  {
        for(j=0;j<nfm;j++) {
          for (i=0; i < nv; i++)  {
            q[EN_dataI((pEntity)face,"MYCT") + j][i]=0.0;
          }
        }
      }
    }
    FIter_delete(fIter);
  }
}
