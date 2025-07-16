/* attach information local to each entity */

#include <stdio.h>
#include "func.h"
#include "parallel.h"
#include "bits.h"

/*********************************************************************
* allocate and initialize a local information array to each entity   *
* in the mesh                                                        *
**********************************************************************/
// local information consists of:
// globalP for each entity
// the mesh is local to the calling proc

extern int globalP;

// only in parallel
// replaces previous MeshSim functions
// declared in SimPMesh.h at least up to
// release 5.3-040830
void PMU_commuArr(void **s, int *ns, void **r, int *nr, MPI_Datatype type,
                  int mult)
{


    int i, m, tag = 0;

    MPI_Request* req = (MPI_Request *)malloc(sizeof(MPI_Request)*2*(PMU_size()-1) ); 
    MPI_Status* stat = (MPI_Status *)malloc(sizeof(MPI_Status)*2*(PMU_size()-1) );
    
/*      MPI_Request* req = new MPI_Request[2*(PMU_size()-1)]; */
/*      MPI_Status*  stat= new MPI_Status[2*(PMU_size()-1)]; */
    
    for (m=0, i=0; i < PMU_size(); i++)  {
        if (i != PMU_rank())  {
            MPI_Irecv(r[i], mult*nr[i], type, i, tag, MPI_COMM_WORLD, &req[m++]);
            MPI_Isend(s[i], mult*ns[i], type, i, tag, MPI_COMM_WORLD, &req[m++]);
        }
    }
    MPI_Waitall(m, req, stat);

/*      delete [] req; */
/*      delete [] stat; */

    free(req);
    free(stat);
}


void PMU_commuInt(int *ns, int *nr)
{
    MPI_Request* req = (MPI_Request *)malloc(sizeof(MPI_Request)*2*(PMU_size()-1) ); 
    MPI_Status* stat = (MPI_Status *)malloc(sizeof(MPI_Status)*2*(PMU_size()-1) );

/*      MPI_Request* req = new MPI_Request[2*(PMU_size()-1)]; */
/*      MPI_Status*  stat= new MPI_Status[2*(PMU_size()-1)]; */
    int i, m, tag = 0;
    for (m=0, i=0; i < PMU_size(); i++)  {
        if (i != PMU_rank())  {
            MPI_Irecv(nr+i, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &req[m++]);
            MPI_Isend(ns+i, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &req[m++]);
        }
    }
    MPI_Waitall(m, req, stat);
    
/*      delete [] req; */
/*      delete [] stat;     */
    free(req);
    free(stat);
}

void PMU_commuDouble(double* ds, double* dr)
{
    MPI_Request* req = (MPI_Request *)malloc(sizeof(MPI_Request)*2*(PMU_size()-1) ); 
    MPI_Status* stat = (MPI_Status *)malloc(sizeof(MPI_Status)*2*(PMU_size()-1) );

/*      MPI_Request* req = new MPI_Request[2*(PMU_size()-1)]; */
/*      MPI_Status*  stat= new MPI_Status[2*(PMU_size()-1)]; */
    int i, m, tag = 0;
    for (m=0, i=0; i < PMU_size(); i++)  {
        if (i != PMU_rank())  {
            MPI_Irecv(dr+i, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &req[m++]);
            MPI_Isend(ds+i, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &req[m++]);
        }
    }
    MPI_Waitall(m, req, stat);
    
/*      delete [] req; */
/*      delete [] stat;     */
    free(req);
    free(stat);
}


// called in initLocalInfo (below)
int periocheck(pEntity ent, globalInfo *info)
{
  pMatch ment;

  // isPeriodic declared in setPeriodic.cc
  if (isPeriodic(ent, &ment))  {
    if (isOnThisProc(ment))
      return 1;
    PList_append(info->perSlv, ent);
    PList_append(info->perMst, ment);
    return 0;
  }
  return 2;
}



//////////////////////////////////////////////////////////////////////////////
// called in mdb2phasta
//////////////////////////////////////////////////////////////////////////////
#ifdef PARALLEL
void initLocalInfo(pParMesh pmesh, globalInfo *ginfo)
#else
void initLocalInfo(pMesh mesh, globalInfo* ginfo)
#endif
{
#ifdef PARALLEL
  pMesh mesh = PM_mesh(pmesh,0);
  void *tmps, *tmpm;
#endif
  pMatch ment;
  pVertex vertex;
  pEdge   edge;
  pFace   face;
  pRegion region;
  pEntity ent;
  VIter vIter;
  EIter eIter;
  FIter fIter;
  RIter rIter;

  int Biggest,i,j,tbit;
  int localP, nen;
  int p; /* typing saver for global/localP */
  int ndof;
  int mybig;
  void *tmp;

  int gid, pid, count = 0, pc;
  pPList myents = PList_new();
#ifdef PARALLEL
  pEntCopies ec;
  int *eqst = (int*) malloc (sizeof(int)*PMU_size());
  int rst, mult = sizeof(P_int)/sizeof(int);
  int *ns = (int*) malloc(sizeof(int)*PMU_size());
  int *nr = (int*) malloc(sizeof(int)*PMU_size());
  int **nsb = (int**) malloc(sizeof(int*)*PMU_size());  /* send back */
  int **nrb = (int**) malloc(sizeof(int*)*PMU_size());  /* recv back */
  void ***mye = (void***) malloc(sizeof(void**)*PMU_size());
  void ***se = (void***) malloc(sizeof(void**)*PMU_size());
  void ***re = (void***) malloc(sizeof(void**)*PMU_size());
  for (i=0; i < PMU_size(); i++)  {
    ns[i] = PM_numBdryEntOwner(pmesh,0,0,i) + PM_numBdryEntOwner(pmesh,0,1,i)
            + PM_numBdryEntOwner(pmesh,0,2,i);
    se[i] = (void**) malloc(sizeof(void*)*2*ns[i]);
    ns[i] = nr[i] = 0;
  }
#endif

  /* initializing the edge /face/ region mode info */
  ginfo->edgeson =0;
  ginfo->faceson =0;
  ginfo->regnson =0;
  ginfo->nshg = 0;
  ginfo->nshgOwn = 0;

  /*******************************************************************/
  /* allocate and attach a local information structure to each       */
  /* entity in the mesh                                              */
//
// what exactly is done here?
  /*******************************************************************/
  /* vertices */
  vIter = M_vertexIter(mesh);
  while (vertex = VIter_next(vIter))  {
    /* setting the number of dofs based on localP */
    localP  = globalP;
    EN_attachDataI(vertex,"POLY",localP);
    EN_attachDataI(vertex,"NDOF",1);
    // MYCT = mycount: local equation number as it occurs in the local linear
    // system 
    // NDOF of specific entity (p=3 ==>NDOF=2)
    // 
    EN_attachDataI(vertex,"MYCT",count++);
    ginfo->nshg++;

    // 
    if ((pc = periocheck(vertex, ginfo)) && EN_isOwnerProc(vertex))  {

        // global  equation number for vertices
      EN_attachDataI(vertex,"EQST", pc==1 ? -1 : ginfo->nshgOwn);


      PList_append(myents, vertex);
#ifdef PARALLEL
      // copies of vertex on partition boundary
      ec = EN_copies(vertex);

      // loop over the multitude of the vertex' copies
      // i.e. all the slaves
      for (i=0; i < EntCopies_size(ec); i++)  {

          // Returns global id of partition on which ec's i'th copy lies.
        gid = EntCopies_gid(ec, i);
        if ((pid = PMU_proc(gid)) != PMU_rank())  {
          se[pid][2*ns[pid]] = EntCopies_ent(ec, i);
          se[pid][2*ns[pid]+1] = (void*) (pc==1 ? -1 : ginfo->nshgOwn);
          ns[pid]++;
        }
      }//for
#endif// PARALLEL
      if (pc != 1){
        ginfo->nshgOwn++;
      }
    }// if ((pc ... periodic stuff
  }//while (vertex
  VIter_delete(vIter);


  /* edges */
  eIter = M_edgeIter(mesh);
  while (edge = EIter_next(eIter))  {

    localP  = globalP;  /* P-change */
    ndof    = (localP) - 1;

    if (localP > 1) ginfo->edgeson=1;// edgemodes only for globalP=2 and higher

    EN_attachDataI(edge,"POLY",localP);

    if (ndof) {

      EN_attachDataI(edge,"NDOF",ndof);
      EN_attachDataI(edge,"MYCT",count);  count += ndof;
      ginfo->nshg += ndof;

      if ((pc = periocheck(edge,ginfo)) && EN_isOwnerProc(edge))  {
        EN_attachDataI(edge,"EQST",ginfo->nshgOwn);
        PList_append(myents, edge);

#ifdef PARALLEL
        ec = EN_copies(edge);
        for (i=0; i < EntCopies_size(ec); i++)  {

          gid = EntCopies_gid(ec, i);
          if ((pid = PMU_proc(gid)) != PMU_rank())  {

            se[pid][2*ns[pid]] = EntCopies_ent(ec, i);
            se[pid][2*ns[pid]+1] = (void*) (pc==1 ? -1 : ginfo->nshgOwn);
            ns[pid]++;
          }
        }
#endif
        if (pc != 1)
          ginfo->nshgOwn += ndof;
      }//periocheck
    }//ndof!=0
  }
  EIter_delete(eIter);

  /* faces */
  fIter = M_faceIter(mesh);
  while (face = FIter_next(fIter))  {
    localP  = globalP; /* P-change */
    switch(F_numEdges(face)){
    case 3:              /* Tri Face */
      if (localP > 2) {
          ginfo->faceson = 1; // facemodes only for globalP=3 and higher
        ndof = ((localP - 1)*(localP - 2))/2;
      } else {
        ndof = 0;
      }
      break;
    case 4:            /* Quad Face */
      if (localP > 3) {
        ginfo->faceson = 1;
        ndof = ((localP - 3)*(localP - 2))/2;
      } else {
        ndof = 0;
      }
      break;
    default:
      fprintf(stderr,"Face with neither 3 nor 4 edges \n");
      exit(1);
    }
    EN_attachDataI(face,"POLY",localP);
    if (ndof) {
      EN_attachDataI(face,"NDOF",ndof);
      EN_attachDataI(face,"MYCT",count);  count += ndof;
      if ((pc = periocheck(face,ginfo)) && EN_isOwnerProc(face))  {
        EN_attachDataI(face,"EQST",ginfo->nshgOwn);
        PList_append(myents, face);
#ifdef PARALLEL
        ec = EN_copies(face);
        for (i=0; i < EntCopies_size(ec); i++)  {
          gid = EntCopies_gid(ec, i);
          if ((pid = PMU_proc(gid)) != PMU_rank())  {
            se[pid][2*ns[pid]] = EntCopies_ent(ec, i);
            se[pid][2*ns[pid]+1] = (void*) (pc==1 ? -1 : ginfo->nshgOwn);
            ns[pid]++;
          }
        }
#endif
        if (pc != 1)
          ginfo->nshgOwn += ndof;
      }
    }
  }
  FIter_delete(fIter);

  /* regions */
  Biggest = 0;
  rIter = M_regionIter(mesh);
  while (region =RIter_next(rIter))  {
    localP  = globalP;  /* P-change, has to be replaced with a
                         function which returns the maximum polynomial
                          order of the all the enclosed entities */
    switch(topology(region)){
    case 1:    /* tets */
      nen  = 4;
      Biggest = setbit(Biggest, 0);
      if ((p = localP) > 3)// region modes  only for globalP=4 and higher
        ndof = (p-1)*(p-2)*(p-3)/6;
      else ndof = 0;
      break;
    case 5:    /* Pyramids */
      nen  = 5;
      Biggest = setbit(Biggest, 1);
      if ((p = localP) > 5)
        ndof = (p-3)*(p-4)*(p-5)/6;
      else ndof = 0;
      break;
    case 3:   /* Wedges */
      nen = 6;
      Biggest = setbit(Biggest, 2);
      if ((p = localP) > 4)
        ndof = (p-2)*(p-3)*(p-4)/6;
      else ndof = 0;
      break;
    case 2:   /* Hexes */
      nen = 8;
      Biggest = setbit(Biggest, 3);
      if ((p = localP) > 5)
        ndof =(p-3)*(p-4)*(p-5)/6;
      else ndof = 0;
      break;
    default:
      fprintf(stderr,"Unknown Element type in LocalInfo\n");
      exit(-1);
    }
    EN_attachDataI(region,"RNEN",nen);
    if (ndof) ginfo->regnson = 1;
    EN_attachDataI(region,"POLY",localP);
    if (ndof) {
      EN_attachDataI(region,"EQST",ginfo->nshgOwn);
      PList_append(myents, region);
      EN_attachDataI(region,"NDOF",ndof);
      EN_attachDataI(region,"MYCT",count);  count += ndof;
      ginfo->nshg += ndof;
      ginfo->nshgOwn += ndof;  /* regions always owned */
    }
  }
  RIter_delete(rIter);

#ifdef PARALLEL
  mybig = Biggest;
  MPI_Allreduce(&mybig, &Biggest, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  tbit = 0;
  i = 3;
  while (!(tbit = getbit(Biggest,i--)) && (i >= 0));
  if (!tbit) {
    fprintf(stderr," Bit test for Topology has problems \n");
    fprintf(stderr," Please stop in func initlocalinfo for debugging \n");
    exit(1);
  } else { tbit = ++i ;}

  switch(tbit){
  case 0:    /* tets */
    ginfo->nen  = 4;
    ginfo->nenb = 3;
    ginfo->nedges = 6;
    ginfo->nfaces = 4;
    break;
  case 1:    /* Pyramids */
    ginfo->nen  = 5;
    ginfo->nenb = 4; /* should not matter, we never let Pyramids
                        get on the boundaries */
    ginfo->nedges = 8;
    ginfo->nfaces = 5;
    break;
  case 2:   /* Wedges */
    ginfo->nen = 6;
    ginfo->nenb = 4;  /* wedge can also have Tri Boundary face,
                         keep this in mind */
    ginfo->nedges = 9;
    ginfo->nfaces = 5;
    break;
  case 3:   /* Hexes */
    ginfo->nen = 8;
    ginfo->nenb = 4;
    ginfo->nedges = 12;
    ginfo->nfaces = 6;
    break;
  default:
    fprintf(stderr,"Unknown Element type in LocalInfo\n");
    exit(-1);
  }
#ifdef PARALLEL


  // now for the communication 
  // needs documentation -JM
  MPI_Scan(&ginfo->nshgOwn, eqst+PMU_rank(), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  eqst[PMU_rank()] -= ginfo->nshgOwn;
  MPI_Allgather(eqst+PMU_rank(), 1, MPI_INT, eqst, 1, MPI_INT, MPI_COMM_WORLD);
  /* so if ginfo->nshgOwn(p0,p1,p2) = 10,12,8; now eqst[0..2] = 0,10,22 */

  PMU_commuInt(ns, nr);
  for (i=0; i < PMU_size(); i++)  {
    re[i] = (void**) malloc(sizeof(void*)*2*nr[i]);
    ns[i] = 0;
  }

  if (eqst[PMU_rank()] > 0)  {
    for (tmp=0; ent = (pEntity) PList_next(myents, &tmp);)
      EN_modifyDataI(ent, "EQST", EN_dataI(ent,"EQST") + eqst[PMU_rank()]);
  }
#endif

  for (tmp=0; ent = (pEntity) PList_next(myents, &tmp);)  {
    if (isPeriodic(ent, &ment) && isOnThisProc(ment))  {

       pc = EN_dataI(Match_ent(ment), "EQST");
 
      EN_modifyDataI(ent, "EQST", pc);
    }
#ifdef PARALLEL
    ec = EN_copies(ent);
    for (i=0; i < EntCopies_size(ec); i++)  {
      gid = EntCopies_gid(ec, i);
      if ((pid = PMU_proc(gid)) != PMU_rank())  {
        if (se[pid][2*ns[pid]+1] == (void*) -1)
          se[pid][2*ns[pid]+1] = (void*) pc;
        ns[pid]++;
      }
    }
#endif
  }

#ifdef PARALLEL
  PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, 2*mult);

  for (i=0; i < PMU_size(); i++)  {
    for (j=0; j < nr[i]; j++)  {
      ent = (pEntity) re[i][2*j];
      rst = (P_int) re[i][2*j+1];
      EN_attachDataI(ent, "EQST", rst + eqst[i]);
    }

    /* Following for setting periodic slave EQST's right */
    free(se[i]);
    free(re[i]);
    ns[i] = (i==PMU_rank()) ? 0 : PList_size(ginfo->perSlv);
    mye[i] = (void**) malloc(sizeof(void*)*ns[i]);
    se[i] = (void**) malloc(sizeof(void*)*ns[i]);
    ns[i] = nr[i] = 0;
  }

  for (tmps=0, tmpm=0; (ent = (pEntity) PList_next(ginfo->perSlv, &tmps))
                    && (ment = (pMatch) PList_next(ginfo->perMst, &tmpm));)  {

       pid = PMU_proc(Match_gid(ment));  /* We are sure ment is off-proc */
       mye[pid][ns[pid]] = ent;

       se[pid][ns[pid]] = Match_ent(ment);
       ns[pid]++;
  }

  PMU_commuInt(ns, nr);
  for (i=0; i < PMU_size(); i++)  {
    re[i] = (void**) malloc(sizeof(void*)*nr[i]);
    nsb[i] = (int*) malloc(sizeof(int)*nr[i]);
    nrb[i] = (int*) malloc(sizeof(int)*ns[i]);
  }
  PMU_commuArr((void**)se, ns, (void**)re, nr, MPI_INT, mult);

  for (i=0; i < PMU_size(); i++)  {
    for (j=0; j < nr[i]; j++)
      nsb[i][j] = EN_dataI((pEntity) re[i][j], "EQST");
  }
  PMU_commuArr((void**)nsb, nr, (void**)nrb, ns, MPI_INT, 1);

  for (i=0; i < PMU_size(); i++)  {
    for (j=0; j < ns[i]; j++)
      EN_attachDataI((pEntity) mye[i][j], "EQST", nrb[i][j]);
    free(mye[i]);
    free(se[i]);
    free(re[i]);
    free(nsb[i]);
    free(nrb[i]);
  }

  free(eqst);
  free(mye);
  free(se);
  free(re);
  free(nsb);
  free(nrb);
  free(ns);
  free(nr);
#endif
  PList_delete(myents);
}
