#include "M2MSolutionTransfer.h"

#include "EN_isPointInside.h"

#include "MeshTypes.h"
#include "SimMeshTools.h"

#include <iostream>
#include <cmath>
#include <stdlib.h>

using std::cout;
using std::endl;
using std::cin;

double distSqTol = 1.e-12; // 1.e-6*1.e-6;
double edgeLenSqFactor = 25.;

void GetAndSetEntTagsInfo(pGModel model, pMesh mesh, pMeshDataId entTaggedId, int nEntTags, int *entTagsInfo) {
  // each unit contains entity-tag, entity-dimension and closure-flag

  cout<<"\n Enter info. for entities to be considered [num. ents.="<<nEntTags<<"] : "<<endl;
  cout<<" (example, for an entity three values are needed - 92 2 0)"<<endl;
  cout<<" (where, 92 is entity-tag for specific entity, 2 is entity-dimension for model-face with tag 92 and 0 is closure-flag set to OFF)"<<endl;

  cout<<" "<<endl;
  for(int iEnt=0; iEnt<3*nEntTags; iEnt++)
    cin>>entTagsInfo[iEnt];

  // collect model ents.
  pGEntity *taggedGEnts = new pGEntity[nEntTags];
  for(int iEnt=0; iEnt<nEntTags; iEnt++) {
    taggedGEnts[iEnt] = GM_entityByTag(model,entTagsInfo[3*iEnt+1],entTagsInfo[3*iEnt+0]);
    if(!taggedGEnts[iEnt]) {
      cout<<"\n Error in GetAndSetEntTagsInfo() : cound not find model ent. [type,tags = "<<entTagsInfo[3*iEnt+1]<<","<<entTagsInfo[3*iEnt+0]<<"] \n"<<endl;
      exit(0);
    }
  }

  for(int iEnt=0; iEnt<nEntTags; iEnt++) {
    int isEntTagged;
    pVertex vtx;
    VIter vit = M_classifiedVertexIter(mesh,taggedGEnts[iEnt],entTagsInfo[3*iEnt+2]);
    while(vtx = VIter_next(vit)) {
      if(!EN_getDataInt(vtx,entTaggedId,&isEntTagged))
        EN_attachDataInt(vtx,entTaggedId,1);
    }
    VIter_delete(vit);
  }

  delete [] taggedGEnts;

}

void M2MSolutionTransfer(pGModel source_model, pGModel dest_model, pMesh source_mesh, pMesh dest_mesh, pMeshDataId solutionID, int numVars, int rfOption, int farOption, int halfToFullOption, int entTagsOption) {
  if(rfOption) {
    // memory is allocated for solution on each vertex
    M2MSolutionTransfer1(source_model,dest_model,source_mesh,dest_mesh,solutionID,numVars,farOption,halfToFullOption,entTagsOption);
  }
  else {
    // memory is allocated for solution on each vertex
    M2MSolutionTransfer0(source_model,dest_model,source_mesh,dest_mesh,solutionID,numVars,farOption,halfToFullOption,entTagsOption);
  }
}

// memory is allocated for solution on each vertex
void M2MSolutionTransfer0(pGModel source_model, pGModel dest_model, pMesh source_mesh, pMesh dest_mesh, pMeshDataId solutionID, int numVars, int farOption, int halfToFullOption, int entTagsOption) {

  cout<<endl;
  cout<<" ...looping over each region (works for mixed topology)"<<endl;

  int signHalfToFullOption = 1;
  if(halfToFullOption<0)
    signHalfToFullOption = -1;
  int halfToFullIndex = signHalfToFullOption*halfToFullOption-1;

  int *entTagsInfo;
  pMeshDataId entTaggedId;
  if(entTagsOption) {
    entTagsInfo = new int[3*entTagsOption]; // each unit contains entity-tag, entity-dimension and closure-flag

    entTaggedId = MD_newMeshDataId("tagged mesh ent. in dest");
    GetAndSetEntTagsInfo(dest_model,dest_mesh,entTaggedId,entTagsOption,entTagsInfo);
  }

  pVertex d_vtx;
  VIter d_vit = M_vertexIter(dest_mesh);
  while(d_vtx = VIter_next(d_vit)) {

    if(entTagsOption) {
      int isEntTagged;
      if(!EN_getDataInt((pEntity)d_vtx,entTaggedId,&isEntTagged)) {
	double *solution = new double[numVars];
	for(int iVar=0; iVar<numVars; iVar++)
	  solution[iVar] = 0.; // default solution value set as zero

	double *d_data;
	if(EN_getDataPtr((pEntity)d_vtx,solutionID,(void **)&d_data)) {
	  cout<<"\n Error in M2MSolutionTransfer0() : data attached to dest. vertex [ID="<<EN_id(d_vtx)<<"]\n"<<endl;
	  exit(0);
	}
	EN_attachDataPtr((pEntity)d_vtx,solutionID,(void *)solution);

	continue;
      }
    }

    double orig_d_xyz[3], d_xyz[3];
    V_coord(d_vtx,orig_d_xyz);
    for(int iComp=0; iComp<3; iComp++)
      d_xyz[iComp] = orig_d_xyz[iComp];

    if(halfToFullOption!=0)
      d_xyz[halfToFullIndex] = signHalfToFullOption*fabs(orig_d_xyz[halfToFullIndex]);

    int foundRegion = 0;
    pRegion s_rgn;
    RIter s_rit = M_regionIter(source_mesh);
    while(s_rgn = RIter_next(s_rit)) {
      // int isPointInside = R_isPointInside(s_rgn,d_xyz);
      double xietazeta[3];
      int isPointInside = inverseMap(s_rgn,d_xyz,xietazeta);

      if(isPointInside) {
	double *solution = new double[numVars];
        double *intrSol = InterpolateSolution(s_rgn,solutionID,xietazeta,numVars);
	for(int iVar=0; iVar<numVars; iVar++)
	  solution[iVar] = intrSol[iVar];
        delete [] intrSol;

	double *d_data;
	if(EN_getDataPtr((pEntity)d_vtx,solutionID,(void **)&d_data)) {
	  cout<<"\n Error in M2MSolutionTransfer0() : data attached to dest. vertex [ID="<<EN_id(d_vtx)<<"]\n"<<endl;
	  exit(0);
	}
	EN_attachDataPtr((pEntity)d_vtx,solutionID,(void *)solution);

	foundRegion = 1;
	break;
      }
    }
    RIter_delete(s_rit);

    if(!foundRegion) {
      // could be because of curved geometries

      double minDistSq = 1.e20;
      pVertex s_vtx, s_minVtx = 0;
      VIter s_vit = M_vertexIter(source_mesh);
      while(s_vtx = VIter_next(s_vit)) {
	double s_xyz[3];
	V_coord(s_vtx,s_xyz);

	double distSq = XYZ_distance2(d_xyz,s_xyz);
	if(minDistSq>distSq) {
	  minDistSq = distSq;
	  s_minVtx = s_vtx;
	}
      }
      VIter_delete(s_vit);

      if(!s_minVtx) {
	cout<<"\n Error in M2MSolutionTransfer0() : for vertex ["<<EN_id(d_vtx)<<"] could not find any vertex within distance of "<<sqrt(minDistSq)<<"\n"<<endl;
	exit(0);
      }

      int longEdgeIndex, s_minVtxNumEdges = V_numEdges(s_minVtx);
      double longEdgeLenSq = 0.;
      pEdge s_minVtxEdge;
      for(int iEdge=0; iEdge<s_minVtxNumEdges; iEdge++) {
	s_minVtxEdge = V_edge(s_minVtx,iEdge);

	double s_minVtxEdgeLenSq = E_lengthSq(s_minVtxEdge);
	if(s_minVtxEdgeLenSq>longEdgeLenSq) {
	  longEdgeIndex = iEdge;
	  longEdgeLenSq = s_minVtxEdgeLenSq;
	}
      }

      // double edgeLenSqFactor = 25.;
      if(farOption && minDistSq>edgeLenSqFactor*longEdgeLenSq) {
	cout<<"\n Error in M2MSolutionTransfer0() : for vertex ["<<EN_id(d_vtx)<<"] could not find any vertex within distance of (for source vertex ["<<EN_id(s_minVtx)<<"]: long edge length. = "<<sqrt(longEdgeLenSq)<<") "<<sqrt(edgeLenSqFactor*longEdgeLenSq)<<", distance = "<<sqrt(minDistSq)<<"\n"<<endl;
	exit(0);
      }

      double *s_data;
      if(!EN_getDataPtr((pEntity)s_minVtx,solutionID,(void **)&s_data)) {
	cout<<"\n Error in M2MSolutionTransfer0() : no data attached to nearest source vertex [ID="<<EN_id(s_minVtx)<<"]\n"<<endl;
	exit(0);
      }

      double *solution = new double[numVars];
      for(int iVar=0; iVar<numVars; iVar++)
	solution[iVar] = s_data[iVar];

      double *d_data;
      if(EN_getDataPtr((pEntity)d_vtx,solutionID,(void **)&d_data)) {
	cout<<"\n Error in M2MSolutionTransfer0() : data attached to dest. vertex [ID="<<EN_id(d_vtx)<<"]\n"<<endl;
	exit(0);
      }
      EN_attachDataPtr((pEntity)d_vtx,solutionID,(void *)solution);

    }
  }
  VIter_delete(d_vit);

  if(entTagsOption) {
    delete [] entTagsInfo;

    int isEntTagged;
    VIter d_vit_clean = M_vertexIter(dest_mesh);
    while(d_vtx = VIter_next(d_vit_clean)) {
      if(EN_getDataInt((pEntity)d_vtx,entTaggedId,&isEntTagged))
        EN_deleteData((pEntity)d_vtx,entTaggedId);
    }
    VIter_delete(d_vit_clean);

    MD_deleteMeshDataId(entTaggedId);
  }

}

// memory is allocated for solution on each vertex
void M2MSolutionTransfer1(pGModel source_mode, pGModel dest_model, pMesh source_mesh, pMesh dest_mesh, pMeshDataId solutionID, int numVars, int farOption, int halfToFullOption, int entTagsOption) {

  cout<<endl;
  cout<<" ...using octree to find regions (works only for tets. in source mesh)"<<endl;

  int dataStructureLevel;
  cout<<endl;
  cout<<" Enter data structure level to be used"<<endl;
  cout<<" (1 is max. resolution leading to less search time)"<<endl;
  cout<<" (typically between 1-5)"<<endl;
  cout<<" ";
  cin>>dataStructureLevel;

  pProgress progress;
  pRegion s_region;
  RIter s_rit = M_regionIter(source_mesh);
  while(s_region=RIter_next(s_rit)) {
    int s_regionTopoType = R_topoType(s_region);

    if(s_regionTopoType==Rtet)
      continue;

    cout<<"\n Error in M2MSolutionTransfer1() : region topology ["<<s_regionTopoType<<"] NOT supported in source mesh\n"<<endl;
    exit(0);
  }
  RIter_delete(s_rit);

  pMeshRegionFinder mrf = MeshRegionFinder_new(source_mesh,dataStructureLevel,0, progress);

  int signHalfToFullOption = 1;
  if(halfToFullOption<0)
    signHalfToFullOption = -1;
  int halfToFullIndex = signHalfToFullOption*halfToFullOption-1;

  int *entTagsInfo;
  pMeshDataId entTaggedId;
  if(entTagsOption) {
    entTagsInfo = new int[3*entTagsOption]; // each unit contains entity-tag, entity-dimension and closure-flag

    entTaggedId = MD_newMeshDataId("tagged mesh ent. in dest");
    GetAndSetEntTagsInfo(dest_model,dest_mesh,entTaggedId,entTagsOption,entTagsInfo);
  }

  pVertex d_vtx;
  VIter d_vit = M_vertexIter(dest_mesh);
  while(d_vtx = VIter_next(d_vit)) {

    if(entTagsOption) {
      int isEntTagged;
      if(!EN_getDataInt((pEntity)d_vtx,entTaggedId,&isEntTagged)) {
	double *solution = new double[numVars];
	for(int iVar=0; iVar<numVars; iVar++)
	  solution[iVar] = 0.; // default solution value set as zero

	double *d_data;
	if(EN_getDataPtr((pEntity)d_vtx,solutionID,(void **)&d_data)) {
	  cout<<"\n Error in M2MSolutionTransfer0() : data attached to dest. vertex [ID="<<EN_id(d_vtx)<<"]\n"<<endl;
	  exit(0);
	}
	EN_attachDataPtr((pEntity)d_vtx,solutionID,(void *)solution);

	continue;
      }
    }

    double orig_d_xyz[3], d_xyz[3];
    V_coord(d_vtx,orig_d_xyz);
    for(int iComp=0; iComp<3; iComp++)
      d_xyz[iComp] = orig_d_xyz[iComp];

    if(halfToFullOption!=0)
      d_xyz[halfToFullIndex] = signHalfToFullOption*fabs(orig_d_xyz[halfToFullIndex]);

    int foundRegion = 0;
    double distance;
    pRegion s_rgn = MeshRegionFinder_find(mrf,d_xyz,0,&distance);

    if(s_rgn) {
      double *solution;

      if(distance*distance>distSqTol) {
	pPList s_rgnEdges = R_edges(s_rgn,1);
	int longEdgeIndex, s_rgnNumEdges = PList_size(s_rgnEdges);
	double longEdgeLenSq = 0.;
	pEdge s_rgnEdge;
	for(int iEdge=0; iEdge<s_rgnNumEdges; iEdge++) {
	  s_rgnEdge = (pEdge)PList_item(s_rgnEdges,iEdge);

	  double s_rgnEdgeLenSq = E_lengthSq(s_rgnEdge);
	  if(s_rgnEdgeLenSq>longEdgeLenSq) {
	    longEdgeIndex = iEdge;
	    longEdgeLenSq = s_rgnEdgeLenSq;
	  }
	}
	PList_delete(s_rgnEdges);

	// double edgeLenSqFactor = 25.;
	if(farOption && distance*distance>edgeLenSqFactor*longEdgeLenSq) {
	  cout<<"\n Error in M2MSolutionTransfer1() : for vertex ["<<EN_id(d_vtx)<<"] could not find any vertex within distance of (for source region ["<<EN_id(s_rgn)<<"]: long edge length. = "<<sqrt(longEdgeLenSq)<<") "<<sqrt(edgeLenSqFactor*longEdgeLenSq)<<", distance = "<<distance<<"\n"<<endl;
	  exit(0);
	}

        solution = new double[numVars];
        for(int iVar=0; iVar<numVars; iVar++)
  	  solution[iVar] = 0.;

        pPList s_rverts = R_vertices(s_rgn,1);
        int s_numRVerts = PList_size(s_rverts);
        for(int iRVert=0; iRVert<s_numRVerts; iRVert++) {
	  double *s_data;
	  pVertex s_rvtx = (pVertex)PList_item(s_rverts,iRVert);
	  if(!EN_getDataPtr((pEntity)s_rvtx,solutionID,(void **)&s_data)) {
	    cout<<"\n Error in M2MSolutionTransfer1() : no data attached to source region's vertex [ID="<<EN_id(s_rvtx)<<"]\n"<<endl;
	    exit(0);
	  }

	  for(int iVar=0; iVar<numVars; iVar++)
  	    solution[iVar] += s_data[iVar];

        }
        PList_delete(s_rverts);

        for(int iVar=0; iVar<numVars; iVar++)
	  solution[iVar] /= s_numRVerts;
      }
      else {
        double xietazeta[3];
        if(!inverseMap(s_rgn,d_xyz,xietazeta)) {
	  cout<<"\n Error in M2MSolutionTransfer1() : point inside region could not be parameterized [pnt. xyz="<<orig_d_xyz[0]<<", "<<orig_d_xyz[1]<<", "<<orig_d_xyz[2]<<", rgn. ID="<<EN_id(s_rgn)<<"]\n"<<endl;
  	  exit(0);
        }

	solution = new double[numVars];
        double *intrSol = InterpolateSolution(s_rgn,solutionID,xietazeta,numVars);
	for(int iVar=0; iVar<numVars; iVar++)
	  solution[iVar] = intrSol[iVar];
        delete [] intrSol;
      }

      double *d_data;
      if(EN_getDataPtr((pEntity)d_vtx,solutionID,(void **)&d_data)) {
	cout<<"\n Error in M2MSolutionTransfer1() : data attached to dest. vertex [ID="<<EN_id(d_vtx)<<"]\n"<<endl;
	exit(0);
      }
      EN_attachDataPtr((pEntity)d_vtx,solutionID,(void *)solution);

      foundRegion = 1;
    }

    if(!foundRegion) {
      // could be because of curved geometries
      // cout<<"\n WARNING : for vertex ["<<EN_id(d_vtx)<<"] could not find source region (using nearest source vertex for solution)\n"<<endl;

      double minDistSq = 1.e20;
      pVertex s_vtx, s_minVtx = 0;
      VIter s_vit = M_vertexIter(source_mesh);
      while(s_vtx = VIter_next(s_vit)) {
	double s_xyz[3];
	V_coord(s_vtx,s_xyz);

	double distSq = XYZ_distance2(d_xyz,s_xyz);
	if(minDistSq>distSq) {
	  minDistSq = distSq;
	  s_minVtx = s_vtx;
	}
      }
      VIter_delete(s_vit);

      if(!s_minVtx) {
	cout<<"\n Error in M2MSolutionTransfer1() : for vertex ["<<EN_id(d_vtx)<<"] could not find any vertex within distance of "<<sqrt(minDistSq)<<"\n"<<endl;
	exit(0);
      }

      int longEdgeIndex, s_minVtxNumEdges = V_numEdges(s_minVtx);
      double longEdgeLenSq = 0.;
      pEdge s_minVtxEdge;
      for(int iEdge=0; iEdge<s_minVtxNumEdges; iEdge++) {
	s_minVtxEdge = V_edge(s_minVtx,iEdge);

	double s_minVtxEdgeLenSq = E_lengthSq(s_minVtxEdge);
	if(s_minVtxEdgeLenSq>longEdgeLenSq) {
	  longEdgeIndex = iEdge;
	  longEdgeLenSq = s_minVtxEdgeLenSq;
	}
      }

      // double edgeLenSqFactor = 25.;
      if(farOption && minDistSq>edgeLenSqFactor*longEdgeLenSq) {
	cout<<"\n Error in M2MSolutionTransfer1() : for vertex ["<<EN_id(d_vtx)<<"] could not find any vertex within distance of (for source vertex ["<<EN_id(s_minVtx)<<"]: long edge length. = "<<sqrt(longEdgeLenSq)<<") "<<sqrt(edgeLenSqFactor*longEdgeLenSq)<<", distance = "<<sqrt(minDistSq)<<"\n"<<endl;
	exit(0);
      }

      double *s_data;
      if(!EN_getDataPtr((pEntity)s_minVtx,solutionID,(void **)&s_data)) {
	cout<<"\n Error in M2MSolutionTransfer1() : no data attached to nearest source vertex [ID="<<EN_id(s_minVtx)<<"]\n"<<endl;
	exit(0);
      }

      double *solution = new double[numVars];
      for(int iVar=0; iVar<numVars; iVar++)
	solution[iVar] = s_data[iVar];

      double *d_data;
      if(EN_getDataPtr((pEntity)d_vtx,solutionID,(void **)&d_data)) {
	cout<<"\n Error in M2MSolutionTransfer1() : data attached to dest. vertex [ID="<<EN_id(d_vtx)<<"]\n"<<endl;
	exit(0);
      }
      EN_attachDataPtr((pEntity)d_vtx,solutionID,(void *)solution);

    }
  }
  VIter_delete(d_vit);

  if(entTagsOption) {
    delete [] entTagsInfo;

    int isEntTagged;
    VIter d_vit_clean = M_vertexIter(dest_mesh);
    while(d_vtx = VIter_next(d_vit_clean)) {
      if(EN_getDataInt((pEntity)d_vtx,entTaggedId,&isEntTagged))
        EN_deleteData((pEntity)d_vtx,entTaggedId);
    }
    VIter_delete(d_vit_clean);

    MD_deleteMeshDataId(entTaggedId);
  }

  MeshRegionFinder_delete(mrf);

}
