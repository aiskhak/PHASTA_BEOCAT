#include "EN_isPointInside.h"

#include "MeshSim.h"
#include "MeshTypes.h"

#include <iostream>
#include <cmath>
#include <stdlib.h>

using std::cout;
using std::endl;

double M_getTolerance() {
  return 1.e-14;
}

int EN_isPointInside(pEntity entity, double *xyz) {

  int entType = EN_type(entity);
  switch(entType) {
  case Tface :
    return F_isPointInside((pFace)entity,xyz);
  case Tregion :
    return R_isPointInside((pRegion)entity,xyz);
  case Tvertex :
  case Tedge :
  default:
    cout<<"\nError in EN_isPointInside()..."<<endl;
    cout<<"entType ["<<entType<<"] NOT supported as of now\n"<<endl;
    exit(0);
  }

}

// see combination : <return value> - <where point lies >
// 1 - purely inside face
// 2-v0, 3-v1, 4-v2
// 5-e0, 6-e1, 7-e2
int F_isPointInside(pFace face, double *xyz) {
  if(F_numEdges(face)!=3) {
    cout<<"\nError in F_isPointInside()..."<<endl;
    cout<<"face is NOT a triangle"<<endl;
    exit(0);
  }

//   double mtol = M_getTolerance();
//   double area, sumAreas = 0., subArea;

  double fxyz[3][3];
  F_coord(face,fxyz);

  return XYZ_isPointInsideTriangle(fxyz,xyz);

//   for(int i=0; i<3; i++) {
//     if(XYZ_distance2(fxyz[i],xyz)<mtol*mtol)
//       return i+2;
//   }

//   double v[3][3];
//   diffVt(fxyz[1],fxyz[0],v[0]);
//   diffVt(fxyz[2],fxyz[1],v[1]);
//   diffVt(fxyz[0],fxyz[2],v[2]);

//   double normal[3];
//   crossProd(v[0],v[1],normal);
//   area = 0.5*sqrt(dotProd(normal,normal));
//   double mtolArea = 1.e-8*area;

//   int counter = 0;
//   double temp[3], nor[3];
//   for(int i=0; i<3; i++) {
//     diffVt(xyz,fxyz[i],temp);
//     crossProd(v[i],temp,nor);

//     subArea = 0.5*sqrt(dotProd(nor,nor));
//     // subArea is always positive
//     if(subArea<mtolArea) {
//       counter = i+1;
//       continue;
//     }

//     // to avoid trouble due to tolerance value
//     // (for example, subArea=2.e-14 and mtol=1.e-14)
//     // hence, using mtolArea, i.e., percentage of area 
//     // might be helpful than mtol, i.e., M_getTolerance()
//     if(dotProd(nor,normal)<0.)
//       return 0;
    
//     sumAreas += subArea;
//   }

//   if(fabs(area-sumAreas)<mtolArea) {
//     if(counter)
//       return 5+counter-1;
      

//     return 1;
//   }

//   return 0;
}

// see combination : <return value> - <where point lies >
// 1 - purely inside face
// 2-v0, 3-v1, 4-v2
// 5-e0, 6-e1, 7-e2
int XYZ_isPointInsideTriangle(double tri_xyz[3][3], double *xyz) {

  double mtol = M_getTolerance();
  double area, sumAreas = 0., subArea;

  for(int i=0; i<3; i++) {
    if(XYZ_distance2(tri_xyz[i],xyz)<mtol*mtol)
      return i+2;
  }

  double v[3][3];
  diffVt(tri_xyz[1],tri_xyz[0],v[0]);
  diffVt(tri_xyz[2],tri_xyz[1],v[1]);
  diffVt(tri_xyz[0],tri_xyz[2],v[2]);

  double normal[3];
  crossProd(v[0],v[1],normal);
  area = 0.5*sqrt(dotProd(normal,normal));
  double mtolArea = 1.e-8*area;

  int counter = 0;
  double temp[3], nor[3];
  for(int i=0; i<3; i++) {
    diffVt(xyz,tri_xyz[i],temp);
    crossProd(v[i],temp,nor);

    subArea = 0.5*sqrt(dotProd(nor,nor));
    // subArea is always positive
    if(subArea<mtolArea) {
      counter = i+1;
      continue;
    }

    // to avoid trouble due to tolerance value
    // (for example, subArea=2.e-14 and mtol=1.e-14)
    // hence, using mtolArea, i.e., percentage of area 
    // might be helpful than mtol, i.e., M_getTolerance()
    if(dotProd(nor,normal)<0.)
      return 0;
    
    sumAreas += subArea;
  }

  if(fabs(area-sumAreas)<mtolArea) {
    if(counter)
      return 5+counter-1;

    return 1;
  }

  return 0;
}

double XYZ_volume2(dArray *X_Y_Z, int rTopoType) {

  int numChildTets;
  int elemType;
  // fit two tets. in a pyramid
  // (diagonal edge of a pyramid is : V0-V2)
  // fit three tets. in a prism/wedge
  // (diagonal edges of a prism are : V0-V4, V1-V5 and V5-V0)
  // int tetVerts[3][4] = {{0,1,2,5}, {0,1,5,4}, {0,4,5,3}};
  int mapNodes[3][3][4] ={{{0,1,2,3}, {-1,-1,-1,-1}, {-1,-1,-1,-1}},
                          {{0,1,2,4}, {0,2,3,4}, {-1,-1,-1,-1}},
			  {{0,1,2,5}, {0,1,5,4}, {0,4,5,3}}};

  switch(rTopoType) {
  case Rtet :
    numChildTets = 1;
    elemType = 0;
    break;
  case Rpyramid :
    numChildTets = 2;
    elemType = 1;
    break;
  case Rwedge :
    numChildTets = 3;
    elemType = 2;
    break;
  default:
    cout<<"\nError in XYZ_volume2()..."<<endl;
    cout<<"region topology ["<<rTopoType<<"] NOT supported"<<endl;
    exit(0);
  }

  double volume = 0.;
  dArray sub_xyz[4];
  for(int iTet=0; iTet<numChildTets; iTet++) {
    for(int iTVert=0; iTVert<4; iTVert++)
      for(int iComp=0; iComp<3; iComp++)
	sub_xyz[iTVert][iComp] = X_Y_Z[mapNodes[elemType][iTet][iTVert]][iComp];

    volume += XYZ_volume(sub_xyz);
  }

  return volume;

}

double R_volume2(pRegion region) {

  int topoType = R_topoType(region);

  int numChildTets;
  int elemType;
  // fit two tets. in a pyramid
  // (diagonal edge of a pyramid is : V0-V2)
  // fit three tets. in a prism/wedge
  // (diagonal edges of a prism are : V0-V4, V1-V5 and V5-V0)
  // int tetVerts[3][4] = {{0,1,2,5}, {0,1,5,4}, {0,4,5,3}};
  int mapNodes[3][3][4] ={{{0,1,2,3}, {-1,-1,-1,-1}, {-1,-1,-1,-1}},
                          {{0,1,2,4}, {0,2,3,4}, {-1,-1,-1,-1}},
			  {{0,1,2,5}, {0,1,5,4}, {0,4,5,3}}};

  switch(topoType) {
  case Rtet :
    numChildTets = 1;
    elemType = 0;
    break;
  case Rpyramid :
    numChildTets = 2;
    elemType = 1;
    break;
  case Rwedge :
    numChildTets = 3;
    elemType = 2;
    break;
  default:
    cout<<"\nError in R_volume2()..."<<endl;
    cout<<"region topology ["<<topoType<<"] NOT supported"<<endl;
    exit(0);
  }

  dArray *X_Y_Z;
  pPList rverts = R_vertices(region,1);
  int numRVerts = PList_size(rverts);
  X_Y_Z = new dArray[numRVerts];
  pVertex vtx;
  for(int iRVert=0; iRVert<numRVerts; iRVert++) {
    vtx = (pVertex)PList_item(rverts,iRVert);
    V_coord(vtx,X_Y_Z[iRVert]);
  }
  PList_delete(rverts);

  double volume = 0.;
  dArray sub_xyz[4];
  for(int iTet=0; iTet<numChildTets; iTet++) {
    for(int iTVert=0; iTVert<4; iTVert++)
      for(int iComp=0; iComp<3; iComp++)
	sub_xyz[iTVert][iComp] = X_Y_Z[mapNodes[elemType][iTet][iTVert]][iComp];

    volume += XYZ_volume(sub_xyz);
  }

  delete [] X_Y_Z;

  return volume;
}

// if point lies inside return :
// 1 - purely inside region
// 2 - purely inside one of the face (for wedge could be diagonal face)
// 3 - purely inside one of the edge (for wedge could be diagonal edge)
// 4 - on one of the vertex
int R_isPointInside(pRegion region, double *xyz) {

  int topoType = R_topoType(region);

  if(!(topoType==Rtet || topoType==Rpyramid || topoType==Rwedge)) {
    cout<<"\nError in R_isPointInside()..."<<endl;
    cout<<"region is NOT a tet. or pyramid or wedge"<<endl;
    exit(0);
  }

  int value = 0;  
  double X_Y_Z[4][3];
  pVertex vtx;
  if(topoType==Rtet) {
    pPList rverts = R_vertices(region,1);
    for(int iVtx=0; iVtx<4; iVtx++) {
      vtx = (pVertex)PList_item(rverts,iVtx);
      V_coord(vtx,X_Y_Z[iVtx]);
    }
    PList_delete(rverts);
    value = XYZ_isPointInside(X_Y_Z,xyz);
  }
  else if(topoType==Rpyramid) {
    // fit two tets in a pyramid
    // (diagonal edge of a pyramid is : V0-V2)
    int tetVerts[2][4] = {{0,1,2,4}, {0,2,3,4}};

    pPList rverts = R_vertices(region,1);
    for(int iTet=0; iTet<2; iTet++) {
      for(int iVtx=0; iVtx<4; iVtx++) {
	vtx = (pVertex)PList_item(rverts,tetVerts[iTet][iVtx]);
	V_coord(vtx,X_Y_Z[iVtx]);
      }
      value = XYZ_isPointInside(X_Y_Z,xyz); 
      if(value)
	break;
    }
    PList_delete(rverts);
  }
  else if(topoType==Rwedge) {
    // fit three tets in a prism/wedge
    // (diagonal edges of a prism are : V0-V4, V1-V5 and V5-V0)
    int tetVerts[3][4] = {{0,1,2,5}, {0,1,5,4}, {0,4,5,3}};

    pPList rverts = R_vertices(region,1);
    for(int iTet=0; iTet<3; iTet++) {
      for(int iVtx=0; iVtx<4; iVtx++) {
	vtx = (pVertex)PList_item(rverts,tetVerts[iTet][iVtx]);
	V_coord(vtx,X_Y_Z[iVtx]);
      }
      value = XYZ_isPointInside(X_Y_Z,xyz); 
      if(value)
	break;
    }
    PList_delete(rverts);
  }

  return value;

//   int counter = 1;
//   double mtol = M_getTolerance();  
//   double volume, sumVols = 0., subVolume;
//   volume = R_volume(region);
//   double mtolVol = 1.e-8*volume;

//   double sub_xyz[4][3];
//   for(int iComp=0; iComp<3; iComp++)
//     sub_xyz[3][iComp] = xyz[iComp];
//   pFace face;
//   pPList verts;  
//   for(int iFace=0; iFace<4; iFace++) {
//     face = R_face(region,iFace);

//     verts = F_vertices(face,1-R_dirUsingFace(region,face));    
//     for(int iVert=0; iVert<3; iVert++)
//       V_coord((pVertex)PList_item(verts,iVert),sub_xyz[iVert]);
//     PList_delete(verts);

//     subVolume = XYZ_volume(sub_xyz);

//     // to avoid trouble due to tolerance value
//     // (for example, subVolume=2.e-14 and mtol=1.e-14)
//     // hence, using mtolVol, i.e., percentage of Vol
//     // might be helpful than mtol, i.e., M_getTolerance()
//     if(fabs(subVolume)<mtolVol) {
//       switch(F_isPointInside(face,xyz)) {
//       case 0 :
// 	return 0;
//       case 1 : // purely inside face
// 	return 2;
//       case 2 : // purely on vertex
//       case 3 : // purely on vertex
//       case 4 : // purely on vertex
// 	return 4;
//       case 5 : // purely inside edge
//       case 6 : // purely insde edge
//       case 7 : // purely inside edge
// 	return 3;
//       }
//     }
//     else if(subVolume<0.0)
//       return 0;
//     sumVols += subVolume;
//   }
 
//   if(fabs(volume-sumVols)<mtolVol)
//     return counter;

//   return 0;
}

// if point lies inside return :
// 1 - purely inside region
// 2 - purely inside one of the face
// 3 - purely inside one of the edge
// 4 - on one of the vertex
// this routine assumes four verts, i.e., tet. element
// (other element topologies can be broken/sub-divided into child tets.)
int XYZ_isPointInside(dArray *X_Y_Z, double *xyz) {
  int value = 1;
  // double mtol = M_getTolerance();
  double volume, sumVols = 0., subVolume;
  volume = XYZ_volume(X_Y_Z);

  // to avoid trouble due to tolerance value
  // (for example, subVolume=2.e-14 and mtol=1.e-14)
  // hence, using mtolVol, i.e., percentage of Vol
  // might be helpful than mtol, i.e., M_getTolerance()
  double mtolVol = 1.e-8*volume;

  double sub_xyz[4][3];
  for(int iComp=0; iComp<3; iComp++)
    sub_xyz[3][iComp] = xyz[iComp];

  int mapChildNode[4][3] = {{0,1,2}, {3,1,0}, {3,2,1}, {3,0,2}};

  for(int iChild=0; iChild<4; iChild++) {
    // 4th node is the point to be checked
    for(int iNode=0; iNode<3; iNode++)
      for(int iDir=0; iDir<3; iDir++)
	sub_xyz[iNode][iDir] = X_Y_Z[mapChildNode[iChild][iNode]][iDir];

    subVolume = XYZ_volume(sub_xyz);

    if(fabs(subVolume)<mtolVol) {
      // value++;
      switch(XYZ_isPointInsideTriangle(sub_xyz,xyz)) {
      case 0 :
   	return 0;
      case 1 : // purely inside face
	return 2;
      case 2 : // purely on vertex
      case 3 : // purely on vertex
      case 4 : // purely on vertex
	return 4;
      case 5 : // purely inside edge
      case 6 : // purely insde edge
      case 7 : // purely inside edge
	return 3;
      }
    }

    if(subVolume<0.0)
      return 0;

    sumVols += subVolume;
  }

  if(fabs(volume-sumVols)<mtolVol)
    return value;                              

  return 0;
}
