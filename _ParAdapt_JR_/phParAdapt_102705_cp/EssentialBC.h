#ifndef H_EssentialBC
#define H_EssentialBC

#include <iostream.h>
#include "MeshSimInternal.h"      // for gType in essentialBC

/* compatibility */
#include "MeshSim.h"
#include "ModelTypes.h"
#include "BoundaryCondition.h"


class EssentialBC : public BoundaryCondition {
public:
  EssentialBC (pGFace gface);
  EssentialBC (pGEdge gedge);
  EssentialBC (pGVertex gvert);
  int eval (pVertex, double *);
  int eval (pEdge, double *, int);      // for hierarchic basis
  int eval (pFace, double *, int);      // for hierarchic basis
  void takeBCfromIC(double *BC, double *qTot, int nshgTot, int GDOFnum);
private:
  gType gtype;                            // Gvertex, Gedge or Gface
  pGFace gf;     // could store in a union since either face, edge or vertex
  pGEdge ge;     //   but too much bother for insignificant memory savings
  pGVertex gv;
  int zxcl;      // Flag to do Axisymmetric Centerline

  void update_inherit ();     // to update dontinherit
  void update_thermo(const double*, double*, int*);  // thermo BC's
  void update_velo(const double*, double*, int*);    // velo BC's
  void update_scalar(const double*, double*, int*);  // scalar BC's
  void update_axisym(const double*, double*, int*);  // axisym centerline
  void pteval (const double*, double *, int *);
  int isZA()
  {
    GEntity* gty;
    switch (gtype){
    case Gface:
      gty = (pGEntity) gf;
      break;
    case Gedge:
      gty = (pGEntity) ge;
      break;
    case Gvertex:
      gty = (pGEntity) gv;
      break;
    default:
      cerr << " Entity does not have type "<<endl;
      exit(1);
    }
    return GEN_dataI(gty,"AXCL");
  }
};

#endif
