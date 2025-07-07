#include "MeshSim.h"

#include <math.h>
#include <iostream>
#include <stdlib.h>

using std::cout;
using std::endl;

#ifdef __cplusplus
extern "C" {
#endif

void ludcmp_( double*, int*, int*, int*, double* );
void lubksb_( double*, int*, int*, int*, double* );

void
display_region( pRegion region ) {

  double xyz[3];
  pPList vertices = R_vertices( region, 1 );
  cout << "-------------------"<< endl ;
  for( int i=0; i<PList_size( vertices ); i++){
      V_coord( (pVertex)PList_item( vertices, i ), xyz );
      cout << xyz[0] <<" "<< xyz[1]<<" "<<xyz[2]<< endl;
  }
  cout << "-------------------"<< endl ;
}

int 
inverseMap( pRegion region, 
            double* qpt,
            double* pt ) {

    // This is the version of Inverse map, which uses the algorithm 
    // by ken from Mesh2Mesh  (MTMHO3)
    // This thing basically, does a linear solution and then tries to
    // get a Iterative Newton correction to it.

    // First to setup the constants of the forward transformation
    //                   x = Ax* xi
    //                   y = Ay* xi
    //                   z = Az* xi
    // Ax,Ay,Az have 8 terms each and can be obtained using the 
    // solution of an 8x8 system,( which is what I am going to do)!!

    double** A;
			   
    double tet[4][4] = { {1., 0., 0., 0.},
                         {1., 1., 0., 0.},
                	 {1., 0., 1., 0.},
                         {1., 0., 0., 1.} };

    double pyramid[5][5] = { {1., -1., -1., -1.,  1.},
                             {1.,  1., -1., -1., -1.},
                             {1.,  1.,  1., -1.,  1.},
                             {1., -1.,  1., -1., -1.},
                             {1.,  0.,  0.,  1.,  0.} };

    double prism[6][6] = { {1., 0., 0., -1.,  0.,  0.},
                           {1., 1., 0., -1., -1.,  0.},
                           {1., 0., 1., -1.,  0., -1.},
                           {1., 0., 0.,  1.,  0.,  0.},
                           {1., 1., 0.,  1.,  1.,  0.},
	       	           {1., 0., 1.,  1.,  0.,  1.} };

    double x = qpt[0], y = qpt[1], z = qpt[2];
    double *xel, *yel, *zel;
    double xisol[3];
    double xyz[3];
    double fnumber;

    pPList verts = R_vertices( region ,0);
    int nshl = PList_size(verts);
    int *indx = new int[nshl];
    xel = new double[nshl]; yel = new double[nshl]; zel = new double[nshl];

    A = new double* [3];
    for(int i=0; i<3; i++) A[i] = new double[nshl];

    double *M = new double[nshl*nshl];

    int rTopoType = R_topoType(region);
    int k=0;
    switch(rTopoType) {
    case Rtet :
      //creating the LHS
      for(int i=0; i<nshl; i++)
        for(int j=0; j<nshl; j++)
	  M[k++]= tet[i][j];
      break;
    case Rpyramid :
      //creating the LHS
      for(int i=0; i<nshl; i++)
        for(int j=0; j<nshl; j++)
	  M[k++]= pyramid[i][j];
      break;
    case Rwedge :
      //creating the LHS
      for(int i=0; i<nshl; i++)
        for(int j=0; j<nshl; j++)
	  M[k++]= prism[i][j];
      break;
    default:
      cout<<"\nError in inverseMap()..."<<endl;
      cout<<"element topology ["<<rTopoType<<"] NOT supported\n"<<endl;
      exit(0);
    }
    
    // LU decompsing the coeff matrix
    ludcmp_( M, &nshl, &nshl, indx, &fnumber);
  
    // Creating the RHS
    for(int i=0; i<nshl; i++){
      V_coord( ( pVertex ) PList_item( verts, i ), xyz );
      xel[i] = xyz[0];
      yel[i] = xyz[1];
      zel[i] = xyz[2];
    }
    PList_delete(verts);

    for(int i=0; i<nshl; i++){
      A[0][i] = xel[i];
      A[1][i] = yel[i];
      A[2][i] = zel[i];
    }

    // Now back substituting to get back the correct set of constants
    // for this element.
  
    lubksb_( M, &nshl, &nshl, indx, A[0] );
    lubksb_( M, &nshl, &nshl, indx, A[1] );
    lubksb_( M, &nshl, &nshl, indx, A[2] );

    // Now we have Ax, Ay and Az (where A is the inverse of matrix Mtemp). 
    // Next, we try to get xi, zeta, eta for a given x, y, z.
    // Ax contains the alpha_x in the form of A[0][0] = alpha_x0, 
    // A[0][1] = alpha_x1, and so on. A[1][0] = alpha_y0, 
    // A[1][1] = alpha_y1, and so on. A[2][0] = alpha_z0, 
    // A[2][1] = alpha_z1, and so on. 

    int indx2[3];
    int three=3;

    double MS[9];
    k = 0;
    for(int i=0; i<3; i++){
      for(int j=1; j<4; j++){
	MS[k++] = A[i][j];
      }
    }

    double xl[3], dxl[3];
    xisol[0] = x - A[0][0];
    xisol[1] = y - A[1][0];
    xisol[2] = z - A[2][0];
    xl[0] = x;
    xl[1] = y;
    xl[2] = z;

    // LU decompsing the coeff matrix and solving for xisol
    ludcmp_( MS, &three, &three, indx2, &fnumber);
    lubksb_( MS, &three, &three, indx2, xisol);

    double tol = 0.000001;
    int flag = 1, loop = 0;

    while(flag) {
      switch(rTopoType) {
      case Rtet :
	k = 0;
	for(int i=0; i<3; i++){
	  dxl[i] = xl[i] - (A[i][0] + A[i][1]*xisol[0] + A[i][2]*xisol[1] + A[i][3]*xisol[2]);
	  MS[k++] = A[i][1];
	  MS[k++] = A[i][2];
	  MS[k++] = A[i][3];
	}
	break;
      case Rpyramid :
	k = 0;
	for(int i=0; i<3; i++){
	  dxl[i] = xl[i] - (A[i][0] + A[i][1]*xisol[0] + A[i][2]*xisol[1] + 
			    A[i][3]*xisol[2] +  A[i][4]*xisol[0]*xisol[1]);
	  MS[k++] = A[i][1] + A[i][4]*xisol[1];
	  MS[k++] = A[i][2] + A[i][4]*xisol[0];
	  MS[k++] = A[i][3];
	}
	break;
      case Rwedge :
	k = 0;
	for(int i=0; i<3; i++){
	  dxl[i] = xl[i] - (A[i][0] + A[i][1]*xisol[0] + A[i][2]*xisol[1] + 
			    A[i][3]*xisol[2] +  A[i][4]*xisol[0]*xisol[2] + A[i][5]*xisol[1]*xisol[2]);
	  MS[k++] = A[i][1] + A[i][4]*xisol[2];
	  MS[k++] = A[i][2] + A[i][5]*xisol[2];
	  MS[k++] = A[i][3] + A[i][4]*xisol[0] + A[i][5]*xisol[1];
	}
	break;
    default:
      cout<<"\nError in inverseMap()..."<<endl;
      cout<<"element topology ["<<rTopoType<<"] NOT supported\n"<<endl;
      exit(0);
      }

      // LU decompsing the coeff matrix and solving for xisol
      ludcmp_( MS, &three, &three, indx2, &fnumber);
      lubksb_( MS, &three, &three, indx2, dxl);

      for(int i=0; i<3; i++) xisol[i]=xisol[i]+dxl[i];
      if( (fabs(dxl[0]) < tol) && (fabs(dxl[1]) < tol) && (fabs(dxl[2]) < tol) )
	flag = 0;
      loop ++;
      
      //maximum recursive correction is 10 steps for linear
      if(loop == 10)
	flag = 0;
    }

    int truth = 1;                        
    if(loop == 10)
      truth = 0;
  
    if (truth) {
      pt[0] = xisol[0];                 
      pt[1] = xisol[1];                 
      pt[2] = xisol[2];  
    }                                    

    delete [] indx;
    delete [] xel;
    delete [] yel;
    delete [] zel;

    for(int i=0; i<3; i++) delete [] A[i] ;
    delete [] A;

    delete [] M;

    return truth;
}

#ifdef __cplusplus
}
#endif
