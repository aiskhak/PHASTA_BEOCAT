#include <stdio.h>
//#include <fstream>
#include "MeshSimInternal.h"
#include "SimPMesh.h"
#include "MeshSim.h"
#include "MeshSimAdapt.h"
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <iostream.h>

#include "phParAdapt.h"

void
display_region( pRegion region ){

  double xyz[3];
  pPList vertices = R_vertices( region, 1 );
  cout << "-------------------"<< endl ;
  for( int i=0; i<PList_size( vertices ); i++) {
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
    double** AA;
    static double M[4][4] ={ {1, 1, 0, 0 }, 
                             {1, 0, 1, 0 },
                             {1, 0, 0, 1 },
                             {1, 0, 0, 0 } };
			   
			   
			   
			   
    int eight = 8;
    int four  = 4;
    static double Mtemp [16];
    double x = qpt[0];
    double y = qpt[1];
    double z = qpt[2];
    double xel[8],yel[8],zel[8];
    double xisol[3];
    double xyz[3];
    int indx[4];
    double fnumber;
    pPList verts = R_vertices( region ,0);

    A = new double* [3];
    AA = new double*[3];
    for(int i =0; i< 3;i++) A[i] = new double [4];
    for(int i =0; i< 3;i++) AA[i] = new double [4];  


    //creating the LHS
    int k=0;
    for(int i =0; i<4; i++)
        for(int j=0; j<4; j++)
            Mtemp[k++]= M[i][j]; 
  
    
  // LU decompsing the coeff matrix
    ludcmp_( Mtemp, &four, &four, indx, &fnumber);
  
    // Creating the RHS
    for(int i=0; i< 4; i++){
        V_coord( ( pVertex ) PList_item( verts, i ), xyz );
        xel[i] = xyz[0];
        yel[i] = xyz[1];
        zel[i] = xyz[2];
    }

    PList_delete(verts);

    for(int i=0; i<4;i++){
        A[0][i]=xel[i];
        A[1][i]=yel[i];
        A[2][i]=zel[i];
    }

    // Now back substituting to get back the correct set of constants
    // for this element.
  
    lubksb_( Mtemp, &four, &four, indx, A[0] );
    lubksb_( Mtemp, &four, &four, indx, A[1] );
    lubksb_( Mtemp, &four, &four, indx, A[2] );


  // Now we have Ax, Ay and Az (where A is the inverse of matrix Mtemp). 
  // Next, we try to get xi, zeta, eta for a given x, y, z.
  // Ax contains the alpha_x in the form of A[0][0] = alpha_x0, 
  // A[0][1] = alpha_x1, and so on. A[1][0] = alpha_y0, 
  // A[1][1] = alpha_y1, and so on. A[2][0] = alpha_z0, 
  // A[2][1] = alpha_z1, and so on. 

  // But first, overwrite the alphas with the solution solved for by
  // paper and pencil.

    AA[0][0] = xel[3];
    AA[0][1] = ( xel[0] - xel[3] );
    AA[0][2] = ( xel[1] - xel[3] );
    AA[0][3] = ( xel[2] - xel[3] );

    AA[1][0] = yel[3];
    AA[1][1] = ( yel[0] - yel[3] );
    AA[1][2] = ( yel[1] - yel[3] );
    AA[1][3] = ( yel[2] - yel[3] );

    AA[2][0] = zel[3];
    AA[2][1] = ( zel[0] - zel[3] );
    AA[2][2] = ( zel[1] - zel[3] );
    AA[2][3] = ( zel[2] - zel[3] );


    int indx2[3];
    int three=3;

    double MS[9];
    k =0;
    for(int i =0;i<3;i++){
        for(int j=1; j<4;j++){
            MS[k++]= A[i][j];
        }
    }

    double xl[3],dxl[3];
    xisol[0] =  x - A[0][0];
    xisol[1] =  y - A[1][0];
    xisol[2] =  z - A[2][0];
    xl[0] = x;
    xl[1] = y;
    xl[2] = z;
    // LU decompsing the coeff matrix and solving for xisol

    ludcmp_( MS, &three, &three, indx2, &fnumber);
    lubksb_( MS, &three, &three, indx2, xisol);

    double tol = 0.000001;
    int truth =1;                        
  
    for( int i=0; i<3 ; i++) {         
        if ( xisol[i] > 1.0+tol || xisol[i] < 0.0-tol ) {
            truth = 0; 
        } 
    }
    double l4 = 1 - xisol[0] - xisol[1] - xisol[2];
    if ( l4 > 1.0+tol || l4 < 0.0-tol ) truth = 0;
  
    if (truth){                          
        pt[0] = xisol[0];                 
        pt[1] = xisol[1];                 
        pt[2] = xisol[2];  
    }                                    

    for(int i =0; i< 3;i++) delete [] A[i] ;
    for(int i =0; i< 3;i++) delete [] AA[i] ;
    delete [] A;
    delete [] AA;


    return truth;
}
