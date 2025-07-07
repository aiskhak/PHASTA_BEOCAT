#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include <stdio.h>
#include "phastaIO.h"
#include "MeshSim.h"
#include "SimModel.h"
#include "SimDiscrete.h"
#include "FortranIO.h"

#define PI 3.1415926535897932384626433832795
/* #define SimDeprecated80 1   */

using std::cout;
using std::endl;
using std::ofstream;
using std::setprecision;
using std::scientific;
using std::setw;
int
main (int argc, char **argv) {
    
    double angleTol=20.0;
    char geomfile[50];
    int header[8];
    
        int version;
            version = 0;

    int displayHelp = 0;
    for(int iArgc=0; iArgc<argc; iArgc++) {
      if(!strcmp(argv[iArgc],"-h")) {
        cout << endl;
        cout << "  HELP requested (by using \"-h\" in arguments) -- " << endl;
        displayHelp = 1;
      }
    }

    if (argc<2 || displayHelp) {
        cout << endl;
        cout <<" Using geombc.dat: "<< endl;
        cout <<"   "<< argv[0] <<" angleTol "<< endl;
        cout <<"   or"<<endl;
        cout <<" Using old geom.dat format: "<< endl;
        cout <<"   "<< argv[0] <<" angleTol geom.dat poly-order optionalCReadEndianFlag"<< endl;
        cout <<"   or"<<endl;
        cout <<" Using: "<<endl;
        cout <<"   "<< argv[0] <<" -f  smsfile [angleTol]"<< endl;
        cout << endl;
        exit(0);
    }
  
  char sim_lic_file[256];
  char *sim_lic_file_env = 0;
  sim_lic_file_env = getenv("SIM_LICENSE_FILE");
  if(sim_lic_file_env) {
    cout << "SIM_LICENSE_FILE is set : " << sim_lic_file_env << endl;
    sprintf(sim_lic_file,sim_lic_file_env);
  }
  else {
    sprintf(sim_lic_file,"/users/SCOREC/public/meshSim/license/license.txt");
  }

  Sim_readLicenseFile(sim_lic_file);

  Sim_logOn("createDM.log");
  MS_init();

  SimDiscrete_start(0);

  pProgress prog=Progress_new();
  Progress_setDefaultCallback(prog);

  pMesh mesh = M_new(0,0);
  pDiscreteModel model = 0;

  int geomkmapflag = 0; 
  pMeshDataId vDataID = MD_newMeshDataId("vertex labels");
  pVertex vertex;

    if( !strcmp(argv[1],"-f")) {

        mesh = M_load(argv[2], 0, prog);
        if ( argc == 4 ) angleTol = atof(argv[3]);
        
    } else {
        geomkmapflag = 0;   // 1;   Igor, 05/2012
        int igeom=1;
        long long int numVerts, numElems, nen, numnodes;
        int iarray[20];
        int ione=1, itwo=2, ithree=3,  iseven=7;
        char* iotype="binary"; 
        double* coords ; 
        int* elementData;
        int* elementType;

	cout <<" argc =   "<< argc <<"; "<< endl;
 
        if (argc>=4) {
        angleTol = atof(argv[1]);
        strcpy(geomfile,argv[2]);
        int poly_order = atoi(argv[3]);
    
        // read the nodes and connectivity
        cout <<" Calling FortranIO geom..  " << endl;
        FortranIO geom(geomfile);
        cout <<" done !  " << endl;

        cout <<" read the nodes and connectivity - started ...  " << endl;

        int cread_endian_flag = 0;
        if(argc>=5)
          cread_endian_flag = atoi(argv[4]);
        // cread_endian_flag of  1 implies C based file reads
        // cread_endian_flag of >1 implies C based file reads and byte swapping
        geom.setCReadEndianFlag(cread_endian_flag);

        geom.read( header, 7 );
    
        numVerts = header[0];
        numElems = header[2];
        nen      = header[4];

        coords   = new double[ 3*numVerts ];
	cout <<" Allocating an array with "<< numElems <<" elements; "<< endl;
        elementType = new int[ numElems ];
	long int nData = numElems*nen;
        cout <<" Allocating an array with "<< nData <<" elements ... ";
        elementData = new int[ nData ];
        cout <<" ... done "<< endl;
        if(!cread_endian_flag) { 
        cout <<" l. 133: first if "<< endl;
          geom.read( coords, numVerts, 3 );
        cout <<" l. 135: coords, numVerts are read "<< endl;
          geom.read( elementData, numElems, nen );
        cout <<" l. 137: elementData, numElems are read "<< endl;
        }
        else {
        cout <<" l. 138: else statement "<< endl;
          double* coordsT   = new double[ 3*numVerts ];
          int* elementDataT = new int[ numElems*nen ];

          geom.read( coordsT, numVerts, 3 );
          geom.read( elementDataT, numElems, nen );

          for ( int i=0; i < numElems; i++ ) {
            for (int j=0; j < nen; j++) {
              elementData[i*nen+j]=elementDataT[i+j*numElems];
            }
          }

          for ( int i=0; i < numVerts; i++ ) {
            for (int j=0; j < 3; j++) {
              coords[i*3+j]=coordsT[i+j*numVerts];
            }
          }

          delete [] coordsT;
          delete [] elementDataT;
        }
        }
    
        cout <<"  ... done " << endl;

          printf("permute nodes (and shift) from ensa mesh to mesh database numbering... \n");
          cout <<"  nen =  " << nen << endl;

        // permute nodes (and shift) from ensa mesh to mesh database numbering
        //
	  numnodes = numElems*nen;

        for ( long long int i=0; i < numnodes; i++ ) elementData[i]--;
    
        switch( nen ) {
        case 4:
            for ( int i=0; i < numElems; i++ ) {
                int tmp = elementData[i*nen+2];
                elementData[i*nen+2] = elementData[i*nen+1];
                elementData[i*nen+1] = tmp;
                elementType[ i ] = 10;
            }
            break;
        case 6:
            for ( int i=0; i < numElems; i++ ) elementType[ i ] = 12;
            break;
        case 8:
            for ( int i=0; i < numElems; i++ ) elementType[ i ] = 13;
            break;
        default:
            cout << "unknown element with "<< nen << "vertices" << endl;
            exit( 1 );
        }

         printf("permute nodes (and shift) from ensa mesh to mesh database numbering... - done ! \n");

        // use vReturn to create geom.kmap and geom.old
        pVertex *vReturn = new pVertex[numVerts];

         printf("create the mesh file from the coordinates and connectivity \n");

        // create the mesh file from the coordinates and connectivity
        printf("numVerts = %ld \n", numVerts);
        printf("numElems = %ld \n", numElems);
        if(M_importFromData( mesh,  numVerts, coords,
                             numElems, elementType, elementData, vReturn, NULL, prog )) { //check for error
          printf("Error importing mesh data\n");
          M_release(mesh);
          return 1;
        }
	printf(" Printing the debug file for simmetrix...  \n");
        const char * Sim_filename = "Debugfile_simmetrix.dat";
//        M_writeImportDataDebugFile(Sim_filename, numVerts, coords,
//                             numElems, elementType, elementData);

         printf("create the mesh file from the coordinates and connectivity - done ! \n");

        delete [] coords;
         printf(" delete [] coords - done! \n");
        delete [] elementData;
         printf("delete [] elementData - done! \n");
        delete [] elementType;
         printf("delete [] elementType - done! \n");

          printf("Importing mesh data - done! \n");

          printf("Skipping the MS_checkMeshIntersections   - fix if needed! \n");

        // check the input mesh for intersections
        // this call must occur before the discrete model is created
//        if(MS_checkMeshIntersections(mesh,0, prog)) {
//          printf("There are intersections in the input mesh\n");
//          M_release(mesh);
//          return 1;
//        }

// Here we will print out the sms file to save the total memory consumption:

        char outputMeshFile[128];

        sprintf(outputMeshFile,"geom.sms");

//        int version;
//        version = 0;

//        cout<<"Writing "<<outputMeshFile<<" (version "<<version<<")..."<<endl;
//        M_write(mesh,outputMeshFile,version, prog);
//        printf("... done! \n");

        printf("create the Discrete model - started! \n");

        // create the Discrete model
        model = (pDiscreteModel) DM_createFromMesh(mesh, 0, prog);   // was "0"; "1" means to "own" the mesh and thus to delete the original copy.
        if(!model) { //check for error
          printf("Error creating Discrete model from mesh\n");
//          M_release(mesh);
          return 1;
        };

        printf("create the Discrete model - done! \n");

//        cout<<"Writing geom.old..."<<endl;
//        ofstream geomoldfile("geom.old");
//        geomoldfile<<numVerts<<" 3 "<<endl;
//        for(int iVert=0; iVert<numVerts; iVert++) {
//          vertex = vReturn[iVert];

//          EN_attachDataInt((pEntity)vertex,vDataID,iVert); // assuming no data attached with vDataID

//          double vxyz[3];
//          V_coord(vertex,vxyz);
//          geomoldfile<<iVert+1<<" 0 "<<setprecision(16)<<scientific<<vxyz[0]<<" "<<vxyz[1]<<" "<<vxyz[2]<<endl;
//        }

//        geomoldfile.close();
        delete [] vReturn;
    }


    cout << " The angle tolerance  : "<< angleTol << endl;
    angleTol = angleTol*PI/180.0; //conversion to radians
    
    // define the Discrete model
    DM_findEdgesByFaceNormalsDegrees(model, angleTol, prog);
    DM_eliminateDanglingEdges(model, prog);
    if(DM_completeTopology(model, prog)) { //check for error
      printf("Error completing Discrete model topology\n");
      M_release(mesh);
      GM_release(model);
      return 1;
    }


    // below should come after DM_xxx APIs as mesh entities may get re-ordered in these calls
    if(geomkmapflag) {
        cout<<"Writing geom.kmap..."<<endl;
        ofstream geomkmapfile("geom.kmap");
        VIter vit = M_vertexIter(mesh);
        while(vertex=VIter_next(vit)) {

          int vID;
          if(!EN_getDataInt((pEntity)vertex,vDataID,&vID)) {
            printf("\n vertex is not labeled\n");
            return 1;
          }
          geomkmapfile<<setw(9)<<vID<<endl;
          EN_deleteData((pEntity)vertex,vDataID);
        }
        VIter_delete(vit);
        geomkmapfile.close();
    }
    MD_deleteMeshDataId(vDataID);
    

  // Print out information about the model
  printf("Number of model vertices: %d\n",GM_numVertices(model));
  printf("Number of model edges: %d\n",GM_numEdges(model));
  printf("Number of model faces: %d\n",GM_numFaces(model));
  printf("Number of model regions: %d\n",GM_numRegions(model));

  char outputModelFile[128];
  char outputMeshFile[128];
  sprintf(outputModelFile,"geom.sdm");
  sprintf(outputMeshFile,"geom.sms");

//  int version;

  version = 0;
  cout<<"Writing "<<outputModelFile<<" (version "<<version<<")..."<<endl;
  DM_write(model,outputModelFile,version, prog); // save the discrete model

//  version = 0;
  cout<<"Writing "<<outputMeshFile<<" (version "<<version<<")..."<<endl;
  M_write(mesh,outputMeshFile,version, prog);

//  char outputMeshFileV2[128];
//  char outputMeshFileV5[128];
//  sprintf(outputMeshFileV2,"geom_ver2.sms");
//  sprintf(outputMeshFileV5,"geom_ver5.sms");
//  version = 2;
//  cout<<"Writing "<<outputMeshFileV2<<" (version "<<version<<")..."<<endl;
//  M_write(mesh,outputMeshFileV2,version);
//  version = 54;
//  cout<<"Writing "<<outputMeshFileV5<<" (version "<<version<<")..."<<endl;
//  M_write(mesh,outputMeshFileV5,version);

  M_release(mesh);
//  GM_release(model);
  SimDiscrete_stop(0);
  MS_exit();
  Sim_logOff();
  Sim_unregisterAllKeys();

  return 0;
}
