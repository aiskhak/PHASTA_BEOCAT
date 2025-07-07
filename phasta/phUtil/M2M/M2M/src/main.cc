#include <iostream>
#include <fstream>
#include <cstdlib>  // for NULL
#include <string.h>
#include <stdio.h>


#include "MeshSim.h"
#include "SimParasolidKrnl.h"
#include "SimDiscrete.h"

#include "attachData.h"
#include "phReadWrite.h"

#include "M2MSolutionTransfer.h"

using std::endl; 
using std::cout;
using std::cin;

using std::ifstream;

int readIndexMap(char *indexMapFileName, int indexMapDir, int* indexMap) {

  ifstream indexMapFile;
  indexMapFile.open(indexMapFileName);

  if(!indexMapFile) {
    cout<<"\n Error : Could not open index map file ["<<indexMapFileName<<"]..."<<endl;
    exit(1);
  }

  int readIndexMapSize;
  indexMapFile>>readIndexMapSize;

  int countIndexMapSize = 0;
  int read_in, read_to, read_junk;
   for(int iIndex=0; iIndex<readIndexMapSize; iIndex++) {
    if(indexMapFile.eof())
      break;

    indexMapFile>>read_in>>read_junk>>read_to;

    if(indexMapDir)
      indexMap[read_in] = read_to;
    else
      indexMap[read_to] = read_in;

    countIndexMapSize++;
  }

  if(countIndexMapSize!=readIndexMapSize) {
    cout<<"\n Error : Size mismatch for index map in file ["<<indexMapFileName<<"]..."<<endl;
    exit(1);
  }

  indexMapFile.close();

  return readIndexMapSize;
}

void mapArrayC(int *indexMap, double *in_var_array, double *out_var_array, int indexMapSize, int numVars) {

  for(int iIndex=0; iIndex<indexMapSize; iIndex++) {
    int inArrayIdx = indexMap[iIndex]*numVars;
    int outArrayIdx = iIndex*numVars;
    for(int iVar=0; iVar<numVars; iVar++) {
      out_var_array[outArrayIdx+iVar] = in_var_array[inArrayIdx+iVar];
    }
  }

}

void mapArrayF(int *indexMap, double *in_var_array, double *out_var_array, int indexMapSize, int numVars) {

  for(int iIndex=0; iIndex<indexMapSize; iIndex++) {
    int inArrayIdx = indexMap[iIndex];
    int outArrayIdx = iIndex;
    for(int iVar=0; iVar<numVars; iVar++) {
      out_var_array[outArrayIdx+iVar*indexMapSize] = in_var_array[inArrayIdx+iVar*indexMapSize];
    }
  }

}

int main(int argc, char *argv[])
{

pProgress progress;
pMConnector connector;

  int displayHelp = 0;
  for(int iArgc=1; iArgc<argc; iArgc++) {
    if(!strcmp(argv[iArgc],"-h")) {
      cout << endl;
      cout << "  HELP requested (by using \"-h\" in arguments) -- " << endl;
      displayHelp = 1;
    }
  }

  if (argc<11 || argc>13 || displayHelp) {
    cout << endl;
    cout << " usage : " << endl;
    cout << "   <executable>" << endl;
    cout << "   <from-model-name.xmt_txt or from-model-name.x_t or from-model-name.sdm> (example, geom.xmt_txt or geom.sdm)" << endl;
    cout << "   <from-mesh-name.sms> (example, geom_from.sms)" << endl;
    cout << "   <from-restart-name.sn.pn> (example, restart_source.20.1)" << endl;
    cout << "   <from-idmap.dat> [use 0 for NULL otherwise filename] (example, idmap_source.dat)" << endl;
    cout << "   <to-model-name.xmt_txt or to-model-name.x_t or to-model-name.sdm> (example, geom-2.xmt_txt or geom-2.sdm)" << endl;
    cout << "   <to-mesh-name.sms> (example, geom_to.sms)" << endl;
    cout << "   <to-restart-name.sn.pn> (example, restart_dest.20.1)" << endl;
    cout << "   <to-idmap.dat> [use 0 for NULL otherwise filename] (example, idmap_dest.dat)" << endl;
    cout << "   <region-finder-option> [0:loop each region-slow (for mixed topolgy), 1:use octree-fast (for tets. only in source mesh)]" << endl;
    cout << "   <outside-point-option> (when true(>1) then exit if a point in to-mesh is far off w.r.t. from-mesh)" << endl;
    cout << "   <half-to-full-domain-option> [optional] [+/-1 for x, +/-2 for y and +/-3 for z, reflection at zero on either axis and direction] (example, use -3 when half/source domain is in negative z)" << endl;
    cout << "   <tags-option> [optional] [n for number of tags - for those model entities to be considered (in dest.)]" << endl;
    cout << endl;
    exit(0);
  }

  char sim_lic_file[256];
  char *sim_lic_file_env = 0;
  sim_lic_file_env = getenv("SIM_LICENSE_FILE");
  if(sim_lic_file_env) {
    cout << " SIM_LICENSE_FILE is set : " << sim_lic_file_env << endl;
    sprintf(sim_lic_file,sim_lic_file_env);
  }
  else {
    strcpy(sim_lic_file,"/net/common/meshSim/license/license.txt");
  }
  Sim_readLicenseFile(sim_lic_file);

  MS_init();

  char *modelExtn = strchr(argv[1],'.');
  // parasolid or discrete
  int parasolid = 1;
  if(!strncmp(modelExtn,".sdm",4))
    parasolid = 0;

  pGModel from_model, to_model;  // declare model objects
  pMesh from_mesh, to_mesh; // declare mesh objects

  if(parasolid) {
    SimParasolid_start(1);
    cout<<"\n ";

    pNativeModel from_nmodel = 0, to_nmodel = 0;
    from_nmodel = ParasolidNM_createFromFile(argv[1],0);
    to_nmodel = ParasolidNM_createFromFile(argv[5],0);

    if(NM_isAssemblyModel(from_nmodel)) {
      pGAModel from_amodel = GAM_createFromNativeModel(from_nmodel, progress);
      NM_release(from_nmodel);
      from_model = GM_createFromAssemblyModel(from_amodel, connector, progress);
      GM_release(from_amodel);
    }
    else {
      from_model = GM_createFromNativeModel(from_nmodel, progress);
      NM_release(from_nmodel);
    }

    if(NM_isAssemblyModel(to_nmodel)) {
      pGAModel to_amodel = GAM_createFromNativeModel(to_nmodel, progress);
      NM_release(to_nmodel);
      to_model = GM_createFromAssemblyModel(to_amodel, connector, progress);
      GM_release(to_amodel);
    }
    else {
      to_model = GM_createFromNativeModel(to_nmodel, progress);
      NM_release(to_nmodel);
    }

    if(!from_model) {
      cout<<"\n Error : Could not load source model..."<<endl;
      exit(1);
    }

    if(!to_model) {
      cout<<"\n Error : Could not load destination model..."<<endl;
      exit(1);
    }

  }
  else {
    SimDiscrete_start(0);
    from_model = (pGModel)DM_load(argv[1], progress);
    to_model = (pGModel)DM_load(argv[5], progress);
  }

  from_mesh = M_new(0,from_model);
  to_mesh = M_new(0,to_model);

  cout<<endl;
  cout<<" Reading mesh from files : "<<endl;
  cout<<"   source mesh : "<<argv[2]<<endl;
  cout<<"   dest.  mesh : "<<argv[6]<<endl;

  from_mesh = M_load(argv[2],from_model, progress);
  to_mesh = M_load(argv[6],to_model, progress);

  pMeshDataId solutionID = MD_newMeshDataId("solution id");

  char fromSolutionFileName[256], toSolutionFileName[256];
  strcpy(fromSolutionFileName,argv[3]);
  strcpy(toSolutionFileName,argv[7]);

  int fromSolutionIndexMapFlag = 0, toSolutionIndexMapFlag = 0;
  int fromSolutionIndexMapDir = 1, toSolutionIndexMapDir = 0; // direction of map
  char fromSolutionIndexMapFileName[256] = "0", toSolutionIndexMapFileName[256] = "0";
  if(strcmp(argv[4],"0")) { // if not NULL/0
    fromSolutionIndexMapFlag = 1;
    strcpy(fromSolutionIndexMapFileName,argv[4]);
  }
  if(strcmp(argv[8],"0")) { // if not NULL/0
    toSolutionIndexMapFlag = 1;
    strcpy(toSolutionIndexMapFileName,argv[8]);
  }
  

  cout<<endl;
  cout<<" Reading solution from file : "<<endl;
  cout<<"   source solution : "<<fromSolutionFileName<<endl;
  if(fromSolutionIndexMapFlag)
    cout<<"   source idmap    : "<<fromSolutionIndexMapFileName<<endl;

  char fieldTag[64];
  sprintf(fieldTag,"solution");
    
  int readFromNumNodes, numSolVars, lstep;
  int fromNumVerts = M_numVertices(from_mesh);

  readParametersFromFile(fromSolutionFileName,fieldTag,
			 readFromNumNodes,numSolVars,lstep);

  if(readFromNumNodes!=fromNumVerts) {
    cout<<"\n Error : mismatch in num. of nodes from solution file ["<<readFromNumNodes<<"] and source mesh ["<<fromNumVerts<<"] "<<endl;
    cout<<endl;
    exit(1);
  }

  int poly = 1; // only supports 1st order as of now

  double *from_solution, *from_solution_mapped;
  readArrayFromFile(fromSolutionFileName,fieldTag,from_solution);

  if(fromSolutionIndexMapFlag) {
    int *fromSolutionIndexMap = new int[readFromNumNodes];
    int fromSolutionIndexMapSize = readIndexMap(fromSolutionIndexMapFileName,fromSolutionIndexMapDir,fromSolutionIndexMap);
    if(fromSolutionIndexMapSize!=readFromNumNodes) {
        cout<<"\n Error : mismatch in size of index map from source idmap ["<<fromSolutionIndexMapSize<<"] and from source restart ["<<readFromNumNodes<<"]..."<<endl;
        exit(1);
    }

    from_solution_mapped = new double[readFromNumNodes*numSolVars]; 
    mapArrayC(fromSolutionIndexMap,from_solution,from_solution_mapped,readFromNumNodes,numSolVars);
    delete [] from_solution;
    delete [] fromSolutionIndexMap;
  }
  else
    from_solution_mapped = from_solution;

  attachArray(from_solution_mapped,from_mesh,solutionID,numSolVars,poly);
  delete [] from_solution_mapped;

  int rfOption = atoi(argv[9]);
  int farOption = atoi(argv[10]);
  int halfToFullOption = 0;
  if(argc>=12)
    halfToFullOption = atoi(argv[11]);

  int entTagsOption = 0;
  if(argc>=13)
    entTagsOption = atoi(argv[12]);

  cout<<endl;
  cout<<" M2M solution transfer begin..."<<endl;
  // memory is allocated for solution on each vertex
  M2MSolutionTransfer(from_model,to_model,from_mesh,to_mesh,solutionID,numSolVars,rfOption,farOption,halfToFullOption,entTagsOption);
  cout<<endl;
  cout<<" M2M solution transfer done..."<<endl;


  cout<<endl;
  cout<<" Writing solution to file : "<<endl;
  cout<<"   dest.  solution   : "<<toSolutionFileName<<endl;
  if(toSolutionIndexMapFlag)
    cout<<"   dest.  idmap used : "<<toSolutionIndexMapFileName<<endl;
  cout<<endl;

  int toNumNodes = M_numVertices(to_mesh);

  double *to_solution, *to_solution_mapped;
  getAttachedArray(to_solution,to_mesh,solutionID,numSolVars,poly);

  if(toSolutionIndexMapFlag) {
    int *toSolutionIndexMap = new int[toNumNodes];
    int toSolutionIndexMapSize = readIndexMap(toSolutionIndexMapFileName,toSolutionIndexMapDir,toSolutionIndexMap);
    if(toSolutionIndexMapSize!=toNumNodes) {
        cout<<"\n Error : mismatch in size of index map from dest. idmap ["<<toSolutionIndexMapSize<<"] and from dest. restart ["<<toNumNodes<<"]..."<<endl;
        exit(1);
    }

    to_solution_mapped = new double[toNumNodes*numSolVars]; 
    mapArrayF(toSolutionIndexMap,to_solution,to_solution_mapped,toNumNodes,numSolVars);
    delete [] to_solution;
    delete [] toSolutionIndexMap;
  }
  else
    to_solution_mapped = to_solution;

  writeArrayToFile(toSolutionFileName,"solution","binary","write",
		   toNumNodes,numSolVars,lstep,to_solution_mapped);
  delete [] to_solution_mapped;

  cleanAttachedData(from_mesh,solutionID,0);
  cleanAttachedData(to_mesh,solutionID,0);
  MD_deleteMeshDataId(solutionID);

  M_release(from_mesh);
  M_release(to_mesh);

  GM_release(from_model);
  GM_release(to_model);

  if(parasolid)
    SimParasolid_stop(1);
  else
    SimDiscrete_stop(0);

  MS_exit();
  return 1;
}
