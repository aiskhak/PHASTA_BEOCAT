#include <stdlib.h>
#include <strings.h>
#include <iostream>

#include "MeshSim.h"
#include "SimParasolidKrnl.h"

using std::cout;
using std::endl;
using std::cerr;

int main(int argc, char *argv[])
{

  int displayHelp = 0;
  for(int iArgc=0; iArgc<argc; iArgc++) {     if(!strcmp(argv[iArgc],"-h")) {
      cout << endl;
      cout << "  HELP requested (by using \"-h\" in arguments) -- " << endl;
      displayHelp = 1;
    }
  }

  if (argc!=5 || displayHelp) {
    cout << endl;
    cout << " usage : " << endl;
    cout << "   <executable-name>" << endl;
    cout << "   <model-name.xmt_txt>" << endl;
    cout << "   <in-mesh-name.sms>" << endl;
    cout << "   <out-mesh-name.sms>" << endl;
    cout << "   <out-mesh-sms-version>" << endl;
    cout << endl;
    exit(0);
  }


  pNativeModel nmodel = 0;
  pGModel model = 0;
  pMesh mesh;     // declare a mesh object

  char sim_lic_file[256];
  char *sim_lic_file_env = 0;
  sim_lic_file_env = getenv("SIM_LICENSE_FILE");
  if(sim_lic_file_env) {
    cout << "SIM_LICENSE_FILE is set : " << sim_lic_file_env << endl;
    sprintf(sim_lic_file,sim_lic_file_env);
    Sim_readLicenseFile(sim_lic_file);
  }
  else {
    Sim_readLicenseFile("/net/common/meshSim/license/license.txt");
  }

  MS_init();      // initialize the MeshSim

  char model_file[1024];
  strcpy(model_file,argv[1]);
  char in_mesh_file[1024];
  strcpy(in_mesh_file,argv[2]);
  char out_mesh_file[1024];
  strcpy(out_mesh_file,argv[3]);

  cout<<endl;
  cout<<" Reading... "<<endl;
  cout<<"  Model from file : "<<model_file<<endl;
  cout<<"  Mesh from file  : "<<in_mesh_file<<endl;
  cout<<endl;

  int modeler = -1;
  char ext[8];
  strcpy(ext,model_file+(strlen(model_file)-4));

  SimParasolid_start(1);
  if(strcmp(ext,"_txt") == 0 || strcmp(ext,"_TXT") == 0 || strcmp(ext,".x_t") == 0) {
    modeler = 0;
    cout << "this is a parasolid model\n";
    nmodel = ParasolidNM_createFromFile(model_file, 0);
  }
  else {
    cerr << "Didn't load a model, check that code was compiled with the correct MODELER specified and that the model file has an extension such as .xmt_txt, or .XMT_TXT, or .x_t" << endl;
    return 0;
  }

   if (NM_isAssemblyModel(nmodel)) {
    pGAModel amodel = GAM_createFromNativeModel(nmodel);
    NM_release(nmodel);
    model = GM_createFromAssemblyModel(amodel);
    GM_release(amodel);
    nmodel = GM_nativeModel(model);
    //NM_write(nmodel, "block-nonman.sat");
  }
  else
    model = GM_createFromNativeModel(nmodel);

  if(!model) {
    cerr << "Didn't load a model, check that code was compiled with the correct MODELER specified and that the model file has an extension such as .xmt_txt, or .XMT_TXT, or .x_t" << endl;
    return 0;
  }

  mesh = M_load(in_mesh_file, model);

  cout<<" Reading model and mesh done..."<<endl;
  cout<<endl;

  // default is latest (i.e., 0)
  // check Simmetric documentation for more details
  int smsversion = atoi(argv[4]);

  cout<<" Writing mesh (in version : "<<smsversion<<") into file : "<<out_mesh_file<<endl;
  cout<<endl;

  M_write(mesh,out_mesh_file,smsversion);

  M_release(mesh);
  GM_release(model);
  // NM_release(nmodel);

  SimParasolid_stop(1);

  MS_exit();

  Sim_unregisterAllKeys();

  return 1;
}
