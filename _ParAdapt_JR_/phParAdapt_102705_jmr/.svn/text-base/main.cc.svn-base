#include <stdio.h>
#include <strstream>
#include <unistd.h>

#include "phParAdapt.h"
#include "func.h"
#include "ccfunc.h"

//  #if ( defined SGI )
//  #include <stream.h>
//  #endif

#if ( defined MODELER_PARASOLID )
#include "SimParasolidKrnl.h"
#elif ( defined MODELER_DISCRETE )
#include "SimDiscrete.h"
#endif

using std::strstream;

extern "C" int readLicenseFile(char*);

// catches the debugger
// to run debugger: setenv (=export) catchDebugger=1
void
catchDebugger() {
    static volatile int debuggerPresent =0;
    while (!debuggerPresent ); // assign debuggerPresent=1
}

#define MkPschema(x) "P_SCHEMA=" #x
#define XMkPschema(x) MkPschema(x)

int main(int argc, char *argv[])
{

    putenv("PARASOLID=/usr/local/parasolid/16.0");
    putenv("P_LIST=/usr/local/parasolid/16.0/lispdata");
    putenv("P_SCHEMA=/usr/local/parasolid/16.0/schema");

    MS_init();

    SimPMesh_start(&argc, &argv);
//    PMU_start(&argc, &argv);

    SimMeshing_start();

    SimModel_start();
    // read the model and mesh
#if  ( defined MODELER_PARASOLID )
    if(SimParasolid_start(1) != 0 ){
        if(PMU_rank()==0){
            printf("\nerror in SimParasolid_start\n");
        }
        exit(0);
    }

#elif ( defined MODELER_DISCRETE )
    // first initialize !!!
    SimDiscrete_start(0);
#endif

    
//     char log_file[128];
//     sprintf(log_file,"MeshSim.%i.log",PMU_gid(PMU_rank(),0)+1);
//     MS_logOn(log_file);
    
    // to run debugger: setenv (=export) catchDebugger=1
    // 
    if ( getenv( "catchDebugger" ) ) {

        int parent_pid = getpid();
        int gdb_child = fork();

        if( gdb_child == 0 ) {
     
            printf("Debugger Process initiating\n");
            strstream exec_string;

#if ( defined decalp )
            exec_string <<"xterm -e idb " 
                        << " -pid "<< parent_pid <<" "<< argv[0] << "\n";
#endif
#if ( defined LINUX )
            exec_string <<"xterm -e idb " 
                        << " -pid "<< parent_pid <<" "<< argv[0] << "\n";
#endif
#if ( defined SUN4 )
            exec_string <<"xterm -e dbx " 
                        << " - "<< parent_pid <<" "<< argv[0] << "\n";
#endif
#if ( defined SGI )
            exec_string <<"xterm -e dbx " 
                        << " -p "<< parent_pid <<" "<< argv[0] << "\n";
#endif
            system( exec_string.str() );
            exit(0);
        }
        catchDebugger();
    }

    // readLicenseFile(NULL);
    MS_registerKey("sci001 surface 20070115 0 vsTpU4h/jQQaHGITz45ZqQ==");
    MS_registerKey("sci001 volume 20070115 0 eaWLTM+m2+GyxmM9EsdFdQ==");
    MS_registerKey("sci001 adv 20070115 0 PBYi0jWj9zpRVicjYBRGAQ==");
    MS_registerKey("sci001 parasolid 20070115 0 7VR32Xbj6a1ZLFcRkMs75A==");
    MS_registerKey("sci001 adapt 20070115 0 oIYECUfjMX1SrefEtHnHlw==");
    MS_registerKey("sci001 attributes 20070115 0 flrqhA3ph/t3KTuDaEcyrQ==");
    MS_registerKey("sci001 discrete 20070115 0 fmQ50HxgMB/UpCjv8j7Ryg==");
    MS_registerKey("sci001 pmesh 20070115 0 iX4nWNeA4GXSm3YQqUf56Q==");

#if defined(PARALLEL) && defined(DEBUG)
  int t;
  if (sscanf(argv[argc-1], "%d", &t) == 1)
    sleep(t);
#endif

    switchAdapt_Preproc(argc, argv);

#ifdef MODELER_DISCRETE 
    SimDiscrete_stop(0);
#endif
            
#ifdef MODELER_PARASOLID
    int sflag ;
    if(( sflag = SimParasolid_stop(1) ) != 0){
        if(PMU_rank()==0){
            printf("\nerror in SimParasolid_stop(1)\n"); 
            printf("\nencountered error %d\n",sflag);
        }   
        exit(0);
         
    }   
#endif
    SimModel_stop(); 

    SimMeshing_stop();

    SimPMesh_stop();

    MS_exit();

    return 0;
}
