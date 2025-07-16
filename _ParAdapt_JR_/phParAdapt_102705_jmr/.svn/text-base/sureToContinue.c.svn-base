#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "SimPMesh.h"
#include "parallel.h"

void
sureToContinue(void){

    char yn[1];
    int yesno=0;
    MPI_Barrier(MPI_COMM_WORLD);
    if (PMU_rank == 0){
        printf("\nare you sure you want to continue (y/n)? ");
        scanf("%d",&yn[0]);
        if(yn[0] == 'y')
            yesno=1;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&yesno, 1, MPI_INT, 0, MPI_COMM_WORLD);
    

    if(yesno==0){
        SimPMesh_stop();
        exit(0);
    }
}
