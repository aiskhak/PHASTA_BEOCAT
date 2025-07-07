#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <iostream>
/* 
   The code calling this has glued together the solution and geometry
   files and has passed in the following:

   nv is the number of degrees of freedom per node
   numel is the total number of elements
   ien is the TRUE global connectivity  
       (index local element node number and element
        number to TRUE global node number)
   numnp is the total number of nodes                          
   x is the coordinates (3*numnp)
   q is the solution (nv*numnp)
   nen is the number of nodes in an element
*/
using namespace std;
void wrtc_(int *nv, int *nverr, int *nvavg,
           int *numel, int ien[], int part[],
           int *numnp, float x[], float q[],
           float a[], float eres[], float edif[],
           float emsf[], float avg[],
           int* nen, int clipPlane, float slope, int procnum)
{
  int i,j,count;
  char *eltype;			/* element type */

  char filename[100];
  bzero((void*)filename,100);
  sprintf(filename,"connect%d.bin",procnum);
  FILE *f1  = fopen(filename,"wb");
  bzero((void*)filename,100);
  sprintf(filename,"points%d.bin" ,procnum);
  FILE *f2  = fopen(filename,"wb");
  bzero((void*)filename,100);
  sprintf(filename,"data%d.bin"   ,procnum);
  FILE *f3  = fopen(filename,"wb");
  bzero((void*)filename,100);
  sprintf(filename,"adata%d.bin"  ,procnum);
  FILE *f3a = fopen(filename,"wb");
  bzero((void*)filename,100);
  sprintf(filename,"eres%d.bin"  ,procnum);
  FILE *f3eres = fopen(filename,"wb");
  bzero((void*)filename,100);
  sprintf(filename,"edif%d.bin"  ,procnum);
  FILE *f3edif = fopen(filename,"wb");
  bzero((void*)filename,100);
  sprintf(filename,"emsf%d.bin"  ,procnum);
  FILE *f3emsf = fopen(filename,"wb");
  bzero((void*)filename,100);
  sprintf(filename,"avrg%d.bin"  ,procnum);
  FILE *f3avg = fopen(filename,"wb");
  bzero((void*)filename,100);
  sprintf(filename,"header%d.dx"  ,procnum);
  FILE *f4  = fopen(filename,"w");
  bzero((void*)filename,100);
  sprintf(filename,"partit%d.bin" ,procnum);
  FILE *fpart  = fopen(filename,"wb");

  /* element type */
  if (*nen == 4) eltype = "tetrahedra";
  else eltype = "cubes";
     
  /* write the positions array */
  count=0;
  for(i=0; i < *numnp; i++){
    for(j=0; j < 3; j++) {
      fwrite((void *)&(x[count]),sizeof(float),1,f2);
      count++;
    }
  }

  /* write the connections array */
  float ctr[3];
  float xyz[3];
  int nn;
  int numelien=0;
  bool firstAbove, coincident, straddles;
  //for(i=0; i < *numel; i++) {
  //  ctr[0]=0.0; ctr[1]=0.0; ctr[2]=0.0;
  //  for(j=0;j<*nen;j++){
  //    nn=ien[i*(*nen)+j];
  //    ctr[0]+=x[nn*3+0];
  //    ctr[1]+=x[nn*3+1];
  //    ctr[2]+=x[nn*3+2];
  //  }
  //  ctr[0]/=(*nen); ctr[1]/=(*nen); ctr[2]/=(*nen);
  //  if(ctr[1]>slope*ctr[0]||!clipPlane){
  //    numelien++;
  //    fwrite((void *)&(part[i]),sizeof(int),1,fpart);
  //    for(j=0; j < *nen; j++) {
  //      fwrite((void *)&(ien[i*(*nen)+j]),sizeof(int),1,f1);
  //    }
  //  }
  //}
  for(i=0; i < *numel; i++) {
    straddles=false; //assume element doesnt straddle plane
    coincident=false; //assume no coincident nodes
    firstAbove=false; //assume first node not above plane
    nn=ien[i*(*nen)+0]; //first node
    if(x[nn*3+1]>slope*x[nn*3+0])firstAbove=true;
    if(x[nn*3+1]==slope*x[nn*3+0])coincident=true;
    for(j=1;j<*nen;j++){
      nn=ien[i*(*nen)+j];
      if(x[nn*3+1]>slope*x[nn*3+0]){
        if(!firstAbove) straddles=true;
      }
      if(x[nn*3+1]<slope*x[nn*3+0]){
        if(firstAbove) straddles=true;
      }
      if(x[nn*3+1]==slope*x[nn*3+0]){
        coincident=true;
      }
    }
    if(straddles||coincident||!clipPlane){
      numelien++;
      fwrite((void *)&(part[i]),sizeof(int),1,fpart);
      for(j=0; j < *nen; j++) {
        fwrite((void *)&(ien[i*(*nen)+j]),sizeof(int),1,f1);
      }
    }
  }

/* write the data array, rho,u1,u2,u3,T */
  count = 0;
  for(i=0; i< *numnp; i++){
    for(j=0; j < *nv; j++){
      fwrite((void *)&(q[count]),sizeof(float),1,f3);
      count++;
    }
  }
  
/* write the acceleration data array */
  count = 0;
  for(i=0; i< *numnp; i++){
    for(j=0; j < *nv; j++){
      fwrite((void *)&(a[count]),sizeof(float),1,f3a);
      count++;
    }
  }
  
/* write the error array */
  count = 0;
  for(i=0; i< *numnp; i++){
    for(j=0; j < 3; j++){
      fwrite((void *)&(eres[count]),sizeof(float),1,f3eres);
      count++;
    }
  }
  count = 0;
  for(i=0; i< *numnp; i++){
    for(j=0; j < 3; j++){
      fwrite((void *)&(edif[count]),sizeof(float),1,f3edif);
      count++;
    }
  }
  count = 0;
  for(i=0; i< *numnp; i++){
    for(j=0; j < 4; j++){
      fwrite((void *)&(emsf[count]),sizeof(float),1,f3emsf);
      count++;
    }
  }

/* write the average array */
  count = 0;
  for(i=0; i< *numnp; i++){
    for(j=0; j < *nvavg; j++){
      fwrite((void *)&(avg[count]),sizeof(float),1,f3avg);
      count++;
    }
  }

  /*  write the header file */
  fprintf(f4,"object 1 class array type float rank 1 shape 3 items %d msb binary\n",*numnp);
  fprintf(f4," data file %s/points%d.bin,0 \n",getcwd(NULL,128),procnum);
  fprintf(f4," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f4,"object 2 class array type int rank 1 shape %d items %d msb binary\n",*nen,numelien);
  fprintf(f4," data file %s/connect%d.bin,0 \n",getcwd(NULL,128),procnum);
  fprintf(f4," attribute \"element type\" string \"%s\" \n",eltype);
  fprintf(f4," attribute \"ref\" string \"positions\" \n\n");

  fprintf(f4,"object 3 class array type float rank 1 shape %d items %d msb binary\n",*nv,*numnp);
  fprintf(f4," data file %s/data%d.bin,0 \n",getcwd(NULL,128),procnum);
  fprintf(f4," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f4,"object 4 class array type int rank 0 items %d msb binary\n", numelien);
  fprintf(f4," data file %s/partit%d.bin \n",getcwd(NULL,128),procnum);
  fprintf(f4," attribute \"dep\" string \"connections\" \n\n");

  fprintf(f4,"object \"irregular positions irregular connections\" class field\n");
  fprintf(f4," component \"positions\" value 1\n");
  fprintf(f4," component \"connections\" value 2\n");
  fprintf(f4," component \"data\" value 3\n");
  fprintf(f4," component \"part\" value 4\n");
  fprintf(f4,"\n end \n");

  printf("\n Files successfully written... \n");
  printf("\n Number of elements = %d \n",*numel);
  printf(" Number of nodes = %d \n\n",*numnp);


  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
}
