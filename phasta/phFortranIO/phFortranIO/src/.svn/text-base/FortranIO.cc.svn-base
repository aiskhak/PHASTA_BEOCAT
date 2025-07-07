#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "FortranIO.h"

#ifdef SGI
#define readFortranInt readfint_
#define readFortranDouble readfdbl_
#define readFortranRestHead readfhd_
#define writeFortranInt writefint_
#define writeFortranDouble writefdbl_
#define writeFortranRestHead writefhd_
#define openFortran openf_
#define closeFortran closef_
#elif SUN4
#define readFortranInt readfint_
#define readFortranDouble readfdbl_
#define readFortranRestHead readfhd_
#define writeFortranInt writefint_
#define writeFortranDouble writefdbl_
#define writeFortranRestHead writefhd_
#define openFortran openf_
#define closeFortran closef_
#elif LINUX
#define readFortranInt readfint_
#define readFortranDouble readfdbl_
#define readFortranRestHead readfhd_
#define writeFortranInt writefint_
#define writeFortranDouble writefdbl_
#define writeFortranRestHead writefhd_
#define openFortran openf_
#define closeFortran closef_
#endif

#define swap_char(A,B) { ucTmp = A; A = B ; B = ucTmp; }

void 
SwapArrayByteOrder_( void*  array, 
                     size_t nbytes, 
                     int    nItems ) {
    /* This swaps the byte order for the array of nItems each
       of size nbytes , This will be called only locally  */
    int i,j;
    unsigned char ucTmp;
    unsigned char* ucDst = (unsigned char*)array;
        
    for(i=0; i < nItems; i++) {
        for(j=0; j < (nbytes/2); j++)
            swap_char( ucDst[j] , ucDst[(nbytes - 1) - j] );
        ucDst += nbytes;
    }
}

extern "C" {
  void readFortranInt(int *unit, int *ar,
			       int *n1, int *n2, int *nth);
  void readFortranDouble(int *unit, double *ar,
				  int *n1, int *n2, int *nth);
  void readFortranRestHead(int *unit, int *ar);

  void writeFortranInt(int *unit,int *ar,int *n1,int *n2);
  void writeFortranDouble(int *unit,double *ar,int *n1,int *n2);
  void writeFortranRestHead(int *unit,int *ar);

  void openFortran(char *fname, int *unit, int len);
  void closeFortran(int *unit);
}

int FortranIO::global_unit = 15;

FortranIO::FortranIO(const char *name)
{
  length = (int)strlen(name);
  fname = new char[length+1];
  strcpy(fname,name);
  CReadFlag = 0;
  EndianFlag = 0;
  inFile = 0;
  unit = global_unit++;
  open();
}

FortranIO::~FortranIO()
{
  delete[]fname;
  close();
  if(inFile)
    fclose(inFile);
}

FortranIO& FortranIO::open()
{
  openFortran(fname,&unit,length);
  
  return *this;
}

void FortranIO::close()
{
  closeFortran(&unit);
}

void FortranIO::setCReadEndianFlag(int flag)
{
  if(flag)
    CReadFlag = 1;
  if(flag>1)
    EndianFlag = 1;
}

FortranIO& FortranIO::read(int *ar, int n1, int n2, int nth)
{
  if(CReadFlag)
    readCInt(ar,n1,n2,nth);
  else
    readFortranInt(&unit,ar,&n1,&n2,&nth);

  return *this;
}

FortranIO& FortranIO::read(double *ar, int n1, int n2, int nth)
{
  if(CReadFlag)
    readCDouble(ar,n1,n2,nth);
  else
    readFortranDouble(&unit,ar,&n1,&n2,&nth);

  return *this;
}

void FortranIO::readCInt(int *ar, int n1, int n2, int nth)
{
  readCData((void *)ar,n1,n2,nth,sizeof(int));
}

void FortranIO::readCDouble(double *ar, int n1, int n2, int nth)
{
  readCData((void *)ar,n1,n2,nth,sizeof(double));
}

void FortranIO::readCData(void *ar, int n1, int n2, int nth, size_t type_size)
{
  int openflag = 0;
  if(!inFile) {
    inFile = fopen(fname,"r");
    openflag = 1;
  }

  // assuming nth-bytes
  fseek(inFile,nth*8,SEEK_CUR);

  // assuming file was written using Fortran IO
  int int_value;
  fread(&int_value,sizeof(int),1,inFile);
  if(!openflag) {
    // another one when file not opened in this call
    fread(&int_value,sizeof(int),1,inFile);
  }

  fread(ar,type_size,n1*n2,inFile);
  if(EndianFlag) SwapArrayByteOrder_(ar,type_size,n1*n2);
}

FortranIO& FortranIO::readRestartHeader(int &num, int &step)
{
  int ar[2];

  readFortranRestHead(&unit,ar);
  num  = ar[0];
  step = ar[1];

  return *this;
}

FortranIO& FortranIO::readRestartHeader(int &num)
{
  int ar[2];

  readFortranRestHead(&unit,ar);
  num  = ar[0];

  return *this;
}

FortranIO& FortranIO::write(int *ar, int n1, int n2)
{
  writeFortranInt(&unit,ar,&n1,&n2);

  return *this;
}

FortranIO& FortranIO::write(double *ar, int n1, int n2)
{
  writeFortranDouble(&unit,ar,&n1,&n2);

  return *this;
}

FortranIO& FortranIO::writeRestartHeader(int num, int step)
{
  int ar[2]={num,step};

  writeFortranRestHead(&unit,ar);

  return *this;

}
