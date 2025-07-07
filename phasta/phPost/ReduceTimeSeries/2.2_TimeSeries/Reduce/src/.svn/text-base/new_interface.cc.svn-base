/* This file provides interface functions for 'partial ' random 
   access into the PHASTA input files 

   Anil Karanam March 2001 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#define Read_Header read_header_
#define Read_Real_Array read_real_array_
#define Read_Integer_Array read_int_array_
#define Initialize_Input  initialize_input_
#define Finalize_Input  finalize_input_
#define Write_Restart write_restart_
#define Write_TimeAvg write_timeavg_
#define swap_char(A,B) { ucTmp = A; A = B ; B = ucTmp; }

static int swapgeombc=-1;
static int swaprestart=-1;
FILE* fgeombc;
FILE* frest;

/* this function is for comparing 2 strings in a case and space insensitive 
   way */ 

static int cscompare( const char* s1, const char* s2)
{
  while( *s1 == ' ') s1++;
  while( *s2 == ' ') s2++;
  while((*s1) && (*s2) && (tolower(*s1++)==tolower(*s2++) ))
    {
      while( *s1 == ' ') s1++;
      while( *s2 == ' ') s2++;
    }
  if (!(*s2)) return 1;
  else return 0;
}

static void SwapArrayByteOrder(void* array, int nbytes, int nItems)
{
  /* This swaps the byte order for the array of nItems each
     of size nbytes , This will be called only locally  */

  int i,j;
  unsigned char  ucTmp;
  unsigned char* ucDst = (unsigned char*)array;

  for(i=0; i < nItems; i++) {
    for(j=0; j < (nbytes/2); j++)
      swap_char( ucDst[j] , ucDst[(nbytes - 1) - j] );
    ucDst += nbytes;
  }
}

void Read_Header(int* If, const char* phrase, int* array, int* expect,
                 int* ierr)
{
  /* Searching for a given phrase and reading the corresponding
     header and setting the file position to the right point */

  FILE* finput;
  static int niblk=0, nbblk=0;
  int lbctr=0, *lbskip, cblk=1,tmpint; /* local block counter */
  char asciiline[255];
  char *header, *token;
  int real_length, *swapbytes, i, found=0 , zero=0, rec=0;
  long start;
  
  switch ( *If ) {
  case 1:
    finput = fgeombc;
    swapbytes = &swapgeombc;
    break;
  case 2:
    finput = frest;
    swapbytes = &swaprestart;
    break;
  }    
  
  if (cscompare(phrase,"connectivity boundary")) lbskip = &nbblk;
  else if (cscompare(phrase,"connectivity interior")) lbskip = &niblk;
  else { lbskip = &zero ; cblk = 0; }
  
  *ierr = 0 ;                /* hopeful , nonzero -> error */
  start = ftell( finput );  /* store current location of the stream */
  if(!fgets(asciiline, 255, finput)) /* get current line */
    if(feof(finput)){
      rewind(finput); rec++;
      fgets(asciiline, 255, finput);
    }
  /* here we implicitly assume that the last binary read has been complete
     and we are currently at an ascii line */
 
  while (!found && (start != ftell( finput ))) { 
    if(asciiline[0]!='\n') {
      if( real_length = strcspn(asciiline,"#")){ /* chopping comments */
	header = (char *)malloc(real_length+1);
	strncpy(header, asciiline, real_length);
	header[real_length] = NULL;
	token = strtok(header,":");
	if( cscompare(token, phrase) && (++lbctr > *lbskip )){
	  /* we found the header we are looking for */
	  found = 1 ;
	  /* if we are  counting blocks we need to store the count */
	  if(cblk) *lbskip++;
	  /* read the skipsize */
	  token = strtok(NULL," ,;<>");
	  /* now we gather all the integer headers we are expecting */
	  for(i=0; i< *expect; i++) {
	    token = strtok(NULL," ,;<>");
	    if(token) array[i]=atoi(token);  /* error */
	  }
	  /* ideally now the stream is positioned just right for reading
	     the following data block binary/ascii*/
	}else if(cscompare(token,"byteorder magic number")){
	  fread((void *)&tmpint,sizeof(int),1,finput);
	  if( 362436 != tmpint ) *swapbytes=1;
	  else *swapbytes = 0;
	  fscanf(finput,"\n");
	}else { /* just need to skip ahead */
	  token = strtok(NULL," ,;<>");
	  i = atoi(token);
	  fseek(finput,i,SEEK_CUR);
	}
      }
    }
    
    if(!found) /*as long as we havent found it yet*/
      if(!fgets(asciiline, 255, finput)) /* get next line */
	if(feof(finput)){
	  if( rec > 0 ) break; /* rewinding more than once */
	  rewind(finput);
          rec++;
	  fgets(asciiline, 255, finput);
	}
  }
  
  if(!found) {
    fprintf(stderr,"data block \"%s\" not found \n", phrase);
    *ierr=1;
    return ;
  }
}

void Read_Real_Array(int* xyzF, double* array, int* size)
{
  int* swapbytes;
  FILE* finput;
  
  switch( *xyzF ){
  case 1:
    finput = fgeombc;
    swapbytes = &swapgeombc;
    break;
  case 2:
    finput = frest;
    swapbytes = &swaprestart;
    break;
  }

  fread((void *)array, sizeof(double), *size, finput); 
  fscanf(finput,"\n");

  if ( *swapbytes ) 
    SwapArrayByteOrder((void *)array, sizeof(double), *size);
}

void Read_Integer_Array(int* xyzF, int* array, int* size)
{
  int* swapbytes;
  FILE* finput;

  switch( *xyzF ){
  case 1:
    finput = fgeombc;
    swapbytes = &swapgeombc;
    break;
  case 2:
    finput = frest;
    swapbytes = &swaprestart;
    break;
  }

  fread((void *)array, sizeof(int), *size, finput); 
  fscanf(finput,"\n");

  if ( *swapbytes ) 
    SwapArrayByteOrder((void *)array, sizeof(int), *size);
}

void Initialize_Input(int* rank, int* stepno)
{
  char gname[30], rname[30];

  sprintf(gname,"geombc.dat.%d",*rank+1);
  sprintf(rname,"restart.%d.%d",*stepno,*rank+1);
  
  fgeombc = fopen(gname,"r");
  frest = fopen(rname,"r");
}

void Finalize_Input(void)
{
  fclose(fgeombc);
  fclose(frest);
}

void Write_Restart(int* pid, int* lstep, int* nshg, int* numVars,
		   double* array1, double* array2)
{
  char rfile[30];
  FILE* fpr;
  int one=1;
  int sizeofarray = (*nshg)*(*numVars);
  int magic_number = 362436;
  int* mptr = &magic_number;
  sprintf(rfile,"restart.%d.%d",*lstep,*pid+1);
  fpr = fopen(rfile,"w");

  /* writing the top ascii header for the restart file */
   
  fprintf(fpr,"# PHASTA restart file version 1.0");
  fprintf(fpr,"\n");
  fprintf(fpr,"# num modes :%d  num var: %d step: %d \n", 
	  *nshg, *numVars, *lstep);
  fprintf(fpr,"byteorder magic number : %d\n", (int)sizeof(int));
  fwrite((void *)mptr, sizeof(int), one, fpr);
  fprintf(fpr,"\n");
  fprintf(fpr,"solution : < %d > %d %d %d \n", 
	  (*nshg)*(*numVars)*((int)sizeof(double))+ (int)sizeof(char),
	  *nshg, *numVars, *lstep);

  /* writing the q array in binary format */
  /* we always write files out in native format and expect the next */
  /* thing to use it,to convert to the byte order format of its choice */
   
  fwrite((void *)array1, sizeof(double), sizeofarray, fpr);
  fprintf(fpr,"\n");

  /* same for time derivative of the solution */

  fprintf(fpr,"time derivative of solution : < %d > %d %d %d \n", 
	  (*nshg)*(*numVars)*((int)sizeof(double))+(int)sizeof(char),
	  *nshg, *numVars, *lstep);

  fwrite((void *)array2, sizeof(double), sizeofarray, fpr);
  fprintf(fpr,"\n");

  fclose(fpr);
}

void Write_TimeAvg(int* start, int* skip, int* stop, int* part,
                   int* nshg,  int* numVars, double* array1)
{
  char rfile[40];
  FILE* fpr;
  int one=1;
  int sizeofarray = (*nshg)*(*numVars);
  int magic_number = 362436;
  int* mptr = &magic_number;
  sprintf(rfile,"restart.%d.%d.%d.%d",*start,*skip,*stop,*part+1);
  //sprintf(rfile,"restart.%d.%d",666,*part+1);
  fpr = fopen(rfile,"w");

  /* writing the top ascii header for the restart file */

  fprintf(fpr,"# PHASTA restart file version 1.0");
  fprintf(fpr,"# num modes :%d  num var: %d step: %d \n",
          *nshg, *numVars, *stop);
  fprintf(fpr,"byteorder magic number : %d\n", sizeof(int));
  fwrite((void *)mptr, sizeof(int), one, fpr);
  fprintf(fpr,"\n");
  fprintf(fpr,"solution : %d %d %d %d \n",
          (*nshg)*(*numVars)*(sizeof(double))+sizeof(char),
          *nshg,  *numVars, *stop );

  /* writing the q array in binary format */
  /* we always write files out in native format and expect the nex\
     t */
  /* thing to use it,to convert to the byte order format of its ch\
     oice */

  fwrite((void *)array1, sizeof(double), sizeofarray, fpr);
  fprintf(fpr,"\n");

  fclose(fpr);
}

