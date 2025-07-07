/* fortran wrapper for C function PyrShapeAndDrv which returns
   the shape functions and their derivatives for pyramid
   */

int PyrShapeAndDrv(int p,double par[3],double N[],double dN[][3]);
/* p:      the order of the element shape function
   par[3]: xi[3]
   N[]:    shape function
   dN[][3]:derivative of shape function
   */

#ifdef sun4_5
shppyr_(int *p, double par[3], double N[], double dN[][3])
#elif LINUX
shppyr_(int *p, double par[3], double N[], double dN[][3])
#elif ibm
shppyr(int *p, double par[3], double N[], double dN[][3])
#elif sgi
void shppyr_(int *p, double par[3], double N[], double dN[][3])
#elif decalp
void shppyr_(int *p, double par[3], double N[], double dN[][3])
#elif intel
void SHPPYR(int *p, double par[3], double N[], double dN[][3])
#endif
{

  PyrShapeAndDrv(*p,par,N,dN);
  
}








