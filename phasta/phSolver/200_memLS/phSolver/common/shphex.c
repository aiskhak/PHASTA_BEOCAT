/* fortran wrapper for C function HexShapeAndDrv which returns
   the shape functions and their derivatives for tets
   */

int HexShapeAndDrv(int p,double par[3],double N[],double dN[][3]);

#ifdef sun4_5
shphex_(int *p, double par[], double N[], double dN[][3])
#elif LINUX
shphex_(int *p, double par[], double N[], double dN[][3])
#elif ibm
shphex(int *p, double par[], double N[], double dN[][3])
#elif sgi
void shphex_(int *p, double par[], double N[], double dN[][3])
#elif decalp
void shphex_(int *p, double par[], double N[], double dN[][3])
#elif intel
void SHPHEX(int *p, double par[], double N[], double dN[][3])
#endif
{
  
  HexShapeAndDrv(*p,par,N,dN);

}
