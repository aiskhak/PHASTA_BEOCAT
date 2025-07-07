/* fortran wrapper for C function WedgeShapeAndDrv which returns
   the shape functions and their derivatives for wedge
   */

int WedgeShapeAndDrv(int p,double par[3],double N[],double dN[][3]);

#ifdef sun4_5
shp6w_(int *p, double par[], double N[], double dN[][3])
#elif LINUX
shp6w_(int *p, double par[], double N[], double dN[][3])
#elif ibm
shp6w(int *p, double par[], double N[], double dN[][3])
#elif sgi
void shp6w_(int *p, double par[], double N[], double dN[][3])
#elif decalp
void shp6w_(int *p, double par[], double N[], double dN[][3])
#elif intel
void SHP6W(int *p, double par[], double N[], double dN[][3])
#endif
{
  WedgeShapeAndDrv(*p,par,N,dN);

}








