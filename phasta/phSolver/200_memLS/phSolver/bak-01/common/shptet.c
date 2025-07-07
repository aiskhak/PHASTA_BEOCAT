/* fortran wrapper for C function TetShapeAndDrv which returns
   the shape functions and their derivatives for tets
   */

int TetShapeAndDrv(int p,double par[3],double N[],double dN[][3]);

#ifdef sun4_5
shptet_(int *p, double par[], double N[], double dN[][3])
#elif LINUX
shptet_(int *p, double par[], double N[], double dN[][3])
#elif ibm
shptet(int *p, double par[], double N[], double dN[][3])
#elif sgi
void shptet_(int *p, double par[], double N[], double dN[][3])
#elif decalp
void shptet_(int *p, double par[], double N[], double dN[][3])
#elif intel
void  SHPTET(int *p, double par[], double N[], double dN[][3])
#endif
{
 TetShapeAndDrv(*p,par,N,dN);

}
