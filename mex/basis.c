#include "basis.h"
#include <math.h>
#include <stdlib.h>


void computeMPMBasis1D (double x, double h, double* f, double* df)
{
   if ( fabs(x) <= h ){
       *f  = 1.0 - fabs(x)/h;
       double sigx = ( x < 0 ) ? -1. : 1.;
       *df = -sigx/h;
   }
   else{
       *f  = 0.;
       *df = 0.;
   }
}

void computeMPMBasis2D (double* x, double* h, double* f, double* df1, double* df2)
{
   /* compute the 1D shape functions */

   double fx,fy,dfx,dfy;

   computeMPMBasis1D ( x[0], h[0], &fx, &dfx );
   computeMPMBasis1D ( x[1], h[1], &fy, &dfy );


   /* compute the 2D shape functions as tensor products*/

   *f     = fx  * fy;
   *df1   = dfx * fy;
   *df2   = fx  * dfy;
}

void computeGIMPBasis1D (double x, double h, double lp, double* f, double* df)
{
   double lp2 = lp/2.;

   if ( fabs(x) < lp2 ){
       *f  = 1.0 - (4*x*x+lp*lp)/(4.*h*lp);
       *df = -8.*x/(4.*h*lp);
   }
   else if  ( fabs(x) < h - lp2){
       *f  = 1-fabs(x)/h;
       double sigx = ( x < 0 ) ? -1. : 1.;
       *df = -1/h*sigx;
   }
   else if  ( fabs(x) < h + lp2){
       double tem =  h+lp2-fabs(x);
       *f  = tem*tem/(2.*h*lp);
       double sigx = ( x < 0 ) ? -1. : 1.;
       *df = -tem/(h*lp)*sigx;
   }
   else{
       *f  = 0;
       *df = 0;
   }
}

void computeGIMPBasis2D (double* x, double* h, double *lp, double* f, double* df1, double* df2)
{
   /* compute the 1D shape functions*/

   double fx,fy,dfx,dfy;

   computeGIMPBasis1D ( x[0], h[0], lp[0], &fx, &dfx );
   computeGIMPBasis1D ( x[1], h[1], lp[1], &fy, &dfy );

   /* compute the 2D shape functions as tensor products */

   *f     = fx  * fy;
   *df1   = dfx * fy;
   *df2   = fx  * dfy;
}
