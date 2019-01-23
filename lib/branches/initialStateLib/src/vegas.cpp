//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/vegas.cpp $
//$LastChangedDate: 2013-10-28 20:27:18 +0100 (Mo, 28. Okt 2013) $
//$LastChangedRevision: 1528 $
//$LastChangedBy: gallmei $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#include <iostream>
#include <math.h>

#include "vegas.h"
#include "random.h"

#define ALPH 1.5
#define NDMX 20
#define MXDIM 10
#define TINY 1.0e-30

#define SQR(x) x*x

using namespace std;


void vegas(int ndim, integrand& fxn, double *tgral, double *sd, double *chi2a)
{

  static int k,it,mds,nd,ndo,ng,npg,ia[MXDIM+1],kg[MXDIM+1];
  static double calls,dv2g,dxg,f,f2,f2b,fb,rc,wgt,xjac,xn,xnd,xo;
  static double d[NDMX+1][MXDIM+1],di[NDMX+1][MXDIM+1],dt[MXDIM+1],
    dx[MXDIM+1],r[NDMX+1],x[MXDIM+1],xin[NDMX+1],
    xi[MXDIM+1][NDMX+1];
  long ncall=100; //set according to ndim below...
  int itmx=2;     //dito
  double schi,si,swgt;


  if(ndim==1){
    ncall=100000;
    itmx=6;
  }
  if(ndim==2){
    //   ncall=200000;
    //   itmx=10;
    ncall=100;
    itmx=2;
  }
  if(ndim==3){
    ncall=10000;
    itmx=8;
  }
  if(ndim==4){
    //     ncall=10000;
    //     itmx=8;
    //     ncall=1000000;
    //     itmx=20;
    ncall=10000;
    itmx=8;
  }

  mds=1;//change to mds=0 to disable stratified sampling.
  ndo=1;
  for (int j=1;j<=ndim;j++){
    dx[j]=1.0;
    xi[j][0]=0.0;
    xi[j][1]=1.0;
  }

  si=swgt=schi=0.0;

  nd=NDMX;
  ng=1;
  if (mds){
    ng=(int)pow(ncall/2.0+0.25,1.0/ndim);
    mds=1;
    if ((2*ng-NDMX) >= 0){
      mds=-1;
      npg=ng/NDMX+1;
      nd=ng/npg;
      ng=npg*nd;
    }
  }
  k=1;
  for (int i=1;i<=ndim;i++) k *= ng;
  if (ncall/k <= 2) npg=2;
  else npg=ncall/k;
  calls=npg*k;
  dxg=1.0/ng;

  dv2g=1;                                 // be         
  for (int i=1;i<=ndim;i++) dv2g *= dxg;  // not        
  dv2g=SQR(calls*dv2g)/npg/npg/(npg-1.0); // understood!

  xnd=nd;
  dxg *= xnd;
  xjac=1.0/calls;

  if (nd != ndo){
    for (int i=1;i<=nd;i++) r[i]=1.0;
    for (int j=1;j<=ndim;j++) rebin(ndo/xnd,nd,r,xin,xi[j]);
    ndo=nd;
  }

  for (it=1;it<=itmx;it++){
    static double ti, tsi;
    ti=tsi=0.0;
    for (int j=1;j<=ndim;j++){
      kg[j]=1;
      for (int i=1;i<=nd;i++) d[i][j]=di[i][j]=0.0;
    }

    for (;;){
      fb=f2b=0.0;
      for (k=1;k<=npg;k++){
        wgt=xjac;
        for (int j=1;j<=ndim;j++){

          xn=(kg[j]-ran2())*dxg+1.0;
          ia[j]=(int)(xn);
          if (ia[j] > NDMX) ia[j]=NDMX;
          else if (ia[j] < 1) ia[j]=1;

          xo=xi[j][ia[j]]-xi[j][ia[j]-1];
          rc=xi[j][ia[j]-1]+(xn-ia[j])*xo;

          x[j]=rc*dx[j];
          wgt *= xo*xnd;
        }

        f=wgt * fxn(x,wgt);

        f2=f*f;
        fb += f;
        f2b += f2;
        for (int j=1;j<=ndim;j++){
          di[ia[j]][j] += f;
          if (mds >= 0) d[ia[j]][j] += f2;
        }
      }

      f2b=sqrt(f2b*npg);
      f2b=(f2b-fb)*(f2b+fb);
      if (f2b <= 0.0) f2b=TINY;
      ti += fb;
      tsi += f2b;
      if (mds < 0){
        for (int j=1;j<=ndim;j++) d[ia[j]][j] += f2b;
      }
      for (k=ndim;k>=1;k--){
        kg[k] %= ng;
        if (++kg[k] != 1) break;
      }
      if (k < 1) break;
    }

    tsi *= dv2g;
    wgt=1.0/tsi;
    si += wgt*ti;
    schi += wgt*ti*ti;
    swgt += wgt;
    *tgral=si/swgt;
    *chi2a=(schi-si*(*tgral))/(it-0.9999);
    if (*chi2a<0.0) *chi2a=0.0;
    *sd=sqrt(1.0/swgt);
    tsi=sqrt(tsi);

    for (int j=1;j<=ndim;j++){
      xo=d[1][j];
      xn=d[2][j];
      d[1][j]=(xo+xn)/2.0;
      dt[j]=d[1][j];
      for (int i=2;i<nd;i++){
        rc=xo+xn;
        xo=xn;
        xn=d[i+1][j];
        d[i][j]=(rc+xn)/3.0;
        dt[j] += d[i][j];
      }
      d[nd][j]=(xo+xn)/2.0;
      dt[j] += d[nd][j];
    }

    for (int j=1;j<=ndim;j++){
      //---------------------
      if(dt[j] > 0.0){
	//---------------------
	rc=0.0;
	for (int i=1;i<=nd;i++){
	  if (d[i][j] < TINY) d[i][j]=TINY;
	  r[i]=pow((1.0-d[i][j]/dt[j])/(log(dt[j])-log(d[i][j])),ALPH);
	  rc += r[i];
	}
	rebin(rc/xnd,nd,r,xin,xi[j]);
	//---------------------
      }
      //---------------------
    }
  }
}

// Utility routine used by VEGAS                           
void rebin(double rc, int nd, double r[], double xin[], double xi[]){

  int i,k=0;
  double dr=0.0,xn=0.0,xo=0.0;

  for (i=1;i<nd;i++){
    while (rc > dr) {
      dr += r[++k];
      xo=xn;
      xn=xi[k];
    }
    dr -= rc;
    xin[i]=xn-(xn-xo)*dr/r[k];
  }
  for (i=1;i<nd;i++) xi[i]=xin[i];
  xi[nd]=1.0;
}

