/* Copyright (c) Colorado School of Mines, 2008.*/
/* All rights reserved.                       */

#include "par.h"
#include "anisotropy.h"
#define EPS_smooth FLT_MIN
#define PII 3.1415926
#define R2D 57.2957795131
/*********************** self documentation **********************/
char *sdoc[]={
"									",
" SUTITOMOSP  - ray-based gridded TI TOMOgraphy using velocity 	        ",
"               perterbation analysis                                   ",  
"									",
" sutitomosp par=cig.par nx= dx= fx= nz= dz= fz= fcdp_cig=		",
" 	 vfile= efile= dfile= nufile= ncdp= dcdp= fcdp= off0= doff=	",
"									",
"									",
" Required Parameters:							",
" ncdp=		number of cdp's		 				",
" dcdp=		horizontal interval of cdp's			 	",
" fcdp=		horizontal coordinate of the first cdp			",
" fcdp_cig=	horizontal coordinate of the first cdp in cig.par	",	
" off0=		first offset in common image gathers 			",
" doff=		offset increment in common image gathers  		",
" offmax=	maximum offset in common image gathers			",
" 									",
" vfile=	file defining VP0 at all grid points from prev. iter.	", 
" efile=	file defining eps at all grid points from prev. iter.	", 
" dfile=	file defining del at all grid points from prev. iter.	", 
" nufile=	file defining tilt at all grid points from prev. iter.	",
" 									",
" In cig.par, there should be information:				",	
" x,z,r2,r1,offmax;...							",
" description of input CIGS for n-th reflector				",
"	x	x-value of a common image point				",
"	z	z-value of a common image point				",
"	r1	hyperbolic component of the residual moveout		",
"	r2	non-hyperbolic component of residual moveout		",
"	offmax  maximun offset of this CIG used for velocity analysis	", 
" nx=100	number of samples in horizontal direction for the model	",
" nz=100	number of samples in vertical direction for the model	",
" dx=10		horizontal grid increment				",
" dz=10		vertical grid increment					",
" fx=0		first horizontal grid point				",
" fz=0		first vertical grid point				",
"									",
" Optional Parameters:							",
" dt=0.008	traveltime increment					",
" nt=500	no. of points on the ray				",
" amax=360	max. angle of emergence					",
" amin=0	min. angle of emergence					",
" Smoothing parameters:							",
" r1=0                  smoothing parameter in the 1 direction          ",
" r2=0                  smoothing parameter in the 2 direction          ",
" win=0,n1,0,n2         array for window range                          ",
" rw=0                  smoothing parameter for window function         ",
" tol=0.1		tolerance in computing the offset (m)		",
" diagnose=0		diagnose=1 with print information of ray tracing",
"									",
" There are four output files,  					",
" b_file = b_file    array for b in ||Ax-b||				",
" The following three are files for sparse storage of A in Matlab format:",
" Av_file = Av_file  array for nonzero values for A			",
" Ai_file = Ai_file  array for row number of nonzero value for A	",
" Aj_file = Aj_file  array for column number of nonzero value for A	",
"									",
" Av_min = 0.001  minimum absolute value in A which is stored		",
"									",
" Notes:								",
" The output of this program are binary files of A and b for the	",
" objective function ||Ax-b||.						",
"									",
" The symmetry axis is assumed to be perpendicular to the interface,	",
" therefore, the reflector dip field is taken the same as the tilt field,",
" which is the input of this program (nufile).				",
"									",
" This program is based on the MVA algorithm developed by 		",
" Xiaoxiang Wang, CWP:2011.						",
"									",
NULL};
/*
 * Author: Xiaoxiang Wang 
 * based on program: velpertan.c written by Debashish Sarkar.
 */
 
/**************** end self doc ***********************************/

/* functions defined and used internally */
/* one step along ray */
typedef struct RayStepStruct {
        float t;                /* time */
        float x,z;              /* x,z coordinates */
        float q1,p1,q2,p2;      /* Cerveny's dynamic ray tracing solution */
        int kmah;               /* KMAH index */
        float c,s;              /* cos(angle) and sin(angle) */
        float v,dvdx,dvdz;      /* velocity and its derivatives */
} RayStep;

/* one ray */
typedef struct RayStruct {
        int nrs;                /* number of ray steps */
        RayStep *rs;            /* array[nrs] of ray steps */
        int nc;                 /* number of circles */
        int ic;                 /* index of circle containing nearest step */
        void *c;                /* array[nc] of circles */
} Ray;

float findangle(float x,float z,float dip, float normal_ang,float amin,float amax,
		int nx,float fx,float dx,int nz,float fz,float dz,int nt,float dt,
		float **c11,float **c13,float **c15,float **c33,float **c35,float **c55,float tol);
		
		
/* functions (velocity interpolation)*/
void* vel2Alloc (int nx, float dx, float fx,
        int nz, float dz, float fz, float **v);
void vel2Free (void *vel2);
void vel2Interp (void *vel2, float x, float z,
        float *v, float *vx, float *vz, float *vxx, float *vxz, float *vzz);

Ray *makeRay (float x0, float z0, float a0, int nt, float dt,
	float **a1111xz, float **a3333xz, float **a1133xz, float **a1313xz, 
	float **a1113xz, float **a3313xz, 
	int nx, float dx, float fx, int nz, float dz, float fz, float amax, 
	float amin);
float zbrentou(float da,float da_2,float tol,float off,float x, float z,
	float a_normal,int nt,float dt,float **a1111,
	float **a3333,float **a1133, float **a1313,float **a1113,
	float **a3313,int nx,float dx,float fx,int nz,float dz,
	float fz,float amax,float amin);
float dvdlam(int para_flag,float VP0,float Vh,float Vn,float nu,float c,float s,float V);
float lagrange_scalar(int point_flag,int ix,int iz,float x,float z,float fx,float fz,float dx,float dz);
void freeRay (Ray *ray);
void smooth2 (int n1, int n2, float r1, float r2, int *win, float rw, float **v);

Stiff2D spar;  /* stiffness tensor*/

int main (int argc, char **argv)
{
	int ix,iz,ncip,icdp,icip,jcip,win[4],NRS_inc,NRS_ref,k,k0,nn,j;
	int ncdp,ioff,nx,nz,nt,is,ixs,izs,iverp1,iverp2,iverp3,iverp4;
	int nnoff,iver,noff,noff_avg,m,diag,row,*mcip,*Ai,*Aj;
	float dcdp,fcdp,fcdp_cig,doff,off0,off,maxoff;
	float fx,fz,dx,dz,amin,amax,tol,zoff_avg,sc;
	float r1,r2,rw,xs,zs,cs,ss,z2,a,scals1,scals2,scals3,scals4;
	float a_normal,a_vert=0.0,ang,X1,X2,t1,t2,T,p1,p2,p,c1,c2;
	float *xcdp,*zoff,*b,*Av,Av_min,dt,dtl;
	float **vp,**e,**d,**nu,**vh,**vn;
        float **c11,**c13,**c15,**c33,**c35,**c55;
	float *xcip,*zcip,*r1cip,*r2cip,*offmax,temp[5];
	float vs,dvdvps,dvdvhs,dvdvns;
	float *dt1dpara,*dt2dpara,**dzdpara,*dzdpara_avg;
	char *vfile,*efile,*dfile,*nufile,*b_file,*Av_file,*Ai_file,*Aj_file;
  //char filename [ FILENAME_MAX ];
  char *pkfile; FILE *pkfp;
	FILE *vf,*ef,*df,*nuf,*bfp,*Avfp,*Aifp,*Ajfp;
	
	float 	vps,dvpsdx,dvpsdz,ddvpsdxdx,ddvpsdxdz,ddvpsdzdz,
		vhs,dvhsdx,dvhsdz,ddvhsdxdx,ddvhsdxdz,ddvhsdzdz,
		vns,dvnsdx,dvnsdz,ddvnsdxdx,ddvnsdxdz,ddvnsdzdz,
		nus,dnusdx,dnusdz,ddnusdxdx,ddnusdxdz,ddnusdzdz;
	
	void *vp2;
	void *vh2;
	void *vn2;
	void *nu2;
	
	Ray *rayinc;
	Ray *rayref;
	
	/* hook up getpar to handle the parameters */
	initargs(argc,argv);
	requestdoc(0);

	/* get required parameters */
	if (!getparint("ncdp",&ncdp)) err("must specify ncdp!\n");
	if (!getparfloat("dcdp",&dcdp)) err("must specify dcdp!\n");
	if (!getparfloat("fcdp",&fcdp)) err("must specify fcdp!\n");
	if (!getparfloat("fcdp_cig",&fcdp_cig)) err("must specify fcdp_cig!\n");
	if (!getparfloat("off0",&off0)) err("must specify off0!\n");
	if (!getparfloat("doff",&doff)) err("must specify doff!\n");
	if (!getparfloat("offmax",&maxoff)) err("must specify offmax!\n");
        if (!getparstring("vfile",&vfile)) err("must specify velocity file!\n");
	if (!getparstring("efile",&efile))  err("must specify epsilon file!\n");

  /*if (!getparstring("picks",&pkfile)) err("must specify cip picks file!\n"); */

	if (!getparstring("dfile",&dfile)) err("must specify delta file!\n");
	if (!getparstring("nufile",&nufile)) err("must specify nu(tilt) file!\n");
        if (!getparint("nx",&nx)) err("must specify nx file!\n");
	if (!getparint("nz",&nz)) err("must specify nz file!\n");
	if (!getparfloat("dx",&dx)) err("must specify dx file!\n");
	if (!getparfloat("dz",&dz)) err("must specify dz file!\n");
	if (!getparfloat("fx",&fx)) err("must specify fx file!\n");
	if (!getparfloat("fz",&fz)) err("must specify fz file!\n");



	/* get optional parameters */
	if(!getparfloat("dt",&dt)) dt=0.002;
	if(!getparint("nt",&nt)) nt=2001;
	if(!getparfloat("amax",&amax)) amax=180.0;
	if(!getparfloat("amin",&amin)) amin=0.0;
	if(!getparfloat("tol",&tol)) tol=0.1;
	if (!getparint("win",win)) {
                win[0] = 0;
                win[1] = nz;
                win[2] = 0;
                win[3] = nx;
                }
	if (!getparfloat("rw",&rw)) rw = 0;
	if (!getparfloat("r1",&r1)) r1 = 0;
	if (!getparfloat("r2",&r2)) r2 = 0;
	if (!getparint("diagnose",&diag)) diag=0;
	if (!getparfloat("Av_min",&Av_min)) Av_min = 1.e-6; //0.001;
	
	if (!getparstring("Av_file",&Av_file)) Av_file = "Av_file";
	Avfp = efopen(Av_file,"w");
	
	if (!getparstring("Ai_file",&Ai_file)) Ai_file = "Ai_file";
	Aifp = efopen(Ai_file,"w");

	if (!getparstring("Aj_file",&Aj_file)) Aj_file = "Aj_file";
	Ajfp = efopen(Aj_file,"w");

	if (!getparstring("b_file",&b_file)) b_file = "b_file";
	bfp = efopen(b_file,"w");
	
	
	/* allocate space */
	vp = alloc2float(nz,nx);
	e  = alloc2float(nz,nx);
	d  = alloc2float(nz,nx);
  nu = alloc2float(nz,nx);   
  vh = alloc2float(nz,nx);
  vn = alloc2float(nz,nx);
        
  c11 = alloc2float(nz,nx);
	c13 = alloc2float(nz,nx);
  c15 = alloc2float(nz,nx);
	c33 = alloc2float(nz,nx);
	c35 = alloc2float(nz,nx);
 	c55 = alloc2float(nz,nx);
        
	xcdp = alloc1float(ncdp);
	mcip = alloc1int(ncdp);
	
	/* dt1dpara = alloc1float(3*nx*nz);
	dt2dpara = alloc1float(3*nx*nz);
	dzdpara_avg = alloc1float(3*nx*nz); */
	
	dt1dpara    = alloc1float(nx*nz);
	dt2dpara    = alloc1float(nx*nz);
	dzdpara_avg = alloc1float(nx*nz);


	
	/* input velocity, epsilon, delta and nu fields */
        if ((vf=fopen(vfile,"r"))==NULL)
		err("error opening vfile=%s!\n",vfile);
	if (fread(vp[0],sizeof(float),nz*nx,vf)!=nz*nx)
		err("error reading vfile=%s!\n",vfile);	

        if ((ef=fopen(efile,"r"))==NULL)
		err("error opening efile=%s!\n",efile);
	if (fread(e[0],sizeof(float),nz*nx,ef)!=nz*nx)
		err("error reading efile=%s!\n",efile);

        if ((df=fopen(dfile,"r"))==NULL)
		err("error opening dfile=%s!\n",dfile);
	if (fread(d[0],sizeof(float),nz*nx,df)!=nz*nx)
		err("error reading dfile=%s!\n",dfile);
		
	if ((nuf=fopen(nufile,"r"))==NULL)
		err("error opening nufile=%s!\n",nufile);
	if (fread(nu[0],sizeof(float),nz*nx,nuf)!=nz*nx)
		err("error reading nufile=%s!\n",nufile);
	
  /*if ((pkfp=fopen(pkfile,"r"))==NULL)
     err("error opening picks file=%s",pkfile); */

	smooth2(nz,nx,r1,r2,win,rw,vp);
	smooth2(nz,nx,r1,r2,win,rw,e);
	smooth2(nz,nx,r1,r2,win,rw,d);
	smooth2(nz,nx,r1,r2,win,rw,nu);	
		
	/* Compute the stiffness tensor */
  for(ix=0; ix<nx; ix++)
		for(iz=0; iz<nz; iz++)
		{
			vn[ix][iz] = vp[ix][iz]*sqrt(1+2*d[ix][iz]);
			vh[ix][iz] = vp[ix][iz]*sqrt(1+2*e[ix][iz]);
				
			thom2stiffTI(vp[ix][iz],0.5*vp[ix][iz],e[ix][iz],d[ix][iz],e[ix][iz],nu[ix][iz],&spar,1);
			
      c11[ix][iz] = (float) spar.a1111;
      c13[ix][iz] = (float) spar.a1133;
			c15[ix][iz] = (float) spar.a1113;
			c33[ix][iz] = (float) spar.a3333;
			c35[ix][iz] = (float) spar.a3313;
			c55[ix][iz] = (float) spar.a1313;
		}
	
	/* allocate and initialize velocities interpolator */
	vp2 = vel2Alloc(nx,dx,fx,nz,dz,fz,vp);
	vh2 = vel2Alloc(nx,dx,fx,nz,dz,fz,vh);
	vn2 = vel2Alloc(nx,dx,fx,nz,dz,fz,vn);
	nu2 = vel2Alloc(nx,dx,fx,nz,dz,fz,nu);
   
	/* compute uniformly sampled cdp */
   /* for(icdp=0; icdp<ncdp; ++icdp) 
     xcdp[icdp] = fcdp+icdp*dcdp; *///warn("%f", xcdp[icdp]);
     
	
	//ncip = 0;
  //warn("%s", filename);

	/* for(icdp=0; icdp<ncdp; ++icdp)
	{ 
		icip = (xcdp[icdp] - fcdp_cig)/dcdp + 1;
	
		sprintf(filename, "cip%d", icip); warn("after");
		mcip[icdp] = countparname(filename);
		mcip[icdp] = 1;
    ncip += mcip[icdp];	//warn("%d",mcip[icdp]);
	}
  warn("%d", ncip); */
		
/*	icip = (xcdp[0] - fcdp_cig)/dcdp + 1;
	sprintf(filename, "cip%d", icip);		
		
	for(jcip=0; jcip<mcip[0]; ++jcip)
 	{	
		getnparfloat(jcip+1,filename,temp);
		xcip[jcip] = temp[0];
		zcip[jcip] = temp[1];
		r2cip[jcip] = temp[2];
		r1cip[jcip] = temp[3];
		offmax[jcip] = temp[4];
	} 
	
	if(fabs(xcip[0]-xcdp[0])>0.1) 
		err("The cdp coordinates in cig.par file do not match the given parameters!\n"); 	*/

/*
  m=0;
	for(icdp=0; icdp<ncdp; ++icdp)
	{ 
		icip = (xcdp[icdp] - fcdp_cig)/dcdp + 1;
		sprintf(filename, "cip%d", icip);		
		
		if (mcip[icdp] > 0.1) {
			for(jcip=0; jcip<mcip[icdp]; ++jcip)
 			{	
				getnparfloat(jcip+1,filename,temp);
				getnparfloat(jcip+1,pkfp,temp);
				fscanf(pkfp,"%f",&temp[0]);
				fscanf(pkfp,"%f",&temp[1]);
				fscanf(pkfp,"%f",&temp[2]);
				fscanf(pkfp,"%f",&temp[3]);
				fscanf(pkfp,"%f\n",&temp[4]);

				xcip[jcip + m] = temp[0]; 
				zcip[jcip + m] = temp[1];
				r2cip[jcip + m] = temp[3];
				r1cip[jcip + m] = temp[2];
				offmax[jcip + m] = temp[4];

        warn("%f\t%f\t%f\t%f\t%f",xcip[jcip + m],zcip[jcip + m],r2cip[jcip + m],r1cip[jcip + m],offmax[jcip + m]);
			}
		 
			m += mcip[icdp];
		}
	}
*/

ncip = 20; fcdp = 1000.; dcdp = 1000.; ncdp = 4;
int nmip = 4, imip; /* points in single CIGs */

float fz0 = 500., dfz0 = 500.; /* depths of the points in single CIGs */
float *dipmax = alloc1float(ncip);
         xcip = alloc1float(ncip);
	       zcip = alloc1float(ncip);
	      r1cip = alloc1float(ncip);
	      r2cip = alloc1float(ncip);
	     offmax = alloc1float(ncip);

icip = 0;

for(icdp=0; icdp<ncdp; icdp++)
 for(imip=0; imip<nmip; imip++)
  { 
    dipmax[icip] = 50.;
    offmax[icip] = 500.;
    xcip[icip]   = fcdp + icdp*dcdp;
    zcip[icip]   = fz0  + imip*dfz0;
    r1cip[icip]  = -10.;
    r2cip[icip]  = -10.;
    icip++;
  }




float dip0 = -50., ddip = 1.;
int ndip, nndip, ndip_avg;
int  idip, dip; 
float *zdip, zdip_avg;
ncip = icip+1; 

/*	for(icip=0; icip<ncip; ++icip)
	{
		noff = (offmax[icip]-off0)/doff+1;
		row += noff;
	}
*/

row=0;
for(icip=0; icip<ncip; ++icip)
	{
		ndip = (int) (dipmax[icip]-dip0)/ddip+1;
		row += ndip;
	}

	b = alloc1float(row);
	
  /*Av = alloc1float(row*3*nx*nz);
	Ai = alloc1int(row*3*nx*nz);
	Aj = alloc1int(row*3*nx*nz); */

  Av = alloc1float(row*nx*nz);
	Ai = alloc1int(row*nx*nz);
	Aj = alloc1int(row*nx*nz);

	memset((void *) b, 0, FSIZE*row);
	
	/*noff=0;
	nnoff=0;
	k=-1;
	nn=0; */
ndip  = 0;
nndip = 0;
k     = -1;
nn    = 0;


  

	for(icip=0; icip<ncip; ++icip)
	{	
  
		/*if(xcip[icip]>(fx+(nx-1)*dx) || xcip[icip]<fx){
		warn("cdp exceeds grid dimensions \n"); break;} 
		
		noff = (offmax[icip]-off0)/doff+1;
		sc = maxoff/offmax[icip];

		warn("%d\t%f",noff,sc);

		zoff = alloc1float(noff);
		memset((void *) zoff, 0, FSIZE*noff);
		
		dzdpara = alloc2float(3*nx*nz,noff);
		memset((void *) dzdpara[0], 0, FSIZE*3*nx*nz*noff); */
    ndip = (int) (dipmax[icip]-dip0)/ddip+1;
    noff = (offmax[icip]-off0)/doff+1;
    
    zdip = alloc1float(ndip);
    zoff = alloc1float(noff);
		memset((void *) zdip, 0, FSIZE*ndip);
		
		dzdpara = alloc2float(nx*nz,ndip);
		memset((void *) dzdpara[0], 0, FSIZE*nx*nz*ndip);

		
		vel2Interp(nu2,xcip[icip],zcip[icip],&nus,&dnusdx,&dnusdz,&ddnusdxdx,&ddnusdxdz,&ddnusdzdz);
    
		if(nus>=0){
			a_normal=(PII-nus)*180/PII;
			a_vert=180-a_normal;
		} else {
			a_normal=-(PII+nus)*180/PII;
			a_vert=180+a_normal;
		}
		if (fabs(a_normal)<=90)
		warn("slope of reflector is infeasible; please pick reflector again");
		
		if (diag>0)
		printf("normal angle is %f\n",a_normal);
		
		noff_avg = 0;
    ndip_avg = 0;
			
		for(idip=0; idip<ndip; idip++) //for(ioff=0;ioff<noff;ioff++)

		{	
			off = off0+ioff*doff;

			dip = dip0+idip*ddip;
			
      z2 = pow(zcip[icip],2) + r1cip[icip]*pow(off/2,2) + 
					r2cip[icip]*pow(off/2,4)/(pow(off/2,2)+pow(zcip[icip],2));
			
      ioff = idip;		 
			if(z2>=0)	zoff[ioff] = sqrt(z2);
			else		zoff[ioff] = 0; 

      zdip[idip] = zoff[ioff];
			
			if (diag>0)	
			printf("the migrated depth of cdp %f and offset %f is %f\n",
				xcip[icip],off,zoff[ioff]);
				
			/* zoff_avg += zoff[nnoff+ioff]/noff; */ 
				
			

			ang = findangle(xcip[icip],zdip[idip],dip,a_normal,amin,amax,
					nx,fx,dx,nz,fz,dz,nt,dt,c11,c13,c15,c33,c35,c55,tol);
			
			
			if(ang<-1) break;
			else
			{ 
			
			/* nnoff += 1;
			noff_avg += 1; */
			
      nndip += 1;
			ndip_avg += 1;
			

			rayinc = makeRay(xcip[icip],zdip[idip],a_normal+ang,nt,dt,c11,c33,c13,c55,
				c15,c35,nx,dx,fx,nz,dz,fz,amax,amin);   
			
			NRS_inc=rayinc->nrs;	

			X1=-((rayinc->rs[NRS_inc-1].x-rayinc->rs[NRS_inc-2].x)/(rayinc->rs[NRS_inc-1].z-
				rayinc->rs[NRS_inc-2].z))*rayinc->rs[NRS_inc-2].z+rayinc->rs[NRS_inc-2].x;
			t1=-((rayinc->rs[NRS_inc-1].t-rayinc->rs[NRS_inc-2].t)/(rayinc->rs[NRS_inc-1].z-
				rayinc->rs[NRS_inc-2].z))*rayinc->rs[NRS_inc-2].z+rayinc->rs[NRS_inc-2].t;
			p1=1/rayinc->rs[0].v;
			if(a_normal>0.0) c1=cos((a_vert-ang)*PII/180);
			else 	 c1=cos((a_vert+ang)*PII/180);
		
			
			rayref = makeRay(xcip[icip],zdip[idip],a_normal-ang,nt,dt,c11,c33,c13,c55,
				c15,c35,nx,dx,fx,nz,dz,fz,amax,amin);   
			
			NRS_ref=rayref->nrs;

			X2=-((rayref->rs[NRS_ref-1].x-rayref->rs[NRS_ref-2].x)/(rayref->rs[NRS_ref-1].z-
				rayref->rs[NRS_ref-2].z))*rayref->rs[NRS_ref-2].z+rayref->rs[NRS_ref-2].x;
			t2=-((rayref->rs[NRS_ref-1].t-rayref->rs[NRS_ref-2].t)/(rayref->rs[NRS_ref-1].z-
				rayref->rs[NRS_ref-2].z))*rayref->rs[NRS_ref-2].z+rayref->rs[NRS_ref-2].t;
			p2=1/rayref->rs[0].v;
			if(a_normal>0.0) c2=cos((a_vert+ang)*PII/180);
			else 	 c2=cos((a_vert-ang)*PII/180);

			T=t1+t2;
			p=1/(p1*c1+p2*c2);
			
			/*if (diag>0){
			printf("est_offset=%f true_offset=%f time=%f cip=%f \n",X2-X1,off,T,xcip[icip]);
			printf("source at %f receiver at %f \n",X1,X2);
			} */
			
			memset((void *) dt1dpara, 0, FSIZE*nx*nz);
			memset((void *) dt2dpara, 0, FSIZE*nx*nz);
				
			for(is=0;is<NRS_inc-1;is++)
			{ xs = rayinc->rs[is].x;
			  zs = rayinc->rs[is].z;
			  cs = rayinc->rs[is].c;
		    ss = rayinc->rs[is].s;
			  vs = rayinc->rs[is].v;
				  
			  ixs = floor((xs-fx)/dx);
			  izs = floor((zs-fz)/dz);
				 
		    /* vel2Interp(vp2,xs,zs,&vps,&dvpsdx,&dvpsdz,&ddvpsdxdx,&ddvpsdxdz,&ddvpsdzdz);	*/
			  vel2Interp(vh2,xs,zs,&vhs,&dvhsdx,&dvhsdz,&ddvhsdxdx,&ddvhsdxdz,&ddvhsdzdz);
			  /* vel2Interp(vn2,xs,zs,&vns,&dvnsdx,&dvnsdz,&ddvnsdxdx,&ddvnsdxdz,&ddvnsdzdz); */
			  vel2Interp(nu2,xs,zs,&nus,&dnusdx,&dnusdz,&ddnusdxdx,&ddnusdxdz,&ddnusdzdz);
				  
		  	/* dvdvps = dvdlam(1,vps,vhs,vns,nus,cs,ss,vs); */
		 	  dvdvhs = dvdlam(2,vps,vhs,vns,nus,cs,ss,vs);
			  /* dvdvns = dvdlam(3,vps,vhs,vns,nus,cs,ss,vs); */
				  
			  iverp1 =  ixs*nz+izs;
			  scals1 = lagrange_scalar(1,ixs,izs,xs,zs,fx,fz,dx,dz);
				  
			  iverp2 =  (ixs+1)*nz+izs;
			  scals2 = lagrange_scalar(2,ixs,izs,xs,zs,fx,fz,dx,dz);
				  
			  iverp3 =  ixs*nz+izs+1;
			  scals3 = lagrange_scalar(3,ixs,izs,xs,zs,fx,fz,dx,dz);
				  
			  iverp4 =  (ixs+1)*nz+izs+1;
			  scals4 = lagrange_scalar(4,ixs,izs,xs,zs,fx,fz,dx,dz);
				  
			  if (is!=NRS_inc-2)
			  { 
			     /* dt1dpara[3*iverp1]   += scals1*dvdvps*dt/vs; */
			     dt1dpara[iverp1+1] += scals1*dvdvhs*dt/vs;
			     /* dt1dpara[3*iverp1+2] += scals1*dvdvns*dt/vs; */
				    
			     /* dt1dpara[3*iverp2]   += scals2*dvdvps*dt/vs; */
			     dt1dpara[iverp2+1] += scals2*dvdvhs*dt/vs;
			     /* dt1dpara[3*iverp2+2] += scals2*dvdvns*dt/vs; */
				    
		       /* dt1dpara[3*iverp3]   += scals3*dvdvps*dt/vs; */
			     dt1dpara[iverp3+1] += scals3*dvdvhs*dt/vs;
		       /* dt1dpara[3*iverp3+2] += scals3*dvdvns*dt/vs; */
				    
			     /* dt1dpara[3*iverp4]   += scals4*dvdvps*dt/vs; */
			     dt1dpara[iverp4+1] += scals4*dvdvhs*dt/vs;
			     /* dt1dpara[3*iverp4+2] += scals4*dvdvns*dt/vs; */
			  }
			  else
			  {
			     dtl = dt*rayinc->rs[NRS_inc-2].z/
			     		(rayinc->rs[NRS_inc-2].z-rayinc->rs[NRS_inc-1].z);
				    
			     /* dt1dpara[3*iverp1]   += scals1*dvdvps*dtl/vs; */
			     dt1dpara[iverp1+1] += scals1*dvdvhs*dtl/vs;
			     /* dt1dpara[3*iverp1+2] += scals1*dvdvns*dtl/vs; */
				    
			     /* dt1dpara[3*iverp2]   += scals2*dvdvps*dtl/vs; */ 
			     dt1dpara[iverp2+1] += scals2*dvdvhs*dtl/vs;
			     /* dt1dpara[3*iverp2+2] += scals2*dvdvns*dtl/vs; */
				    
			    /*  dt1dpara[3*iverp3]   += scals3*dvdvps*dtl/vs; */
			     dt1dpara[iverp3+1] += scals3*dvdvhs*dtl/vs;
			    /* dt1dpara[3*iverp3+2] += scals3*dvdvns*dtl/vs; */
				    
			    /* dt1dpara[3*iverp4]   += scals4*dvdvps*dtl/vs; */
			     dt1dpara[iverp4+1] += scals4*dvdvhs*dtl/vs;
			    /* dt1dpara[3*iverp4+2] += scals4*dvdvns*dtl/vs; */ 
			  }
			}  
				
				
			for(is=0; is<NRS_ref-1; is++)
			{ xs = rayref->rs[is].x;
			  zs = rayref->rs[is].z;
			  cs = rayref->rs[is].c;
			  ss = rayref->rs[is].s;
			  vs = rayref->rs[is].v;
			  
			  ixs = floor((xs-fx)/dx);
			  izs = floor((zs-fz)/dz);
				 
			 /* vel2Interp(vp2,xs,zs,&vps,&dvpsdx,&dvpsdz,&ddvpsdxdx,&ddvpsdxdz,&ddvpsdzdz);*/	
			  vel2Interp(vh2,xs,zs,&vhs,&dvhsdx,&dvhsdz,&ddvhsdxdx,&ddvhsdxdz,&ddvhsdzdz);
			 /* vel2Interp(vn2,xs,zs,&vns,&dvnsdx,&dvnsdz,&ddvnsdxdx,&ddvnsdxdz,&ddvnsdzdz);*/
			  vel2Interp(nu2,xs,zs,&nus,&dnusdx,&dnusdz,&ddnusdxdx,&ddnusdxdz,&ddnusdzdz);
				  
			 /* dvdvps = dvdlam(1,vps,vhs,vns,nus,cs,ss,vs);*/
			  dvdvhs = dvdlam(2,vps,vhs,vns,nus,cs,ss,vs);
			/*  dvdvns = dvdlam(3,vps,vhs,vns,nus,cs,ss,vs);*/
				  
			  iverp1 =  ixs*nz+izs;
			  scals1 = lagrange_scalar(1,ixs,izs,xs,zs,fx,fz,dx,dz);
			  
			  iverp2 =  (ixs+1)*nz+izs;
			  scals2 = lagrange_scalar(2,ixs,izs,xs,zs,fx,fz,dx,dz);
			  
			  iverp3 =  ixs*nz+izs+1;
			  scals3 = lagrange_scalar(3,ixs,izs,xs,zs,fx,fz,dx,dz);
				  
			  iverp4 =  (ixs+1)*nz+izs+1;
			  scals4 = lagrange_scalar(4,ixs,izs,xs,zs,fx,fz,dx,dz);
				  
			  if (is!=NRS_ref-2)
			  { 
			   /* dt2dpara[3*iverp1]   += scals1*dvdvps*dt/vs;*/
			    dt2dpara[iverp1+1] += scals1*dvdvhs*dt/vs;
			   /* dt2dpara[3*iverp1+2] += scals1*dvdvns*dt/vs;*/
			    
			   /* dt2dpara[3*iverp2]   += scals2*dvdvps*dt/vs;*/
			    dt2dpara[iverp2+1] += scals2*dvdvhs*dt/vs;
			  /*  dt2dpara[3*iverp2+2] += scals2*dvdvns*dt/vs;*/
			    
			   /* dt2dpara[3*iverp3]   += scals3*dvdvps*dt/vs;*/
			    dt2dpara[iverp3+1] += scals3*dvdvhs*dt/vs;
			  /*  dt2dpara[3*iverp3+2] += scals3*dvdvns*dt/vs;*/
				    
			   /* dt2dpara[3*iverp4]   += scals4*dvdvps*dt/vs;*/
			    dt2dpara[iverp4+1] += scals4*dvdvhs*dt/vs;
			   /* dt2dpara[3*iverp4+2] += scals4*dvdvns*dt/vs;*/ 
			  }
			  else
			  {
			    dtl = dt*rayinc->rs[NRS_ref-2].z/
				    		(rayinc->rs[NRS_ref-2].z-rayinc->rs[NRS_ref-1].z);
				    
			  /*  dt2dpara[3*iverp1]   += scals1*dvdvps*dtl/vs;*/
			    dt2dpara[iverp1+1] += scals1*dvdvhs*dtl/vs;
			   /* dt2dpara[3*iverp1+2] += scals1*dvdvns*dtl/vs;*/
				    
			   /* dt2dpara[3*iverp2]   += scals2*dvdvps*dtl/vs;*/
			    dt2dpara[iverp2+1] += scals2*dvdvhs*dtl/vs;
			   /* dt2dpara[3*iverp2+2] += scals2*dvdvns*dtl/vs;*/
				    
			   /* dt2dpara[3*iverp3]   += scals3*dvdvps*dtl/vs;*/
			    dt2dpara[iverp3+1] += scals3*dvdvhs*dtl/vs;
			  /*  dt2dpara[3*iverp3+2] += scals3*dvdvns*dtl/vs;*/
			    
			  /*  dt2dpara[3*iverp4]   += scals4*dvdvps*dtl/vs;*/
			    dt2dpara[iverp4+1] += scals4*dvdvhs*dtl/vs;
			  /*  dt2dpara[3*iverp4+2] += scals4*dvdvns*dtl/vs;*/  
			  }
			}				
				
			for (iver=0;iver<nx*nz;iver++)
			{
				dzdpara[idip][iver] = p*(dt1dpara[iver] + dt2dpara[iver]);	
			} 
				
				
			freeRay(rayinc);
			freeRay(rayref);

			}		
		}
			
		memset((void *) dzdpara_avg, 0, FSIZE*nx*nz);
		
		zoff_avg = 0;
		zdip_avg = 0;

		for(idip=0;idip<ndip_avg;idip++)
		{
			zdip_avg += zdip[idip]/ndip_avg;
			
			for (iver=0;iver<nx*nz;iver++)
			{
				dzdpara_avg[iver] += dzdpara[idip][iver]/ndip_avg; warn("%f",dzdpara_avg[iver]);
			}
		}
			
		
		for (idip=0;idip<ndip_avg;idip++)
		 {	
      k0=k;
		  for (iver=0;iver<nx*nz;iver++)
			 {
				/*A[iver][ioff+nnoff-noff_avg] = dzdpara[ioff][iver] - dzdpara_avg[iver];*/

				 a = sc*(dzdpara[idip][iver] - dzdpara_avg[iver]);
        					
				if (fabs(a)>Av_min)
				 {	
          k += 1;
					Av[k] = a; warn("%f",a);
					Aj[k] = iver+1;			
				 }
			 }
			
			if (k==k0)
			   nn += 1;
			else 
      {
				b[idip+nndip-ndip_avg-nn] = sc*(zdip[idip] - zdip_avg); warn("%f",b[idip+nndip-ndip_avg-nn]);
				for (j=k0+1;j<k+1;j++) Ai[j] = idip+nndip-ndip_avg-nn+1;			
			}
		}
			
		free1float(zoff); free1float(zdip);
    free2float(dzdpara);
        
	
	}

	if (diag>0){
        printf("%d offsets are used to construct the inverse matrix.\n",nnoff);
  	}
	
	
	/*efwrite(*A,sizeof(float),nnoff*3*nx*nz,Afp);*/
	
	efwrite(b,sizeof(float),nndip,bfp);
	efwrite(Av,sizeof(float),k+1,Avfp);
	efwrite(Aj,sizeof(int),k+1,Ajfp);
	efwrite(Ai,sizeof(int),k+1,Aifp);
	

        free2float(c11);
        free2float(c13);
        free2float(c15);
        free2float(c33);
        free2float(c35);
        free2float(c55);
        
        free2float(vp);
        free2float(e);
        free2float(d);
        free2float(nu);
        free2float(vh);
        free2float(vn);

	      free1float(xcip);
      	free1float(zcip);
        free1float(r1cip);
        free1float(r2cip);
	      free1float(offmax); free1float(dipmax);
	      free1float(xcdp);
	      free1int(mcip);
        
        free1float(dt1dpara);
        free1float(dt2dpara);
        free1float(dzdpara_avg);
        
        free1float(b);
  	    free1float(Av);
	      free1int(Aj);
	      free1int(Ai);
        
	      vel2Free(vp2);
	      vel2Free(vh2);
	      vel2Free(vn2);
	      vel2Free(nu2);
	
return EXIT_SUCCESS;
}





float findangle(float x,float z,float dip, float normal_ang,float amin,float amax,
		int nx,float fx,float dx,int nz,float fz,float dz,int nt,float dt,
		float **c11,float **c13,float **c15,float **c33,float **c35,float **c55,float tol)
		
/*****************************************************************************
***************** ray tracing for a given offset and cdp point ***************
******************************************************************************	
Input:
x		x coordinate of takeoff point
z		z coordinate of takeoff point
off		offset for the specific event
normal_ang	angle between normal vector and z (clockwise rotation is negative)
amin            minimum emergence angle	
amax            maximum emergence angle
tol 		tolerance in searching for the specific event	
*****************************************************************************/
{	
	int NRS;
	//float da,X,da_2,angle,X1,X2;
	float da, angle;
  float psx, psz, prx, prz, px, pz, dip_c, mag;
	Ray *ray;
	
	/************ find the time for a given offset and diffractor point **********/
	/* find the closest two rays that bracket the desired offset */
	do{
	da=da+1;
	ray = makeRay(x,z,normal_ang+da,nt,dt,c11,c33,c13,c55,
			c15,c35,nx,dx,fx,nz,dz,fz,amax,amin);          
	NRS=ray->nrs;	
	if(ray->rs[NRS-1].z>0) {angle=-2; break;}
  psx = ray->rs[0].p1;
  psz = ray->rs[0].p2;
	/*X1=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x; */
	freeRay(ray);
	ray = makeRay(x,z,normal_ang-da,nt,dt,c11,c33,c13,c55,
				c15,c35,nx,dx,fx,nz,dz,fz,amax,amin);   
	NRS=ray->nrs;
	if(ray->rs[NRS-1].z>0) {angle=-2; break;}	
  prx = ray->rs[0].p1;
  prz = ray->rs[0].p2;
	/*X2=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x; */
	freeRay(ray);
	/* X=X2-X1; */ 
  px = psx+prx; pz = psz+prz; 
  mag = sqrt(px*px+pz*pz);
  dip_c = asin(px/(mag+1e-10))*R2D;

	if (fabs(dip_c-dip) < dip*0.02) break;
	}while(dip_c<dip);
	//da_2=da-1;

/*	X=0; angle=0;da=-1;
  do{
	da=da+1;
	ray = makeRay(x,z,normal_ang+da,nt,dt,c11,c33,c13,c55,
			c15,c35,nx,dx,fx,nz,dz,fz,amax,amin);          
	NRS=ray->nrs;	
	if(ray->rs[NRS-1].z>0) {angle=-2; break;}
	X1=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x;
	freeRay(ray);
	ray = makeRay(x,z,normal_ang-da,nt,dt,c11,c33,c13,c55,
				c15,c35,nx,dx,fx,nz,dz,fz,amax,amin);   
	NRS=ray->nrs;
	if(ray->rs[NRS-1].z>0) {angle=-2; break;}	
	X2=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x;
	freeRay(ray);
	X=X2-X1;
	if (fabs(X-off)< off*0.002) break;
	}while(X<off);
	da_2=da-1;

	if(angle>=0)
	{
		if(fabs(X-off)>=off*0.002) 
		angle = zbrentou(da,da_2,tol,off,x,z,normal_ang,nt,dt,c11,c33,c13,c55,c15,c35,nx,dx,fx,nz,dz,fz,amax,amin);
		else
		angle = da;
	} */
	
	return da;
}




/**************************************************************************/
/****** the numerical algorithm to calculate the right takeoff angle ******/
/**************************************************************************/
float zbrentou(float da,float da_2,float tol,float off,float x, float z,
		float a_normal,int nt,float dt,float **a1111,
		float **a3333,float **a1133, float **a1313,float **a1113,
		float **a3313,int nx,float dx,float fx,int nz,float dz,
		float fz,float amax,float amin)
{
 int NRS,iter,ITMAX;
 float a=da,b=da_2,c=da_2,d=0,e=0,min1,min2,X1,X2;
 float fa,fb,fc,p,q,r,s,tol1,xm,EPS;
 Ray *ray;

ITMAX=1000; EPS=3e-8;

ray = makeRay(x,z,a_normal+a,nt,dt,a1111,a3333,a1133,a1313,
                        a1113,a3313,nx,dx,fx,nz,dz,fz,amax,amin);
NRS=ray->nrs;
	X1=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x;
freeRay(ray);
ray = makeRay(x,z,a_normal-a,nt,dt,a1111,a3333,a1133,a1313,
                        a1113,a3313,nx,dx,fx,nz,dz,fz,amax,amin);
NRS=ray->nrs;
	X2=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x;
freeRay(ray);
fa=X2-X1-off;
/*printf("%f \t",X2-X1);*/
ray = makeRay(x,z,a_normal+b,nt,dt,a1111,a3333,a1133,a1313,
                        a1113,a3313,nx,dx,fx,nz,dz,fz,amax,amin);
NRS=ray->nrs;
	X1=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x;
freeRay(ray);
ray = makeRay(x,z,a_normal-b,nt,dt,a1111,a3333,a1133,a1313,
                        a1113,a3313,nx,dx,fx,nz,dz,fz,amax,amin);
NRS=ray->nrs;
	X2=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x;
freeRay(ray);
fb=X2-X1-off;
/*printf("%f \n",X2-X1);*/
if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) 
    printf("Root must be bracketed in zbrent\n");

fc=fb;
for (iter = 1;iter<=ITMAX;iter++){
   if ((fb > 0.0 && fc > 0.0)||(fb < 0.0 && fc < 0.0)) {
       c=a;
       fc=fa;
       e=d=b-a;
}
if (fabs(fc) < fabs(fb)) {
       a=b;
       b=c;
       c=a;
       fa=fb;
       fb=fc;
       fc=fa;
}

tol1=2.0*EPS*fabs(b)+0.5*tol;
xm=0.5*(c-b);
if ((fabs(xm)<=tol1)||(fb==0) ){
return b;
}
if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
       s=fb/fa;
       if (a==c) {
          p=2.0*xm*s;
          q=1.0-s;
       } else {
          q=fa/fc;
          r=fb/fc;
          p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1));
          q=(q-1.0)*(r-1.0)*(s-1.0);
}
if (p>0.0) q = -q;
p=fabs(p);
min1=3.0*xm*q-fabs(tol1*q);
min2=fabs(e*q);
if(2.0*p < (min1-min2 ? min1 : min2)){
     e=d;
     d=p/q;
} else {
     d=xm;
     e=d;
}
}else{
   d=xm;
   e=d;
}
a=b;
fa=fb;
if (fabs(d) > tol1)
    b+=d;
else
    b += fabs(tol1)*xm/fabs(xm);
ray = makeRay(x,z,a_normal+b,nt,dt,a1111,a3333,a1133,a1313,
                        a1113,a3313,nx,dx,fx,nz,dz,fz,amax,amin);
NRS=ray->nrs;
	X1=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x;
freeRay(ray);
ray = makeRay(x,z,a_normal-b,nt,dt,a1111,a3333,a1133,a1313,
                        a1113,a3313,nx,dx,fx,nz,dz,fz,amax,amin);
NRS=ray->nrs;
	X2=-((ray->rs[NRS-1].x-ray->rs[NRS-2].x)/(ray->rs[NRS-1].z-
		ray->rs[NRS-2].z))*ray->rs[NRS-2].z+ray->rs[NRS-2].x;
freeRay(ray);
fb=X2-X1-off;
}
printf("Maximum number of iterations exceeded in zbrent\n");
return 0.0;
}



/****************************************************************************/
/***************************  Tariq's ray tracer ****************************/
/****************************************************************************/
Ray *makeRay (float x0, float z0, float a0, int nt, float dt,
	float **a1111xz, float **a3333xz, float **a1133xz, float **a1313xz, 
	float **a1113xz, float **a3313xz, 
	int nx, float dx, float fx, int nz, float dz, float fz, float amax, 
	float amin)
/*****************************************************************************
Trace a ray for uniformly sampled v(x,z).
******************************************************************************
Input:
x0		x coordinate of takeoff point
z0		z coordinate of takeoff point
a0		takeoff angle (radians)
nt		number of time samples
dt		time sampling interval
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
amax            maximum emergence angle
amin            minimum emergence angle
a1111		array[nx][nz] of uniformly sampled density normalized elastic coef.
a3333		array[nx][nz] of uniformly sampled density normalized elastic coef.
a1133           array[nx][nz] of uniformly sampled density normalized elastic coef.
a1313           array[nx][nz] of uniformly sampled density normalized elastic coef.
a1113           array[nx][nz] of uniformly sampled density normalized elastic coef.
a3313           array[nx][nz] of uniformly sampled density normalized elastic coef.
******************************************************************************
Returned:	pointer to ray parameters sampled at discrete ray steps
******************************************************************************
Notes:
The ray ends when it runs out of time (after nt steps) or with the first 
step that is out of the (x,z) bounds of the velocity function v(x,z).
*****************************************************************************
Technical Reference:

Cerveny, V., 1972, Seismic rays and ray intensities 
	in inhomogeneous anisotropic media: 
	Geophys. J. R. Astr. Soc., 29, 1--13.

*****************************************************************************
 Credits: CWP: Tariq Alkhalifah
*****************************************************************************/
{
	int it,kmah;
	float t,x,z,c,s,p1,q1,p2,q2,px,pz,px2,pz2,pxz,
		lx,lz,cc,ss,
		a1111,da1111dx,da1111dz,dda1111dxdx,dda1111dxdz,dda1111dzdz,
		a3333,da3333dx,da3333dz,dda3333dxdx,dda3333dxdz,dda3333dzdz,
		a1133,da1133dx,da1133dz,dda1133dxdx,dda1133dxdz,dda1133dzdz,
		a1313,da1313dx,da1313dz,dda1313dxdx,dda1313dxdz,dda1313dzdz,
		a1113,da1113dx,da1113dz,dda1113dxdx,dda1113dxdz,dda1113dzdz,
		a3313,da3313dx,da3313dz,dda3313dxdx,dda3313dxdz,dda3313dzdz,
		da1111dn,dda1111dndn,da3333dn,dda3333dndn,da1133dn,dda1133dndn,
		da1313dn,dda1313dndn,da1113dn,dda1113dndn,da3313dn,dda3313dndn,
		gamm11,gamm13,gamm33,vp2,vp,ovp,sqr;
	Ray *ray;
	RayStep *rs;
	void *a11112;
	void *a33332;
	void *a11332;
	void *a13132;
	void *a11132;
	void *a33132;

	/*Convert degrees to radians*/
	a0=a0*PII/180; amax=amax*PII/180; amin=amin*PII/180;
	/* allocate and initialize velocities interpolator */
	a11112 = vel2Alloc(nx,dx,fx,nz,dz,fz,a1111xz);
	a33332 = vel2Alloc(nx,dx,fx,nz,dz,fz,a3333xz);
	a11332 = vel2Alloc(nx,dx,fx,nz,dz,fz,a1133xz);
	a13132 = vel2Alloc(nx,dx,fx,nz,dz,fz,a1313xz);
	a11132 = vel2Alloc(nx,dx,fx,nz,dz,fz,a1113xz);
	a33132 = vel2Alloc(nx,dx,fx,nz,dz,fz,a3313xz);
	
	/* last x and z in velocity model */
	lx = fx+(nx-1)*dx;
	lz = fz+(nz-1)*dz;

	/* ensure takeoff point is within model */
	if (x0<fx || x0>lx || z0<fz || z0>lz ){
	warn("The CRP lies outside the defined grid");
	return NULL;
	}

	/* allocate space for ray and raysteps */
	ray = (Ray*)alloc1(1,sizeof(Ray));
	rs = (RayStep*)alloc1(nt,sizeof(RayStep));

	/* cosine and sine of takeoff angle */
	c = cos(a0);
	s = sin(a0);
	cc = c*c;
	ss = s*s;
	
	/* velocities and derivatives at takeoff point */
	vel2Interp(a11112,x0,z0,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
		&dda1111dxdz,&dda1111dzdz);
	da1111dn    = da1111dx*c-da1111dz*s;
	dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

	vel2Interp(a33332,x0,z0,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
	da3333dn    = da3333dx*c-da3333dz*s;
	dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
	vel2Interp(a11332,x0,z0,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
		&dda1133dxdz,&dda1133dzdz);
	da1133dn    = da1133dx*c-da1133dz*s;
	dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

	vel2Interp(a13132,x0,z0,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
		&dda1313dxdz,&dda1313dzdz);
	da1313dn    = da1313dx*c-da1313dz*s;
	dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

	vel2Interp(a11132,x0,z0,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
		&dda1113dxdz,&dda1113dzdz);
	da1113dn    = da1113dx*c-da1113dz*s;
	dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

	vel2Interp(a33132,x0,z0,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
		&dda3313dxdz,&dda3313dzdz);
	da3313dn    = da3313dx*c-da3313dz*s;
	dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;

	/*computing the phase velocity for a0 angle */
	gamm11 = a1111*ss+ a1313*cc +2*a1113*s*c;
	gamm33 = a3333*cc + a1313*ss+2*a3313*s*c;
	gamm13 = (a1133+a1313)*s*c+ a1113*ss+ a3313*cc;
	sqr    = sqrt((gamm11+gamm33)*(gamm11+gamm33)-
			4*(gamm11*gamm33-gamm13*gamm13));
	vp2    = gamm11+gamm33+sqr;
	vp     = sqrt(vp2*.5);
	ovp    = 1/vp;
	px     = s*ovp;
	pz     = c*ovp;

	/* first ray step */
	rs[0].t = t = 0;
	rs[0].x = x = x0;
	rs[0].z = z = z0;
	rs[0].q1 = q1 = 1.0;
	rs[0].p1 = p1 = 0.0;
	rs[0].q2 = q2 = 0.0;
	rs[0].p2 = p2 = 1.0;
	rs[0].kmah = kmah = 0;
	rs[0].c = c;
	rs[0].s = s;
	rs[0].v = vp;
	rs[0].dvdx = .5*da3333dx*vp/a3333;
	rs[0].dvdz = .5*da3333dz*vp/a3333;

	/* loop over time steps */
	for (it=1; it<nt; ++it) {

		/* variables used for Runge-Kutta integration */
		float h=dt,hhalf=dt/2.0,hsixth=dt/6.0,
			q2old,xt,zt,p1t,q1t,p2t,q2t,
			dx,dz,dp1,dq1,dp2,dq2,
			dxt,dzt,dp1t,dq1t,dp2t,dq2t,
			dxm,dzm,dp1m,dq1m,dp2m,dq2m,
			gamma11,gamma33,gamma13,g11,g13,g33,den,
			sxfactor,szfactor,snfact,dpx,dpz,pxt,pzt,dpxt,
			dpzt,dpxm,dpzm,dxdn,dzdn,snfactor,
			dxx,dzz,dcdp1,dcdp3,dcdp13,ddcdnn,ddcdqq,
			ddcdpn,dgamma11dpx,dgamma11dpz,dgamma33dpx,
			dgamma33dpz,dgamma13dpx,dgamma13dpz,dg11dpx,
			dg11dpz,dg33dpx,dg33dpz,dg13dpx,dg13dpz,ddxdpx,
			ddzdpz,ddxdpz,dgamma11dn,dgamma33dn,dgamma13dn,
			dg11dn,dg33dn,dg13dn,dsnfactordn,ddxdn,ddzdn;

		/* if ray is out of bounds, break */
		if (x<fx || x>lx || z<fz || z>lz || c>(cos(amin)+0.01) || c<(cos(amax))-0.01) break;

		/* remember old q2 */
		q2old = q2;
		
	        /* step 1 of 4th-order Runge-Kutta */
		px2   = px*px;
		pz2   = pz*pz;
		pxz   = px*pz;

		/*anisotropy parameters*/
		gamma11 = a1111*px2+ a1313*pz2 +2*a1113*pxz;
		gamma33 = a3333*pz2 + a1313*px2+2*a3313*pxz;
		gamma13 = (a1133+a1313)*pxz+ a1113*px2+ a3313*pz2;
		den     = 1/(gamma11+gamma33-2);
		g11     = (gamma33-1)*den;
		g33     = (gamma11-1)*den;
		g13     = -gamma13*den;
		sxfactor = da1111dx*px2*g11+da3333dx*pz2*g33+
			2*(da1133dx+da1313dx)*pxz*g13+da1313dx*(px2*g33+pz2*g11)+
			2*da3313dx*(pz2*g13+pxz*g33)+2*da1113dx*(pxz*g11+px2*g13);
		szfactor = da1111dz*px2*g11+da3333dz*pz2*g33+
			2*(da1133dz+da1313dz)*pxz*g13+da1313dz*(px2*g33+pz2*g11)+
			2*da3313dz*(pz2*g13+pxz*g33)+2*da1113dz*(pxz*g11+px2*g13);
		snfact = sxfactor*c-szfactor*s;
		
		/*computing ray velocities and derivatives*/
		dx =  (a1111*px*g11+(a1133+a1313)*pz*g13+a3313*pz*g33
			+a1113*(pz*g11+2*px*g13)+a1313*g33*px);
		dz =  (a3333*pz*g33+(a1133+a1313)*px*g13+a1113*px*g11
			+a3313*(px*g33+2*pz*g13)+a1313*g11*pz);

		dgamma11dpx = 2*a1111*px+2*a1113*pz;
		dgamma11dpz = 2*a1313*pz+2*a1113*px;
		dgamma33dpx = 2*a1313*px+2*a3313*pz;
		dgamma33dpz = 2*a3333*pz+2*a3313*px;
		dgamma13dpx= (a1133+a1313)*pz+2*a1113*px;
		dgamma13dpz= (a1133+a1313)*px+2*a3313*pz;
		dgamma11dn = da1111dn*px2+ da1313dn*pz2 +2*da1113dn*pxz;
		dgamma33dn = da3333dn*pz2 + da1313dn*px2+2*da3313dn*pxz;
		dgamma13dn = (da1133dn+da1313dn)*pxz+ da1113dn*px2+ da3313dn*pz2;
		dg11dpx = -(gamma33-1)*(dgamma11dpx+dgamma33dpx-4*dx)*den*den+
			(dgamma33dpx-2*dx)*den;
		dg11dpz = -(gamma33-1)*(dgamma11dpz+dgamma33dpz-4*dz)*den*den+
			(dgamma33dpz-2*dz)*den;
		dg33dpx = -(gamma11-1)*(dgamma11dpx+dgamma33dpx-4*dx)*den*den+
			(dgamma11dpx-2*dx)*den;
		dg33dpz = -(gamma11-1)*(dgamma11dpz+dgamma33dpz-4*dz)*den*den+
			(dgamma11dpz-2*dz)*den;
		dg13dpx = gamma13*(dgamma11dpx+dgamma33dpx-4*dx)*den*den-
			dgamma13dpx*den;
		dg13dpz = gamma13*(dgamma11dpz+dgamma33dpz-4*dz)*den*den-
			dgamma13dpz*den;
		dg11dn = -(gamma33-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma33dn-snfact)*den;
		dg33dn = -(gamma11-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma11dn-snfact)*den;
		dg13dn = gamma13*(dgamma11dn+dgamma33dn-2*snfact)*den*den-
			dgamma13dn*den;
		ddxdpx=   a1111*px*dg11dpx+(a1133+a1313)*pz*dg13dpx+
			a3313*pz*dg33dpx+a1113*(pz*dg11dpx+2*px*dg13dpx)
			+a1313*dg33dpx*px;
		ddzdpz= a3333*pz*dg33dpz+(a1133+a1313)*px*dg13dpz+
			a1113*px*dg11dpz+a3313*(px*dg33dpz+2*pz*dg13dpz)+
			a1313*dg11dpz*pz;
		ddxdpz= a1111*px*dg11dpz+(a1133+a1313)*pz*dg13dpz+
			a3313*pz*dg33dpz+a1113*(pz*dg11dpz+2*px*dg13dpz)+
			a1313*dg33dpz*px;
		dsnfactordn = da1111dn*px2*dg11dn+da3333dn*pz2*dg33dn+
			2*(da1133dn+da1313dn)*pxz*dg13dn+da1313dn*(px2*dg33dn+pz2*dg11dn)+
			2*da3313dn*(pz2*dg13dn+pxz*dg33dn)+2*da1113dn*(pxz*dg11dn+px2*dg13dn);
		ddxdn =  (a1111*px*dg11dn+(a1133+a1313)*pz*dg13dn+a3313*pz*dg33dn
			+a1113*(pz*dg11dn+2*px*dg13dn)+a1313*dg33dn*px);
		ddzdn =  (a3333*pz*dg33dn+(a1133+a1313)*px*dg13dn+a1113*px*dg11dn
			+a3313*(px*dg33dn+2*pz*dg13dn)+a1313*dg11dn*pz);
		

		/*evaluating change in slowness and amplitude along ray*/
		dpx = -0.5*sxfactor;
		dpz = -0.5*szfactor;

		dcdp1  = a1111*g11+a1313*g33+2*a1113*g13+ddxdpx-dx*dx;
		dcdp3  = a3333*g33+a1313*g11+2*a3313*g13+ddzdpz-dz*dz;
		dcdp13 = a1133*g13+a1313*g13+a1113*g11+a3313*g33+ddxdpz-dx*dz;
		ddcdqq = dcdp1*cc-2.0*dcdp13*s*c+dcdp3*ss;
		dxdn   =  (da1111dn*px*g11+(da1133dn+da1313dn)*pz*g13+da3313dn*pz*g33
			+da1113dn*(pz*g11+2*px*g13)+da1313dn*g33*px);
		dzdn   =  (da3333dn*pz*g33+(da1133dn+da1313dn)*px*g13+da1113dn*px*g11
			+da3313dn*(px*g33+2*pz*g13)+da1313dn*g11*pz);
		ddcdpn = dxdn*c-dzdn*s-.5*dx*sxfactor*cc+
			.5*(dx*szfactor+dz*sxfactor)*s*c-.5*dz*szfactor*ss
			+ddxdn*c-ddzdn*s;
		snfactor = dda1111dndn*px2*g11+dda3333dndn*pz2*g33+
			2*(dda1133dndn+dda1313dndn)*pxz*g13+
			dda1313dndn*(px2*g33+pz2*g11)+
			2*dda3313dndn*(pz2*g13+pxz*g33)+
			2*dda1113dndn*(pxz*g11+px2*g13);
		ddcdnn = 0.5*snfactor-.25*sxfactor*sxfactor*cc+
			.5*sxfactor*szfactor*s*c-.25*szfactor*szfactor*ss
			+.5*dsnfactordn;

		dp1 = -ddcdnn*q1-ddcdpn*p1;
		dq1 = ddcdqq*p1+ddcdpn*q1;
		dp2 = -ddcdnn*q2-ddcdpn*p2;
		dq2 = ddcdqq*p2+ddcdpn*q2;
		xt = x+hhalf*dx;
		zt = z+hhalf*dz;
		pxt = px+hhalf*dpx;
		pzt = pz+hhalf*dpz;
		p1t = p1+hhalf*dp1;
		q1t = q1+hhalf*dq1;
		p2t = p2+hhalf*dp2;
		q2t = q2+hhalf*dq2;
		vp  = 1/sqrt(pxt*pxt+pzt*pzt);
		s   = pxt*vp;
		c   = pzt*vp;
		ss  = s*s;
		cc  = c*c;
		
		vel2Interp(a11112,xt,zt,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
			&dda1111dxdz,&dda1111dzdz);
		da1111dn    = da1111dx*c-da1111dz*s;
		dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

		vel2Interp(a33332,xt,zt,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
		da3333dn    = da3333dx*c-da3333dz*s;
		dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
		vel2Interp(a11332,xt,zt,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
			&dda1133dxdz,&dda1133dzdz);
		da1133dn    = da1133dx*c-da1133dz*s;
		dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

		vel2Interp(a13132,xt,zt,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
			&dda1313dxdz,&dda1313dzdz);
		da1313dn    = da1313dx*c-da1313dz*s;
		dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

		vel2Interp(a11132,xt,zt,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
			&dda1113dxdz,&dda1113dzdz);
		da1113dn    = da1113dx*c-da1113dz*s;
		dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

		vel2Interp(a33132,xt,zt,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
			&dda3313dxdz,&dda3313dzdz);
		da3313dn    = da3313dx*c-da3313dz*s;
		dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;
		
        	/* step 2 of 4th-order Runge-Kutta */
		px2   = pxt*pxt;
		pz2   = pzt*pzt;
		pxz   = pxt*pzt;

		/*anisotropy parameters*/
		gamma11 = a1111*px2+ a1313*pz2 +2*a1113*pxz;
		gamma33 = a3333*pz2 + a1313*px2+2*a3313*pxz;
		gamma13 = (a1133+a1313)*pxz+ a1113*px2+ a3313*pz2;
		den     = 1/(gamma11+gamma33-2);
		g11     = (gamma33-1)*den;
		g33     = (gamma11-1)*den;
		g13     = -gamma13*den;
		sxfactor = da1111dx*px2*g11+da3333dx*pz2*g33+
			2*(da1133dx+da1313dx)*pxz*g13+da1313dx*(px2*g33+pz2*g11)+
			2*da3313dx*(pz2*g13+pxz*g33)+2*da1113dx*(pxz*g11+px2*g13);
		szfactor = da1111dz*px2*g11+da3333dz*pz2*g33+
			2*(da1133dz+da1313dz)*pxz*g13+da1313dz*(px2*g33+pz2*g11)+
			2*da3313dz*(pz2*g13+pxz*g33)+2*da1113dz*(pxz*g11+px2*g13);
		snfact = sxfactor*c-szfactor*s;
		
		/*computing ray velocities and derivatives*/
		dxt =  (a1111*pxt*g11+(a1133+a1313)*pzt*g13+a3313*pzt*g33
			+a1113*(pzt*g11+2*pxt*g13)+a1313*g33*pxt);
		dzt =  (a3333*pzt*g33+(a1133+a1313)*pxt*g13+a1113*pxt*g11
			+a3313*(pxt*g33+2*pzt*g13)+a1313*g11*pzt);
		dpxt = -0.5*sxfactor;
		dpzt = -0.5*szfactor;

		dgamma11dpx = 2*a1111*pxt+2*a1113*pzt;
		dgamma11dpz = 2*a1313*pzt+2*a1113*pxt;
		dgamma33dpx = 2*a1313*pxt+2*a3313*pzt;
		dgamma33dpz = 2*a3333*pzt+2*a3313*pxt;
		dgamma13dpx= (a1133+a1313)*pzt+2*a1113*pxt;
		dgamma13dpz= (a1133+a1313)*pxt+2*a3313*pzt;
		dgamma11dn = da1111dn*px2+ da1313dn*pz2 +2*da1113dn*pxz;
		dgamma33dn = da3333dn*pz2 + da1313dn*px2+2*da3313dn*pxz;
		dgamma13dn = (da1133dn+da1313dn)*pxz+ da1113dn*px2+ da3313dn*pz2;
		dg11dpx = -(gamma33-1)*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den+
			(dgamma33dpx-2*dxt)*den;
		dg11dpz = -(gamma33-1)*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den+
			(dgamma33dpz-2*dzt)*den;
		dg33dpx = -(gamma11-1)*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den+
			(dgamma11dpx-2*dxt)*den;
		dg33dpz = -(gamma11-1)*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den+
			(dgamma11dpz-2*dzt)*den;
		dg13dpx = gamma13*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den-
			dgamma13dpx*den;
		dg13dpz = gamma13*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den-
			dgamma13dpz*den;
		dg11dn = -(gamma33-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma33dn-snfact)*den;
		dg33dn = -(gamma11-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma11dn-snfact)*den;
		dg13dn = gamma13*(dgamma11dn+dgamma33dn-2*snfact)*den*den-
			dgamma13dn*den;
		ddxdpx=   a1111*pxt*dg11dpx+(a1133+a1313)*pzt*dg13dpx+
			a3313*pzt*dg33dpx+a1113*(pzt*dg11dpx+2*pxt*dg13dpx)
			+a1313*dg33dpx*pxt;
		ddzdpz= a3333*pzt*dg33dpz+(a1133+a1313)*pxt*dg13dpz+
			a1113*pxt*dg11dpz+a3313*(pxt*dg33dpz+2*pzt*dg13dpz)+
			a1313*dg11dpz*pzt;
		ddxdpz= a1111*pxt*dg11dpz+(a1133+a1313)*pzt*dg13dpz+
			a3313*pzt*dg33dpz+a1113*(pzt*dg11dpz+2*pxt*dg13dpz)+
			a1313*dg33dpz*pxt;
		dsnfactordn = da1111dn*px2*dg11dn+da3333dn*pz2*dg33dn+
			2*(da1133dn+da1313dn)*pxz*dg13dn+da1313dn*(px2*dg33dn+pz2*dg11dn)+
			2*da3313dn*(pz2*dg13dn+pxz*dg33dn)+2*da1113dn*(pxz*dg11dn+px2*dg13dn);
		ddxdn =  (a1111*pxt*dg11dn+(a1133+a1313)*pzt*dg13dn+a3313*pzt*dg33dn
			+a1113*(pzt*dg11dn+2*pxt*dg13dn)+a1313*dg33dn*pxt);
		ddzdn =  (a3333*pzt*dg33dn+(a1133+a1313)*pxt*dg13dn+a1113*pxt*dg11dn
			+a3313*(pxt*dg33dn+2*pzt*dg13dn)+a1313*dg11dn*pzt);
		
		dcdp1  = a1111*g11+a1313*g33+2*a1113*g13+ddxdpx-dxt*dxt;
		dcdp3  = a3333*g33+a1313*g11+2*a3313*g13+ddzdpz-dzt*dzt;
		dcdp13 = a1133*g13+a1313*g13+a1113*g11+a3313*g33+ddxdpz-dxt*dzt;
		ddcdqq = dcdp1*cc-2.0*dcdp13*s*c+dcdp3*ss;
		dxdn   =  (da1111dn*pxt*g11+(da1133dn+da1313dn)*pzt*g13+
			da3313dn*pzt*g33+da1113dn*(pzt*g11+2*pxt*g13)+
			da1313dn*g33*pxt);
		dzdn   =  (da3333dn*pzt*g33+(da1133dn+da1313dn)*pxt*g13+
			da1113dn*pxt*g11+da3313dn*(pxt*g33+2*pzt*g13)+
			da1313dn*g11*pzt);
		ddcdpn = dxdn*c-dzdn*s-.5*dxt*sxfactor*cc+
			.5*(dxt*szfactor+dzt*sxfactor)*s*c-.5*dzt*szfactor*ss
			+ddxdn*c-ddzdn*s;
		snfactor = dda1111dndn*px2*g11+dda3333dndn*pz2*g33+
			2*(dda1133dndn+dda1313dndn)*pxz*g13+
			dda1313dndn*(px2*g33+pz2*g11)+
			2*dda3313dndn*(pz2*g13+pxz*g33)+
			2*dda1113dndn*(pxz*g11+px2*g13);
		ddcdnn = 0.5*snfactor-.25*sxfactor*sxfactor*cc+
			.5*sxfactor*szfactor*s*c-.25*szfactor*szfactor*ss
			+.5*dsnfactordn;

		dp1t = -ddcdnn*q1t-ddcdpn*p1t;
		dq1t = ddcdqq*p1t+ddcdpn*q1t;
		dp2t = -ddcdnn*q2t-ddcdpn*p2t;
		dq2t = ddcdqq*p2t+ddcdpn*q2t;
		xt = x+hhalf*dxt;
		zt = z+hhalf*dzt;
		pxt = px+hhalf*dpxt;
		pzt = pz+hhalf*dpzt;
		p1t = p1+hhalf*dp1t;
		q1t = q1+hhalf*dq1t;
		p2t = p2+hhalf*dp2t;
		q2t = q2+hhalf*dq2t;
		vp  = 1/sqrt(pxt*pxt+pzt*pzt);
		s   = pxt*vp;
		c   = pzt*vp;
		ss  = s*s;
		cc  = c*c;
		
		vel2Interp(a11112,xt,zt,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
			&dda1111dxdz,&dda1111dzdz);
		da1111dn    = da1111dx*c-da1111dz*s;
		dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

		vel2Interp(a33332,xt,zt,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
		da3333dn    = da3333dx*c-da3333dz*s;
		dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
		vel2Interp(a11332,xt,zt,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
			&dda1133dxdz,&dda1133dzdz);
		da1133dn    = da1133dx*c-da1133dz*s;
		dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

		vel2Interp(a13132,xt,zt,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
			&dda1313dxdz,&dda1313dzdz);
		da1313dn    = da1313dx*c-da1313dz*s;
		dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

		vel2Interp(a11132,xt,zt,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
			&dda1113dxdz,&dda1113dzdz);
		da1113dn    = da1113dx*c-da1113dz*s;
		dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

		vel2Interp(a33132,xt,zt,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
			&dda3313dxdz,&dda3313dzdz);
		da3313dn    = da3313dx*c-da3313dz*s;
		dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;
		
	/* step 3 of 4th-order Runge-Kutta */
		px2   = pxt*pxt;
		pz2   = pzt*pzt;
		pxz   = pxt*pzt;

		/*anisotropy parameters*/
		gamma11 = a1111*px2+ a1313*pz2 +2*a1113*pxz;
		gamma33 = a3333*pz2 + a1313*px2+2*a3313*pxz;
		gamma13 = (a1133+a1313)*pxz+ a1113*px2+ a3313*pz2;
		den     = 1/(gamma11+gamma33-2);
		g11     = (gamma33-1)*den;
		g33     = (gamma11-1)*den;
		g13     = -gamma13*den;
		sxfactor = da1111dx*px2*g11+da3333dx*pz2*g33+
			2*(da1133dx+da1313dx)*pxz*g13+da1313dx*(px2*g33+pz2*g11)+
			2*da3313dx*(pz2*g13+pxz*g33)+2*da1113dx*(pxz*g11+px2*g13);
		szfactor = da1111dz*px2*g11+da3333dz*pz2*g33+
			2*(da1133dz+da1313dz)*pxz*g13+da1313dz*(px2*g33+pz2*g11)+
			2*da3313dz*(pz2*g13+pxz*g33)+2*da1113dz*(pxz*g11+px2*g13);
		snfact = sxfactor*c-szfactor*s;
		
		/*computing ray velocities and derivatives*/
		dxm =  (a1111*pxt*g11+(a1133+a1313)*pzt*g13+a3313*pzt*g33
			+a1113*(pzt*g11+2*pxt*g13)+a1313*g33*pxt);
		dzm =  (a3333*pzt*g33+(a1133+a1313)*pxt*g13+a1113*pxt*g11
			+a3313*(pxt*g33+2*pzt*g13)+a1313*g11*pzt);
		dpxm = -0.5*sxfactor;
		dpzm = -0.5*szfactor;

		dgamma11dpx = 2*a1111*pxt+2*a1113*pzt;
		dgamma11dpz = 2*a1313*pzt+2*a1113*pxt;
		dgamma33dpx = 2*a1313*pxt+2*a3313*pzt;
		dgamma33dpz = 2*a3333*pzt+2*a3313*pxt;
		dgamma13dpx= (a1133+a1313)*pzt+2*a1113*pxt;
		dgamma13dpz= (a1133+a1313)*pxt+2*a3313*pzt;
		dgamma11dn = da1111dn*px2+ da1313dn*pz2 +2*da1113dn*pxz;
		dgamma33dn = da3333dn*pz2 + da1313dn*px2+2*da3313dn*pxz;
		dgamma13dn = (da1133dn+da1313dn)*pxz+ da1113dn*px2+ da3313dn*pz2;
		dg11dpx = -(gamma33-1)*(dgamma11dpx+dgamma33dpx-4*dxm)*den*den+
			(dgamma33dpx-2*dxm)*den;
		dg11dpz = -(gamma33-1)*(dgamma11dpz+dgamma33dpz-4*dzm)*den*den+
			(dgamma33dpz-2*dzm)*den;
		dg33dpx = -(gamma11-1)*(dgamma11dpx+dgamma33dpx-4*dxm)*den*den+
			(dgamma11dpx-2*dxm)*den;
		dg33dpz = -(gamma11-1)*(dgamma11dpz+dgamma33dpz-4*dzm)*den*den+
			(dgamma11dpz-2*dzm)*den;
		dg13dpx = gamma13*(dgamma11dpx+dgamma33dpx-4*dxm)*den*den-
			dgamma13dpx*den;
		dg13dpz = gamma13*(dgamma11dpz+dgamma33dpz-4*dzm)*den*den-
			dgamma13dpz*den;
		dg11dn = -(gamma33-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma33dn-snfact)*den;
		dg33dn = -(gamma11-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma11dn-snfact)*den;
		dg13dn = gamma13*(dgamma11dn+dgamma33dn-2*snfact)*den*den-
			dgamma13dn*den;
		ddxdpx=   a1111*pxt*dg11dpx+(a1133+a1313)*pzt*dg13dpx+
			a3313*pzt*dg33dpx+a1113*(pzt*dg11dpx+2*pxt*dg13dpx)
			+a1313*dg33dpx*pxt;
		ddzdpz= a3333*pzt*dg33dpz+(a1133+a1313)*pxt*dg13dpz+
			a1113*pxt*dg11dpz+a3313*(pxt*dg33dpz+2*pzt*dg13dpz)+
			a1313*dg11dpz*pzt;
		ddxdpz= a1111*pxt*dg11dpz+(a1133+a1313)*pzt*dg13dpz+
			a3313*pzt*dg33dpz+a1113*(pzt*dg11dpz+2*pxt*dg13dpz)+
			a1313*dg33dpz*pxt;
		dsnfactordn = da1111dn*px2*dg11dn+da3333dn*pz2*dg33dn+
			2*(da1133dn+da1313dn)*pxz*dg13dn+da1313dn*(px2*dg33dn+pz2*dg11dn)+
			2*da3313dn*(pz2*dg13dn+pxz*dg33dn)+2*da1113dn*(pxz*dg11dn+px2*dg13dn);
		ddxdn =  (a1111*pxt*dg11dn+(a1133+a1313)*pzt*dg13dn+a3313*pzt*dg33dn
			+a1113*(pzt*dg11dn+2*pxt*dg13dn)+a1313*dg33dn*pxt);
		ddzdn =  (a3333*pzt*dg33dn+(a1133+a1313)*pxt*dg13dn+a1113*pxt*dg11dn
			+a3313*(pxt*dg33dn+2*pzt*dg13dn)+a1313*dg11dn*pzt);

		
		dcdp1  = a1111*g11+a1313*g33+2*a1113*g13+ddxdpx-dxm*dxm;
		dcdp3  = a3333*g33+a1313*g11+2*a3313*g13+ddzdpz-dzm*dzm;
		dcdp13 = a1133*g13+a1313*g13+a1113*g11+a3313*g33+ddxdpz-dxm*dzm;
		ddcdqq = dcdp1*cc-2.0*dcdp13*s*c+dcdp3*ss;
		dxdn   =  (da1111dn*pxt*g11+(da1133dn+da1313dn)*pzt*g13+
			da3313dn*pzt*g33+da1113dn*(pzt*g11+2*pxt*g13)+
			da1313dn*g33*pxt);
		dzdn   =  (da3333dn*pzt*g33+(da1133dn+da1313dn)*pxt*g13+
			da1113dn*pxt*g11+da3313dn*(pxt*g33+2*pzt*g13)+
			da1313dn*g11*pzt);
		ddcdpn = dxdn*c-dzdn*s-.5*dxm*sxfactor*cc+
			.5*(dxm*szfactor+dzm*sxfactor)*s*c-.5*dzm*szfactor*ss
			+ddxdn*c-ddzdn*s;
		snfactor = dda1111dndn*px2*g11+dda3333dndn*pz2*g33+
			2*(dda1133dndn+dda1313dndn)*pxz*g13+
			dda1313dndn*(px2*g33+pz2*g11)+
			2*dda3313dndn*(pz2*g13+pxz*g33)+
			2*dda1113dndn*(pxz*g11+px2*g13);
		ddcdnn = 0.5*snfactor-.25*sxfactor*sxfactor*cc+
			.5*sxfactor*szfactor*s*c-.25*szfactor*szfactor*ss
			+.5*dsnfactordn;

		dp1m = -ddcdnn*q1t-ddcdpn*p1t;
		dq1m = ddcdqq*p1t+ddcdpn*q1t;
		dp2m = -ddcdnn*q2t-ddcdpn*p2t;
		dq2m = ddcdqq*p2t+ddcdpn*q2t;
		xt = x+hhalf*dx;
		zt = z+hhalf*dz;
		pxt = px+h*dpxm;
		pzt = pz+h*dpzm;
		p1t = p1+h*dp1m;
		q1t = q1+h*dq1m;
		p2t = p2+h*dp2m;
		q2t = q2+h*dq2m;
		dxm += dxt;
		dzm += dzt;
		dpxm += dpxt;
		dpzm += dpzt;
		dp1m += dp1t;
		dq1m += dq1t;
		dp2m += dp2t;
		dq2m += dq2t;
		vp  = 1/sqrt(pxt*pxt+pzt*pzt);
		s   = pxt*vp;
		c   = pzt*vp;
		ss  = s*s;
		cc  = c*c;
		
		vel2Interp(a11112,xt,zt,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
			&dda1111dxdz,&dda1111dzdz);
		da1111dn    = da1111dx*c-da1111dz*s;
		dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

		vel2Interp(a33332,xt,zt,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
		da3333dn    = da3333dx*c-da3333dz*s;
		dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
		vel2Interp(a11332,xt,zt,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
			&dda1133dxdz,&dda1133dzdz);
		da1133dn    = da1133dx*c-da1133dz*s;
		dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

		vel2Interp(a13132,xt,zt,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
			&dda1313dxdz,&dda1313dzdz);
		da1313dn    = da1313dx*c-da1313dz*s;
		dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

		vel2Interp(a11132,xt,zt,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
			&dda1113dxdz,&dda1113dzdz);
		da1113dn    = da1113dx*c-da1113dz*s;
		dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

		vel2Interp(a33132,xt,zt,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
			&dda3313dxdz,&dda3313dzdz);
		da3313dn    = da3313dx*c-da3313dz*s;
		dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;
		
	/* step 4 of 4th-order Runge-Kutta */
		px2   = pxt*pxt;
		pz2   = pzt*pzt;
		pxz   = pxt*pzt;

		/*anisotropy parameters*/
		gamma11 = a1111*px2+ a1313*pz2 +2*a1113*pxz;
		gamma33 = a3333*pz2 + a1313*px2+2*a3313*pxz;
		gamma13 = (a1133+a1313)*pxz+ a1113*px2+ a3313*pz2;
		den     = 1/(gamma11+gamma33-2);
		g11     = (gamma33-1)*den;
		g33     = (gamma11-1)*den;
		g13     = -gamma13*den;
		sxfactor = da1111dx*px2*g11+da3333dx*pz2*g33+
			2*(da1133dx+da1313dx)*pxz*g13+da1313dx*(px2*g33+pz2*g11)+
			2*da3313dx*(pz2*g13+pxz*g33)+2*da1113dx*(pxz*g11+px2*g13);
		szfactor = da1111dz*px2*g11+da3333dz*pz2*g33+
			2*(da1133dz+da1313dz)*pxz*g13+da1313dz*(px2*g33+pz2*g11)+
			2*da3313dz*(pz2*g13+pxz*g33)+2*da1113dz*(pxz*g11+px2*g13);
		snfact = sxfactor*c-szfactor*s;
		
		/*computing ray velocities and derivatives*/
		dxt =  (a1111*pxt*g11+(a1133+a1313)*pzt*g13+a3313*pzt*g33
			+a1113*(pzt*g11+2*pxt*g13)+a1313*g33*pxt);
		dzt =  (a3333*pzt*g33+(a1133+a1313)*pxt*g13+a1113*pxt*g11
			+a3313*(pxt*g33+2*pzt*g13)+a1313*g11*pzt);
		dpxt = -0.5*sxfactor;
		dpzt = -0.5*szfactor;

		dgamma11dpx = 2*a1111*pxt+2*a1113*pzt;
		dgamma11dpz = 2*a1313*pzt+2*a1113*pxt;
		dgamma33dpx = 2*a1313*pxt+2*a3313*pzt;
		dgamma33dpz = 2*a3333*pzt+2*a3313*pxt;
		dgamma13dpx= (a1133+a1313)*pzt+2*a1113*pxt;
		dgamma13dpz= (a1133+a1313)*pxt+2*a3313*pzt;
		dgamma11dn = da1111dn*px2+ da1313dn*pz2 +2*da1113dn*pxz;
		dgamma33dn = da3333dn*pz2 + da1313dn*px2+2*da3313dn*pxz;
		dgamma13dn = (da1133dn+da1313dn)*pxz+ da1113dn*px2+ da3313dn*pz2;
		dg11dpx = -(gamma33-1)*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den+
			(dgamma33dpx-2*dxt)*den;
		dg11dpz = -(gamma33-1)*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den+
			(dgamma33dpz-2*dzt)*den;
		dg33dpx = -(gamma11-1)*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den+
			(dgamma11dpx-2*dxt)*den;
		dg33dpz = -(gamma11-1)*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den+
			(dgamma11dpz-2*dzt)*den;
		dg13dpx = gamma13*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den-
			dgamma13dpx*den;
		dg13dpz = gamma13*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den-
			dgamma13dpz*den;
		dg11dn = -(gamma33-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma33dn-snfact)*den;
		dg33dn = -(gamma11-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma11dn-snfact)*den;
		dg13dn = gamma13*(dgamma11dn+dgamma33dn-2*snfact)*den*den-
			dgamma13dn*den;
		ddxdpx=   a1111*pxt*dg11dpx+(a1133+a1313)*pzt*dg13dpx+
			a3313*pzt*dg33dpx+a1113*(pzt*dg11dpx+2*pxt*dg13dpx)
			+a1313*dg33dpx*pxt;
		ddzdpz= a3333*pzt*dg33dpz+(a1133+a1313)*pxt*dg13dpz+
			a1113*pxt*dg11dpz+a3313*(pxt*dg33dpz+2*pzt*dg13dpz)+
			a1313*dg11dpz*pzt;
		ddxdpz= a1111*pxt*dg11dpz+(a1133+a1313)*pzt*dg13dpz+
			a3313*pzt*dg33dpz+a1113*(pzt*dg11dpz+2*pxt*dg13dpz)+
			a1313*dg33dpz*pxt;
		dsnfactordn = da1111dn*px2*dg11dn+da3333dn*pz2*dg33dn+
			2*(da1133dn+da1313dn)*pxz*dg13dn+da1313dn*(px2*dg33dn+pz2*dg11dn)+
			2*da3313dn*(pz2*dg13dn+pxz*dg33dn)+2*da1113dn*(pxz*dg11dn+px2*dg13dn);
		ddxdn =  (a1111*pxt*dg11dn+(a1133+a1313)*pzt*dg13dn+a3313*pzt*dg33dn
			+a1113*(pzt*dg11dn+2*pxt*dg13dn)+a1313*dg33dn*pxt);
		ddzdn =  (a3333*pzt*dg33dn+(a1133+a1313)*pxt*dg13dn+a1113*pxt*dg11dn
			+a3313*(pxt*dg33dn+2*pzt*dg13dn)+a1313*dg11dn*pzt);

		
		dcdp1  = a1111*g11+a1313*g33+2*a1113*g13+ddxdpx-dxt*dxt;
		dcdp3  = a3333*g33+a1313*g11+2*a3313*g13+ddzdpz-dzt*dzt;
		dcdp13 = a1133*g13+a1313*g13+a1113*g11+a3313*g33+ddxdpz-dxt*dzt;
		ddcdqq = dcdp1*cc-2.0*dcdp13*s*c+dcdp3*ss;
		dxdn   =  (da1111dn*pxt*g11+(da1133dn+da1313dn)*pzt*g13+
			da3313dn*pzt*g33+da1113dn*(pzt*g11+2*pxt*g13)+
			da1313dn*g33*pxt);
		dzdn   =  (da3333dn*pzt*g33+(da1133dn+da1313dn)*pxt*g13+
			da1113dn*pxt*g11+da3313dn*(pxt*g33+2*pzt*g13)+
			da1313dn*g11*pzt);
		ddcdpn = dxdn*c-dzdn*s-.5*dxt*sxfactor*cc+
			.5*(dxt*szfactor+dzt*sxfactor)*s*c-.5*dzt*szfactor*ss
			+ddxdn*c-ddzdn*s;
		snfactor = dda1111dndn*px2*g11+dda3333dndn*pz2*g33+
			2*(dda1133dndn+dda1313dndn)*pxz*g13+
			dda1313dndn*(px2*g33+pz2*g11)+
			2*dda3313dndn*(pz2*g13+pxz*g33)+
			2*dda1113dndn*(pxz*g11+px2*g13);
		ddcdnn = 0.5*snfactor-.25*sxfactor*sxfactor*cc+
			.5*sxfactor*szfactor*s*c-.25*szfactor*szfactor*ss
			+.5*dsnfactordn;

		dp1t = -ddcdnn*q1t-ddcdpn*p1t;
		dq1t = ddcdqq*p1t+ddcdpn*q1t;
		dp2t = -ddcdnn*q2t-ddcdpn*p2t;
		dq2t = ddcdqq*p2t+ddcdpn*q2t;
		dxx  = hsixth*(dx+dxt+2.0*dxm);
		dzz  = hsixth*(dz+dzt+2.0*dzm);
		x += dxx;
		z += dzz;
		px += hsixth*(dpx+dpxt+2.0*dpxm);
		pz += hsixth*(dpz+dpzt+2.0*dpzm);
		p1 += hsixth*(dp1+dp1t+2.0*dp1m);
		q1 += hsixth*(dq1+dq1t+2.0*dq1m);
		p2 += hsixth*(dp2+dp2t+2.0*dp2m);
		q2 += hsixth*(dq2+dq2t+2.0*dq2m);
		vp  = 1/sqrt(px*px+pz*pz);
		s   = px*vp;
		c   = pz*vp;
		ss  = s*s;
		cc  = c*c;
		
		vel2Interp(a11112,x,z,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
			&dda1111dxdz,&dda1111dzdz);
		da1111dn    = da1111dx*c-da1111dz*s;
		dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

		vel2Interp(a33332,x,z,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
		da3333dn    = da3333dx*c-da3333dz*s;
		dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
		vel2Interp(a11332,x,z,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
			&dda1133dxdz,&dda1133dzdz);
		da1133dn    = da1133dx*c-da1133dz*s;
		dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

		vel2Interp(a13132,x,z,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
			&dda1313dxdz,&dda1313dzdz);
		da1313dn    = da1313dx*c-da1313dz*s;
		dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

		vel2Interp(a11132,x,z,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
			&dda1113dxdz,&dda1113dzdz);
		da1113dn    = da1113dx*c-da1113dz*s;
		dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

		vel2Interp(a33132,x,z,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
			&dda3313dxdz,&dda3313dzdz);
		da3313dn    = da3313dx*c-da3313dz*s;
		dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;


		/* update time */
		t  += dt;

		/* update kmah index */
                if ((q2<=0.0 && q2old>0.0) || (q2>=0.0 && q2old<0.0)) kmah++;

		/* save ray parameters */
		rs[it].t = t;
		rs[it].x = x;
		rs[it].z = z;
		rs[it].c = c;
		rs[it].s = s;
		rs[it].q1 = q1;
		rs[it].p1 = p1;
		rs[it].q2 = q2;
		rs[it].p2 = p2;
		rs[it].kmah = kmah;
		rs[it].v = vp;
		rs[it].dvdx = .5*da3333dx*vp/a3333;
		rs[it].dvdz = .5*da3333dz*vp/a3333;
	/*printf("%d %f %f %f %f %f %f %f \n",it,t,x,z,fx,fz,lx,lz);*/
	/*printf("%f %f \n",z,x);*/
	}


	/* free velocity interpolator */
	vel2Free(a11112);
	vel2Free(a33332);
	vel2Free(a11332);
	vel2Free(a13132);
	vel2Free(a11132);
	vel2Free(a33132);

	/* return ray */
	ray->nrs = it;
	ray->rs = rs;
	ray->nc = 0;
	ray->c = NULL;
	return ray;
}



float dvdlam(int para_flag,float VP0,float Vh,float Vn,float nu,float c,float s,float V)

/*****************************************************************************
**  compute the derivative of phase velocity with respect to one parameter  **
******************************************************************************	
Input:
para_flag	flag of specific parameter, =1 for VP0, =2 for Vh, =3 for Vn
VP0		symmetry-direction P-wave velocity at the ray step
Vh		horizontal velocity at the ray step
Vn		NMO velocity at the ray step
nu 		tilt value at the ray step
c		cosine of phase angle with respect to the vertical 
s		sine of phase angle with respect to the vertical 
V		phase velocity at the ray step
*****************************************************************************/
{	
	float pha_an,pha_v,vpe2,ss,cc,dvdvp,dvdvh,dvdvn,sqroot;
	
	pha_an = -atan(s/c)-nu;
	
	ss = sin(pha_an);
	cc = cos(pha_an);
	
	vpe2 = pow(Vh,2)*pow(ss,2) + pow(VP0,2)*pow(cc,2); 
	
	pha_v = sqrt(vpe2*(1+sqrt(1+4*pow(VP0,2)*(pow(Vn,2)-pow(Vh,2))*pow(ss,2)*pow(cc,2)/pow(vpe2,2)))/2);
	
	if (abs((pha_v-V)/V)>0.02)
		printf("Because the difference between analytical value of phase velocity and the value obtained from ray tracing is too large, the errors in the computation of derivatives may be large. \n");
	
	sqroot = sqrt(pow(vpe2,2)+4*pow(VP0,2)*(pow(Vn,2)-pow(Vh,2))*pow(ss,2)*pow(cc,2));
	
	switch(para_flag)
	{	
		case 1:	dvdvp = (VP0*pow(cc,2)+(vpe2*VP0*pow(cc,2)+2*VP0*(pow(Vn,2)-pow(Vh,2))*pow(ss,2)*pow(cc,2))/sqroot)/(2*pha_v); return dvdvp; break;
		case 2: dvdvh = (Vh*pow(ss,2)+(vpe2*Vh*pow(ss,2)-2*Vh*pow(VP0,2)*pow(ss,2)*pow(cc,2))/sqroot)/(2*pha_v); return dvdvh; break;
		case 3: dvdvn = Vn*pow(VP0,2)*pow(ss,2)*pow(cc,2)/(pha_v*sqroot); return dvdvn; break;	
	}

}



float lagrange_scalar(int point_flag,int ix,int iz,float x,float z,float fx,float fz,float dx,float dz)

/*****************************************************************************
***************  compute the scalar of Lagrange interpolation  ***************
******************************************************************************	
Input:
point_flag	flag of specific vertex, 1 for (ix,iz), 2 for (ix+1,iz),
		3 for (ix,iz+1), and 4 for (ix+1,iz+1)
ix, iz		the index of upper left vertex
x, z 		the coordinates of the ray step
fx, fz 		the first sample of x and z coordinates
dx, dz 		the interval of x and z sampling 
*****************************************************************************/

{	
	float ls,x1,z1,x2,z2,x3,z3,x4,z4;
	
	x1 = fx + ix*dx;
	z1 = fz + iz*dz;
	
	x2 = x1+dx;
	z2 = z1;
	
	x3 = x1;
	z3 = z1+dz;
	
	x4 = x2;
	z4 = z3;
	
	switch(point_flag)
	{
		case 1: ls = (sqrt(pow(x-x2,2)+pow(z-z2,2))/dx)*(sqrt(pow(x-x3,2)+pow(z-z3,2))/dz)*(sqrt(pow(x-x4,2)+pow(z-z4,2))/sqrt(pow(dx,2)+pow(dz,2))); return ls; break;
		
		case 2: ls = (sqrt(pow(x-x1,2)+pow(z-z1,2))/dx)*(sqrt(pow(x-x4,2)+pow(z-z4,2))/dz)*(sqrt(pow(x-x3,2)+pow(z-z3,2))/sqrt(pow(dx,2)+pow(dz,2))); return ls; break;
		
		case 3: ls = (sqrt(pow(x-x4,2)+pow(z-z4,2))/dx)*(sqrt(pow(x-x1,2)+pow(z-z1,2))/dz)*(sqrt(pow(x-x2,2)+pow(z-z2,2))/sqrt(pow(dx,2)+pow(dz,2))); return ls; break;
		
		case 4: ls = (sqrt(pow(x-x3,2)+pow(z-z3,2))/dx)*(sqrt(pow(x-x2,2)+pow(z-z2,2))/dz)*(sqrt(pow(x-x1,2)+pow(z-z1,2))/sqrt(pow(dx,2)+pow(dz,2))); return ls; break;	
	}
				 
}


void freeRay (Ray *ray)
/*****************************************************************************
Free a ray.
******************************************************************************
Input:
ray		ray to be freed
*****************************************************************************/
{
	if (ray->c!=NULL) free1((void*)ray->c);
	free1((void*)ray->rs);
	free1((void*)ray);
}


/*****************************************************************************
Functions to support interpolation of velocity and its derivatives.
******************************************************************************
Functions:
vel2Alloc	allocate and initialize an interpolator for v(x,z)
vel2Interp	interpolate v(x,z) and its derivatives
******************************************************************************
Notes:
Interpolation is performed by piecewise cubic Hermite polynomials,
so that velocity and first derivatives are continuous.  Therefore,
velocity v, first derivatives dv/dx and dv/dz, and the mixed
derivative ddv/dxdz are continuous.  However, second derivatives
ddv/dxdx and ddv/dzdz are discontinuous.
*****************************************************************************
Technical Reference:
	Hale, D., 1992, Migration by the Kirchhoff, 
	slant stack, and Gaussian beam methods:
	Colorado School of Mines.
*****************************************************************************
 Credits: 	CWP: Dave Hale
*****************************************************************************/

/* number of pre-computed, tabulated interpolators */
#define NTABLE 101

/* length of each interpolator in table (4 for piecewise cubic) */
#define LTABLE 4

/* table of pre-computed interpolators, for 0th, 1st, and 2nd derivatives */
static float tbl[3][NTABLE][LTABLE];

/* constants */
static int ix=1-LTABLE/2-LTABLE,iz=1-LTABLE/2-LTABLE;
static float ltable=LTABLE,ntblm1=NTABLE-1;

/* indices for 0th, 1st, and 2nd derivatives */
static int kx[6]={0,1,0,2,1,0};
static int kz[6]={0,0,1,0,1,2};

/* function to build interpolator tables; sets tabled=1 when built */
static void buildTables (void);
static int tabled=0;

/* interpolator for velocity function v(x,z) of two variables */
typedef struct Vel2Struct {
	int nx;		/* number of x samples */
	int nz;		/* number of z samples */
	int nxm;	/* number of x samples minus LTABLE */
	int nzm;	/* number of x samples minus LTABLE */
	float xs,xb,zs,zb,sx[3],sz[3],**vu;
} Vel2;

void* vel2Alloc (int nx, float dx, float fx,
	int nz, float dz, float fz, float **v)
/*****************************************************************************
Allocate and initialize an interpolator for v(x,z) and its derivatives.
******************************************************************************
Input:
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
v		array[nx][nz] of uniformly sampled v(x,z)

*****************************************************************************
Returned:	pointer to interpolator
*****************************************************************************/
{
	Vel2 *vel2;

	/* allocate space */
	vel2 = (Vel2*)alloc1(1,sizeof(Vel2));

	/* set state variables used for interpolation */
	vel2->nx = nx;
	vel2->nxm = nx-LTABLE;
	vel2->xs = 1.0/dx;
	vel2->xb = ltable-fx*vel2->xs;
	vel2->sx[0] = 1.0;
	vel2->sx[1] = vel2->xs;
	vel2->sx[2] = vel2->xs*vel2->xs;
	vel2->nz = nz;
	vel2->nzm = nz-LTABLE;
	vel2->zs = 1.0/dz;
	vel2->zb = ltable-fz*vel2->zs;
	vel2->sz[0] = 1.0;
	vel2->sz[1] = vel2->zs;
	vel2->sz[2] = vel2->zs*vel2->zs;
	vel2->vu = v;
	
	/* if necessary, build interpolator coefficient tables */
	if (!tabled) buildTables();

	return vel2;
}

void vel2Free (void *vel2)
/*****************************************************************************
Free an interpolator for v(x,z) and its derivatives.
******************************************************************************
Input:
vel2		pointer to interpolator as returned by vel2Alloc()
*****************************************************************************/
{
	free1(vel2);
}

void vel2Interp (void *vel2, float x, float z,
	float *v, float *vx, float *vz, float *vxx, float *vxz, float *vzz)
/*****************************************************************************
Interpolation of a velocity function v(x,z) and its derivatives.
******************************************************************************
Input:
vel2		pointer to interpolator as returned by vel2Alloc()
x		x coordinate at which to interpolate v(x,z) and derivatives
z		z coordinate at which to interpolate v(x,z) and derivatives
*****************************************************************************
Output:
v		v(x,z)
vx		dv/dx
vz		dv/dz
vxx		ddv/dxdx
vxz		ddv/dxdz
vzz		ddv/dzdz
*****************************************************************************/
{
	Vel2 *v2=vel2;
	int nx=v2->nx,nz=v2->nz,nxm=v2->nxm,nzm=v2->nzm;
	float xs=v2->xs,xb=v2->xb,zs=v2->zs,zb=v2->zb,
		*sx=v2->sx,*sz=v2->sz,**vu=v2->vu;
	int i,jx,lx,mx,jz,lz,mz,jmx,jmz,mmx,mmz;
	float ax,bx,*px,az,bz,*pz,sum,vd[6];

	/* determine offsets into vu and interpolation coefficients */
	ax = xb+x*xs;
	jx = (int)ax;
	bx = ax-jx;
	lx = (bx>=0.0)?bx*ntblm1+0.5:(bx+1.0)*ntblm1-0.5;
	lx *= LTABLE;
	mx = ix+jx;
	az = zb+z*zs;
	jz = (int)az;
	bz = az-jz;
	lz = (bz>=0.0)?bz*ntblm1+0.5:(bz+1.0)*ntblm1-0.5;
	lz *= LTABLE;
	mz = iz+jz;
	
	/* if totally within input array, use fast method */
	if (mx>=0 && mx<=nxm && mz>=0 && mz<=nzm) {
		for (i=0; i<6; ++i) {
			px = &(tbl[kx[i]][0][0])+lx;
			pz = &(tbl[kz[i]][0][0])+lz;
			vd[i] = sx[kx[i]]*sz[kz[i]]*(
				vu[mx][mz]*px[0]*pz[0]+
				vu[mx][mz+1]*px[0]*pz[1]+
				vu[mx][mz+2]*px[0]*pz[2]+
				vu[mx][mz+3]*px[0]*pz[3]+
				vu[mx+1][mz]*px[1]*pz[0]+
				vu[mx+1][mz+1]*px[1]*pz[1]+
				vu[mx+1][mz+2]*px[1]*pz[2]+
				vu[mx+1][mz+3]*px[1]*pz[3]+
				vu[mx+2][mz]*px[2]*pz[0]+
				vu[mx+2][mz+1]*px[2]*pz[1]+
				vu[mx+2][mz+2]*px[2]*pz[2]+
				vu[mx+2][mz+3]*px[2]*pz[3]+
				vu[mx+3][mz]*px[3]*pz[0]+
				vu[mx+3][mz+1]*px[3]*pz[1]+
				vu[mx+3][mz+2]*px[3]*pz[2]+
				vu[mx+3][mz+3]*px[3]*pz[3]);
		}
		
	/* else handle end effects with constant extrapolation */
	} else {
		for (i=0; i<6; ++i) {
			px = &(tbl[kx[i]][0][0])+lx;
			pz = &(tbl[kz[i]][0][0])+lz;
			for (jx=0,jmx=mx,sum=0.0; jx<4; ++jx,++jmx) {
				mmx = jmx;
				if (mmx<0) mmx = 0;
				else if (mmx>=nx) mmx = nx-1;
				for (jz=0,jmz=mz; jz<4; ++jz,++jmz) {
					mmz = jmz;
					if (mmz<0) mmz = 0;
					else if (mmz>=nz) mmz = nz-1;
					sum += vu[mmx][mmz]*px[jx]*pz[jz];
				}
			}
			vd[i] = sx[kx[i]]*sz[kz[i]]*sum;
		}
	}

	/* set output variables */
	*v = vd[0];
	*vx = vd[1];
	*vz = vd[2];
	*vxx = vd[3];
	*vxz = vd[4];
	*vzz = vd[5];
}

/* hermite polynomials */
static float h00 (float x) {return 2.0*x*x*x-3.0*x*x+1.0;}
static float h01 (float x) {return 6.0*x*x-6.0*x;}
static float h02 (float x) {return 12.0*x-6.0;}
static float h10 (float x) {return -2.0*x*x*x+3.0*x*x;}
static float h11 (float x) {return -6.0*x*x+6.0*x;}
static float h12 (float x) {return -12.0*x+6.0;}
static float k00 (float x) {return x*x*x-2.0*x*x+x;}
static float k01 (float x) {return 3.0*x*x-4.0*x+1.0;}
static float k02 (float x) {return 6.0*x-4.0;}
static float k10 (float x) {return x*x*x-x*x;}
static float k11 (float x) {return 3.0*x*x-2.0*x;}
static float k12 (float x) {return 6.0*x-2.0;}

/* function to build interpolation tables */
static void buildTables(void)
{
	int i;
	float x;
	
	/* tabulate interpolator for 0th derivative */
	for (i=0; i<NTABLE; ++i) {
		x = (float)i/(NTABLE-1.0);
		tbl[0][i][0] = -0.5*k00(x);
		tbl[0][i][1] = h00(x)-0.5*k10(x);
		tbl[0][i][2] = h10(x)+0.5*k00(x);
		tbl[0][i][3] = 0.5*k10(x);
		tbl[1][i][0] = -0.5*k01(x);
		tbl[1][i][1] = h01(x)-0.5*k11(x);
		tbl[1][i][2] = h11(x)+0.5*k01(x);
		tbl[1][i][3] = 0.5*k11(x);
		tbl[2][i][0] = -0.5*k02(x);
		tbl[2][i][1] = h02(x)-0.5*k12(x);
		tbl[2][i][2] = h12(x)+0.5*k02(x);
		tbl[2][i][3] = 0.5*k12(x);
	}
	
	/* remember that tables have been built */
	tabled = 1;
}

#ifdef TEST2
#include "cwp.h"
#define NZ 2
#define NX 2
#define NXOUT 21
#define NZOUT 21
main()
{
	int nx=NX,nz=NZ,nxout=NXOUT,nzout=NZOUT,i,j;
	float dx=2.0,fx=-1.0,dxout=0.2,fxout=-2.0;
	float dz=4.0,fz=-2.0,dzout=0.4,fzout=-4.0;
	float x,z,v,vx,vz,vxx,vxz,vzz,**vu,**vo;
	void *vel2;

	vu = alloc2float(nz,nx);
	vo = alloc2float(nzout,nxout);

	vu[0][0] = 1.0;
	vu[1][0] = 2.0;
	vu[0][1] = 1.0;
	vu[1][1] = 2.0;

	vel2 = vel2Alloc(nx,dx,fx,nz,dz,fz,vu);

	for (i=0; i<nxout; ++i) {
		x = fxout+i*dxout;
		for (j=0; j<nzout; ++j) {
			z = fzout+j*dzout;
			vel2Interp(vel2,x,z,&v,&vx,&vz,&vxx,&vxz,&vzz);
			vo[i][j] = vz;
		}
	}
	vel2Free(vel2);
	fwrite(vo[0],sizeof(float),nxout*nzout,stdout);
}
#endif



void smooth2 (int n1, int n2, float r1, float r2, int *win, float rw, float **v)
/**************************************************************************
smooth2 --> smooths the parameter field **v depending on the values of r1 and r2
***************************************************************************
Input:
v       parameter field 
n1      number of samples in the fast direction 
n2	number of samples in the slow direction
r1	smoothing parameter in the fast direction
r2	smoothing parameter in the slow direction
win	1d array defining the corners of smoothing window 
rw	smoothing parameter for window 

*************************************************************************
*************************************************************************
Credits: CWP:  Zhenyue Liu, 1992 
************************************************************************/
{
	int nmax;	/* max of n1 and n2 */
	int ix, iz;	/* counters */
	float **w;	/* intermediate array */
	float *errz;	/* array of error estimates as a function of x1 */
	float *d, *e;	/* input arrays for subroutine tripd */
	float *f;	/* intermediate array */

	/* scale the smoothing parameter */
	r1 = r1*r1*0.25;
	r2 = r2*r2*0.25;

	/* allocate space */
	nmax = (n1<n2)?n2:n1;
	w = alloc2float(n1,n2);
	errz = alloc1float(nmax);
	d = alloc1float(nmax);
	e = alloc1float(nmax);
	f = alloc1float(nmax);

	rw = rw*rw*0.25;
 
	/* define the window function */
	for(ix=0; ix<n2; ++ix)
	 	for(iz=0; iz<n1; ++iz)
			w[ix][iz] = 0;	
	for(ix=win[2]; ix<win[3]; ++ix)
	 	for(iz=win[0]; iz<win[1]; ++iz)
			w[ix][iz] = 1;	

	if(win[0]>0 || win[1]<n1 || win[2]>0 || win[3]<n2){
	/*	smooth the window function */
         	for(iz=0; iz<n1; ++iz){
	 		for(ix=0; ix<n2; ++ix){
				d[ix] = 1.0+2.0*rw;
				e[ix] = -rw;
				f[ix] = w[ix][iz];
			}
        		d[0] -= rw;
         		d[n2-1] -= rw;
         		tripd(d,e,f,n2);
	 		for(ix=0; ix<n2; ++ix)
				w[ix][iz] = f[ix];
		}
         	for(ix=0; ix<n2; ++ix){
	 		for(iz=0; iz<n1; ++iz){
				d[iz] = 1.0+2.0*rw;
				e[iz] = -rw;
				f[iz] = w[ix][iz];
		}
        		d[0] -= rw;
         		d[n1-1] -= rw;
         		tripd(d,e,f,n1);
	 		for(iz=0; iz<n1; ++iz)
				w[ix][iz] = f[iz];
		}
	}

	/*      solving for the smoothing velocity */
        for(iz=0; iz<n1; ++iz){
	 	for(ix=0; ix<n2-1; ++ix){
			d[ix] = 1.0+r2*(w[ix][iz]+w[ix+1][iz]);
			e[ix] = -r2*w[ix+1][iz];
			f[ix] = v[ix][iz];
		}
        	d[0] -= r2*w[0][iz];
         	d[n2-1] = 1.0+r2*w[n2-1][iz];
		f[n2-1] = v[n2-1][iz];
         	tripd(d,e,f,n2);
	 	for(ix=0; ix<n2; ++ix)
			v[ix][iz] = f[ix];
	}
         for(ix=0; ix<n2; ++ix){
	 	for(iz=0; iz<n1-2; ++iz){
			d[iz] = 1.0+r1*(w[ix][iz+1]+w[ix][iz+2]);
			e[iz] = -r1*w[ix][iz+2];
			f[iz] = v[ix][iz+1];
		}
		f[0] += r1*w[ix][1]*v[ix][0];
         	d[n1-2] = 1.0+r1*w[ix][n1-1];
		f[n1-2] = v[ix][n1-1];
         	tripd(d,e,f,n1-1);
	 	for(iz=0; iz<n1-1; ++iz)
			v[ix][iz+1] = f[iz];
	}
	free2float(w);
	free1float(errz);
	free1float(d);
	free1float(e);
	free1float(f);

}
	
