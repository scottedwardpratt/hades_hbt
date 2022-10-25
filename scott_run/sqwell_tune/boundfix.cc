#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include "msu_commonutils/commonutils.h"

using namespace std;

int NQMAX=40;
double DELQ=2.0;
double m1=938.28;
double m2=2.0*m1-2.224;
int Q1Q2=0;

void SolveTanXoverX(double target,double &x);

void FixA(double BE,int nwells,vector<double> &V0,vector<double> &a);

void SquareWell_Init(int ell,int nwells,vector<double> &V0,vector<double> &a,vector<double> &delta,double &scatt_length);

int main(int argc,char *argv[]){
	if(argc!=7){
		printf("argc=%d, but must be 7\n",argc);
		exit(1);
	}
	int nwells=3,ell=0;
	vector<double> V0;
	vector<double> a;
	vector<double> delta;
	double scatt_length;
	V0.resize(nwells);
	a.resize(nwells+1);
	delta.resize(NQMAX);
	V0[0]=atof(argv[1]);
	V0[1]=atof(argv[2]);
	V0[2]=atof(argv[3]);
	a[0]=atof(argv[4]);
	a[1]=atof(argv[5]);
	a[2]=atof(argv[6]);
	a[3]=1.0E20;
	printf("V=(%g,%g,%g), a=(%g,%g,%g)\n",V0[0],V0[1],V0[2],a[0],a[1],a[2]);
	
	double target,x;
	for(int itest=0;itest<10;itest++){
		printf("Enter target: ");
		scanf("%lf",&target);
		SolveTanXoverX(target,x);
	}
	
	
	FixA(7.8,nwells,V0,a);
	return 0;
}

void SquareWell_Init(int ell,int nwells,vector<double> &V0,vector<double> &a,vector<double> &delta,double &scatt_length){
	// To set up the wave functions and phase shifts
	int nqmax=40;
	
	CGSLMatrix_Complex *cmatrix;
	double q,mu_coulomb,E,mu;
	double beta;
	double F1b,G1b,F1bprime,G1bprime;
	double F2a,G2a,F2aprime,G2aprime;
	double F2b,G2b,F2bprime,G2bprime;
	double F3a,G3a,F3aprime,G3aprime;
	double F3b,G3b,F3bprime,G3bprime;
	double F,G,Fprime,Gprime;
	double F1,G1,F1prime,G1prime;
	
	complex<double> eta0,eta1,eta2,eta3;
	complex<double> x1b,x2a,x2b,x3a,x3b,x,q1,q2,q3;
	
	vector<vector<complex<double>>> M;
	vector<complex<double>> Y;
	vector<vector<complex<double>>> A;
	
	complex<double> x1,x2;
	double F2,G2,F2prime,G2prime,qsquared,r;
	int i,j,iq,ir;
	mu=m1*m2/(m1+m2);
	
	if(nwells==1){
		for(iq=0;iq<NQMAX;iq++){
			q=(iq+0.5)*DELQ;
			E=sqrt(q*q+m1*m1)+sqrt(q*q+m2*m2);
			mu_coulomb=0.25*(E-pow(m1*m1-m2*m2,2)/pow(E,3));
			eta0=Q1Q2*mu_coulomb*ALPHA/q;
				
			qsquared=q*q-2.0*mu*V0[0];
			if(qsquared>0)
				q1=sqrt(qsquared);
			else
				q1=ci*sqrt(abs(qsquared));
			x1=q1*a[0]/HBARC;
			eta1=Q1Q2*mu_coulomb*ALPHA/q1;
			CoulWave::GetFGprime_ComplexQ(ell,x1,eta1,&F1,&G1,&F1prime,&G1prime);
			x2=q*a[0]/HBARC;
			CoulWave::GetFGprime_ComplexQ(ell,x2,eta0,&F2,&G2,&F2prime,&G2prime);     
			beta=(abs(q1)/q)*F1prime/F1;
			delta[iq]=-atan2(beta*F2-F2prime,beta*G2-G2prime);
			A[iq][0]=0.5*(exp(-2.0*ci*delta[iq])
				*(F2+ci*G2)+(F2-ci*G2))/F1;
			A[iq][1]=exp(-2.0*ci*delta[iq]);
				
		}
	}
	else if(nwells==2){
		Eigen::VectorXcd Y(4),A(4);
		Eigen::MatrixXcd M(4,4);
			
		for(iq=0;iq<NQMAX;iq++){
			q=(iq+0.5)*DELQ;
			E=sqrt(q*q+m1*m1)+sqrt(q*q+m2*m2);
			mu_coulomb=0.25*(E-pow(m1*m1-m2*m2,2)/pow(E,3));
				
			q1=sqrt(abs(q*q-2.0*mu*V0[0]));
			if(q*q-2.0*mu*V0[0]<0.0)
				q1=ci*q1;
			q2=sqrt(abs(q*q-2.0*mu*V0[1]));
			if(q*q-2.0*mu*V0[1]<0.0)
				q2=ci*q2;
			x1b=a[0]*q1/HBARC;
			x2a=a[0]*q2/HBARC;
			x2b=a[1]*q2/HBARC;
			x=a[1]*q/HBARC;
			eta1=Q1Q2*mu_coulomb*ALPHA/q1;
			eta2=Q1Q2*mu_coulomb*ALPHA/q2;
			eta0=Q1Q2*mu_coulomb*ALPHA/q;
			CoulWave::GetFGprime_ComplexQ(ell,x1b,eta1,&F1b,&G1b,&F1bprime,&G1bprime);
			CoulWave::GetFGprime_ComplexQ(ell,x2a,eta2,&F2a,&G2a,&F2aprime,&G2aprime);
			CoulWave::GetFGprime_ComplexQ(ell,x2b,eta2,&F2b,&G2b,&F2bprime,&G2bprime);
			CoulWave::GetFGprime_ComplexQ(ell,x,eta0,&F,&G,&Fprime,&Gprime);
			
			for(i=0;i<4;i++){
				Y[i]=0.0;
				for(j=0;j<4;j++)
					M(i,j)=0.0;
			}
			M(0,0)=F1b;              M(0,1)=-F2a;              M(0,2)=-G2a;
			M(1,0)=abs(q1)*F1bprime; M(1,1)=-abs(q2)*F2aprime; M(1,2)=-abs(q2)*G2aprime;
			M(2,1)=F2b;              M(2,2)=G2b;               M(2,3)=-0.5*(F+ci*G);
			M(3,1)=abs(q2)*F2bprime; M(3,2)=abs(q2)*G2bprime;  M(3,3)=-0.5*q*(Fprime+ci*Gprime);
				
			Y[2]=0.5*(F-ci*G); Y[3]=0.5*q*(Fprime-ci*Gprime);
			A=M.colPivHouseholderQr().solve(Y);
				
			delta[iq]=-0.5*atan2(imag(A(  3)),real(A(3)));
			if(delta[iq]<-PI)
				delta[iq]=delta[iq]+PI;
				
		}
	}
	else if(nwells==3){
		Eigen::VectorXcd Y(6),A(6);
		Eigen::MatrixXcd M(6,6);
			
		for(iq=0;iq<NQMAX;iq++){
			q=(iq+0.5)*DELQ;
			E=sqrt(q*q+m1*m1)+sqrt(q*q+m2*m2);
			mu_coulomb=0.25*(E-pow(m1*m1-m2*m2,2)/pow(E,3));
				
			q1=sqrt(abs(q*q-2.0*mu*V0[0]));
			if(q*q-2.0*mu*V0[0]<0.0) q1=ci*q1;
			q2=sqrt(abs(q*q-2.0*mu*V0[1]));
			if(q*q-2.0*mu*V0[1]<0.0) q2=ci*q2;
			q3=sqrt(abs(q*q-2.0*mu*V0[2]));
			if(q*q-2.0*mu*V0[2]<0.0) q3=ci*q3;
			x1b=a[0]*q1/HBARC;
			x2a=a[0]*q2/HBARC;
			x2b=a[1]*q2/HBARC;
			x3a=a[1]*q3/HBARC;
			x3b=a[2]*q3/HBARC;
			x=a[2]*q/HBARC;
			eta1=Q1Q2*mu_coulomb*ALPHA/q1;
			eta2=Q1Q2*mu_coulomb*ALPHA/q2;
			eta3=Q1Q2*mu_coulomb*ALPHA/q3;
			eta0=Q1Q2*mu_coulomb*ALPHA/q;
			CoulWave::GetFGprime_ComplexQ(ell,x1b,eta1,&F1b,&G1b,&F1bprime,&G1bprime);
			CoulWave::GetFGprime_ComplexQ(ell,x2a,eta2,&F2a,&G2a,&F2aprime,&G2aprime);
			CoulWave::GetFGprime_ComplexQ(ell,x2b,eta2,&F2b,&G2b,&F2bprime,&G2bprime);
			CoulWave::GetFGprime_ComplexQ(ell,x3a,eta3,&F3a,&G3a,&F3aprime,&G3aprime);
			CoulWave::GetFGprime_ComplexQ(ell,x3b,eta3,&F3b,&G3b,&F3bprime,&G3bprime);
			CoulWave::GetFGprime_ComplexQ(ell,x,eta0,&F,&G,&Fprime,&Gprime);
			
				
			for(i=0;i<6;i++){
				Y(i)=0.0;
				for(j=0;j<6;j++)
					M(i,j)=0.0;
			}
			M(0,0)=F1b;              M(0,1)=-F2a;              M(0,2)=-G2a;
			M(1,0)=abs(q1)*F1bprime; M(1,1)=-abs(q2)*F2aprime; M(1,2)=-abs(q2)*G2aprime;
			M(2,1)=F2b;              M(2,2)=G2b;               M(2,3)=-F3a;               M(2,4)=-G3a;
			M(3,1)=abs(q2)*F2bprime; M(3,2)=abs(q2)*G2bprime;  M(3,3)=-abs(q3)*F3aprime;  M(3,4)=-abs(q3)*G3aprime;
			M(4,3)=F3b;              M(4,4)=G3b;               M(4,5)=-0.5*(F+ci*G);
			M(5,3)=abs(q3)*F3bprime; M(5,4)=abs(q3)*G3bprime;  M(5,5)=-0.5*q*(Fprime+ci*Gprime);
			
				
			Y[4]=0.5*(F-ci*G); Y[5]=0.5*q*(Fprime-ci*Gprime);
			A=M.colPivHouseholderQr().solve(Y);	
			
			delta[iq]=-0.5*atan2(imag(A(5)),real(A(5)));
			if(delta[iq]<-PI)
				delta[iq]=delta[iq]+PI;
		}
	}
	
	iq=0;
	q=(iq+0.5)*DELQ;
	scatt_length=(delta[0]/(0.5*DELQ))*HBARC;
	printf("scattering length=%g\n",scatt_length);
	for(iq=0;iq<NQMAX;iq++){
		q=(iq+0.5)*DELQ;
		printf("%6.3f %8.3f\n",q,delta[iq]*180.0/PI);
	}

	/*
	for(iq=0;iq<NQMAX;iq++){
		DelPhiArray[iq][0]=0.0;
		for(ir=1;ir<=DelPhiArray_NRMAX;ir++){
			r=ir*DelPhiArray_DELR;
			SquareWell_CalcDelPhi(iq,r,DelPhiArray[iq][ir]);
		}
	}
	*/
	
}

void FixA(double BE,int nwells,vector<double> &V0,vector<double> &a){
	double mu=m1*m2/(m1+m2);
	double Q,lnderiv;
	int iwell;
	vector<double> q(nwells+1),KE(nwells+1),delta(nwells+1),A(nwells+1),psiborder(nwells+1);
	A[nwells]=1.0;
	KE[nwells]=-BE;
	q[nwells]=sqrt(2.0*mu*fabs(KE[nwells]))/HBARC;
	delta[0]=0.0;
	A[0]=1.0;
	delta[nwells]=0.0;
	for(iwell=0;iwell<nwells;iwell++){
		KE[iwell]=-BE-V0[iwell];
		q[iwell]=sqrt(2.0*mu*fabs(KE[iwell]))/HBARC;
	}
	for(iwell=0;iwell<nwells-1;iwell++){
		if(KE[iwell]>0.0){
			psiborder[iwell]=A[iwell]*sin(q[iwell]*a[iwell]+delta[iwell]);
			lnderiv=q[iwell]*cos(q[iwell]*a[iwell]+delta[iwell])/sin(q[iwell]*a[iwell]+delta[iwell]);
		}
		else{
			psiborder[iwell]=A[iwell]*sinh(q[iwell]*a[iwell]+delta[iwell]);
			lnderiv=q[iwell]*cosh(q[iwell]*a[iwell]+delta[iwell])/sinh(q[iwell]*a[iwell]+delta[iwell]);
		}
		printf("iwell=%d, lnderiv=%g\n",iwell,lnderiv);
		if(KE[iwell+1]>0){
			delta[iwell+1]=-q[iwell+1]*a[iwell]+atan(q[iwell]/lnderiv);
			A[iwell+1]=psiborder[iwell]/sin(q[iwell+1]*a[iwell]+delta[iwell+1]);
		}
		else{
			printf("x=%g,lnderiv=%g\n",q[iwell]/lnderiv,lnderiv);
			delta[iwell+1]=-q[iwell+1]*a[iwell]+atanh(q[iwell]/lnderiv);
			A[iwell+1]=psiborder[iwell]/sin(q[iwell+1]*a[iwell]+delta[iwell+1]);
		}
	}
	printf("iwell=%d, lnderiv=%g\n",nwells-1,lnderiv);
	if(KE[nwells-1]>0.0){
		A[nwells-1]=(-1.0/q[nwells-1])*atan(q[nwells-1]/q[nwells])-delta[nwells-1]/q[nwells-1];
	}
	else{
		A[nwells-1]=(-1.0/q[nwells-1])*atanh(q[nwells-1]/q[nwells])-delta[nwells-1]/q[nwells-1];
	}
	A[nwells]=psiborder[nwells-1];
	
	for(iwell=0;iwell<=nwells;iwell++){
		printf("%d: a=%6.2f KE=%7.3f, q=%7.3f, A=%7.3f, delta=%7.3f\n",iwell,a[iwell],KE[iwell],q[iwell],A[iwell],delta[iwell]);
	}
	
	printf("new amax=%g\n",a[nwells+1]);
	double r,delr=0.1,psi;
	for(r=0.0;r<5.0;r+=delr){
		iwell=0;
		while(r>a[iwell] && iwell<nwells)
			iwell+=1;
		if(iwell<nwells){
			if(KE[iwell]<0.0){
				psi=A[iwell]*sin(q[iwell]*r+delta[iwell]);
			}
			else{
				psi=A[iwell]*sinh(q[iwell]*r+delta[iwell]);
			}
		}
		else{
			psi=A[nwells]*exp(q[nwells]*(r-a[nwells-1]));
		}
		printf("%7.3f %d %7.3f\n",r,iwell,psi);
	}	
}

void SolveTanXoverX(double target,double &x){
	double guess,delx,y,dydx,missby=1.0E10;
	int ntry=0;
	if(target>=0.0 && target<=1.0){
		x=0.5;
		y=tanh(x)/x;
		while(fabs(missby)>1.0E-10 && ntry<100){
			ntry+=1;
			dydx=(1.0/x)*pow(cosh(x),-2)-y/x;
			delx=(target-y)/dydx;
			if(fabs(delx)>0.5*fabs(x))
				delx=(0.5*fabs(x))*fabs(delx)/delx;
			x=x+delx;
			y=tanh(x)/x;
			missby=(y-target);
		}
	}
	else if (target>=1.0){
		x=0.25*PI;
		y=tan(x)/x;
		missby=y-target;
		while(fabs(missby)>1.0E-10 && ntry<100){
			ntry+=1;
			dydx=(1.0/x)*pow(cos(x),-2)-y/x;
			delx=(target-y)/dydx;
			if(x+delx>0.5*PI)
				delx=0.75*(0.5*PI-x);
			if(x+delx<0.0)
				delx=-0.25*x;
			x=x+delx;
			
			y=tan(x)/x;
			missby=(y-target);
		}
	}
	else if (target<=0.0){
		x=0.75*PI;
		y=tan(x)/x;
		missby=y-target;
		while(fabs(missby)>1.0E-10 && ntry<100){
			ntry+=1;
			dydx=(1.0/x)*pow(cos(x),-2)-y/x;
			delx=(target-y)/dydx;
			if(x+delx>PI){
				delx=0.75*(PI-x);
			}
			if(x+delx<0.5*PI){
				delx=-0.25*(x-0.5*PI);
			}
			x=x+delx;
			y=tan(x)/x;
			missby=(y-target);
		}
	}
	if(ntry==100){
		printf("Ntry=100!!!\n");
		exit(1);
	}
	printf("----- missby=%g, x=%g, ntry=%d\n",missby,x,ntry);
}
	

