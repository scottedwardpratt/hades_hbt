#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;
double HBARC=197.3269718;
double PI=4.0*atan(1.0);

int main(){
	const int NK=50;
	const double delk=2.0;
	double k,q,delta;
	double V0,R,a,q0;
	double mproton=938.27,md,mu,Ep;
	md=2.0*mproton-2.224;
	mu=mproton*md/(mproton+md);

	printf("Enter R: ");
	scanf("%lf",&R);
	q0=0.5*PI*HBARC/R;
	V0=0.5*q0*q0/mu;
	printf("V0 for bound state=%g, Enter V0: ",V0);
	scanf("%lf",&V0);
		
	
	/*for(V0=21.0;V0<22;V0+=0.05){
		q0=sqrt(2.0*mu*V0)/HBARC;
		a=R-tan(q0*R)/q0;
		printf("V0=%g, a=%g\n",V0,a);
	}
	exit(1);*/
	//R=1.0; V0=81.3;
	//R=0.5; V0=316.0;
	//R=2.0; V0=21.6;
	
	
	q0=sqrt(2.0*mu*V0)/HBARC;
	a=R-tan(q0*R)/q0;
	printf("a=%g\n",a);
	
	for(k=delk;k<=(0.5+NK)*delk;k+=delk){
		q=sqrt(q0*q0*HBARC*HBARC+k*k);
		delta=-(k*R/HBARC)+atan((k/q)*tan(q*R/HBARC));
		delta=delta*180.0/PI;
		Ep=0.5*k*k/mproton;
		printf("%8.4f %8.5f %8.2f\n",k,Ep,delta);		
	}
	

	
	return 0;
}

void SquareWell_Init(int ell,int nchannels){
  // To set up the wave functions and phase shifts
  CGSLMatrix_Complex *cmatrix;
  double q,mu_coulomb,E;
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
  complex<double> **M,*Y;
  complex<double> x1,x2;
  double F2,G2,F2prime,G2prime,qsquared,r;
  int i,j,iq,ir,ichannel;
  mu=m1*m2/(m1+m2);
	
	for(ichannel=0;ichannel<nchannels;ichannel++){
		if(nwells[ichannel]==1){
			for(iq=0;iq<nqmax;iq++){
				q=GetQ(iq);
				E=sqrt(q*q+m1*m1)+sqrt(q*q+m2*m2);
				mu_coulomb=0.25*(E-pow(m1*m1-m2*m2,2)/pow(E,3));
				eta0=q1q2*mu_coulomb*ALPHA/q;
				
				qsquared=q*q-2.0*mu*V0[ichannel][0];
				if(qsquared>0) q1=sqrt(qsquared);
				else q1=ci*sqrt(abs(qsquared));
				x1=q1*a[ichannel][0]/HBARC;
				eta1=q1q2*mu_coulomb*ALPHA/q1;
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x1,eta1,&F1,&G1,&F1prime,&G1prime);
				x2=q*a[ichannel][0]/HBARC;
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2,eta0,&F2,&G2,&F2prime,&G2prime);     
				beta=(abs(q1)/q)*F1prime/F1;
				delta[ichannel][iq]=-atan2(beta*F2-F2prime,beta*G2-G2prime);
				A[ichannel][iq][0]=0.5*(exp(-2.0*ci*delta[ichannel][iq])
					*(F2+ci*G2)+(F2-ci*G2))/F1;
				A[ichannel][iq][1]=exp(-2.0*ci*delta[ichannel][iq]);
				
			}
		}
		else if(nwells[ichannel]==2){
			cmatrix=new CGSLMatrix_Complex(4);
			Y=new complex<double>[4];
			M=new complex<double> *[4];
			for(i=0;i<4;i++) M[i]=new complex<double>[4];
			
			for(iq=0;iq<nqmax;iq++){
				q=GetQ(iq);
				E=sqrt(q*q+m1*m1)+sqrt(q*q+m2*m2);
				mu_coulomb=0.25*(E-pow(m1*m1-m2*m2,2)/pow(E,3));
				
				q1=sqrt(abs(q*q-2.0*mu*V0[ichannel][0]));
				if(q*q-2.0*mu*V0[ichannel][0]<0.0) q1=ci*q1;
				q2=sqrt(abs(q*q-2.0*mu*V0[ichannel][1]));
				if(q*q-2.0*mu*V0[ichannel][1]<0.0) q2=ci*q2;
				x1b=a[ichannel][0]*q1/HBARC;
				x2a=a[ichannel][0]*q2/HBARC;
				x2b=a[ichannel][1]*q2/HBARC;
				x=a[ichannel][1]*q/HBARC;
				eta1=q1q2*mu_coulomb*ALPHA/q1;
				eta2=q1q2*mu_coulomb*ALPHA/q2;
				eta0=q1q2*mu_coulomb*ALPHA/q;
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x1b,eta1,&F1b,&G1b,&F1bprime,&G1bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2a,eta2,&F2a,&G2a,&F2aprime,&G2aprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2b,eta2,&F2b,&G2b,&F2bprime,&G2bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x,eta0,&F,&G,&Fprime,&Gprime);
				
				for(i=0;i<4;i++){
					Y[i]=0.0;
					for(j=0;j<4;j++) M[i][j]=0.0;
				}
				M[0][0]=F1b;              M[0][1]=-F2a;              M[0][2]=-G2a;
				M[1][0]=abs(q1)*F1bprime; M[1][1]=-abs(q2)*F2aprime; M[1][2]=-abs(q2)*G2aprime;
				M[2][1]=F2b;              M[2][2]=G2b;               M[2][3]=-0.5*(F+ci*G);
				M[3][1]=abs(q2)*F2bprime; M[3][2]=abs(q2)*G2bprime;  M[3][3]=-0.5*q*(Fprime+ci*Gprime);
				
				Y[2]=0.5*(F-ci*G); Y[3]=0.5*q*(Fprime-ci*Gprime);
				cmatrix->SolveLinearEqs(Y,M,A[ichannel][iq]);
				
				delta[ichannel][iq]=-0.5*atan2(imag(A[ichannel][iq][3]),real(A[ichannel][iq][3]));
				if(delta[ichannel][iq]<0.0) delta[ichannel][iq]=delta[ichannel][iq]+PI;
				
			}
			delete(cmatrix);
			delete [] Y;
			for(i=0;i<4;i++) delete [] M[i];
			delete [] M;
		}
		else if(nwells[ichannel]==3){
			cmatrix=new CGSLMatrix_Complex(6);
			Y=new complex<double>[6];
			M=new complex<double> *[6];
			for(i=0;i<6;i++) M[i]=new complex<double>[6];
			
			for(iq=0;iq<nqmax;iq++){
				q=GetQ(iq);
				E=sqrt(q*q+m1*m1)+sqrt(q*q+m2*m2);
				mu_coulomb=0.25*(E-pow(m1*m1-m2*m2,2)/pow(E,3));
				
				q1=sqrt(abs(q*q-2.0*mu*V0[ichannel][0]));
				if(q*q-2.0*mu*V0[ichannel][0]<0.0) q1=ci*q1;
				q2=sqrt(abs(q*q-2.0*mu*V0[ichannel][1]));
				if(q*q-2.0*mu*V0[ichannel][1]<0.0) q2=ci*q2;
				q3=sqrt(abs(q*q-2.0*mu*V0[ichannel][2]));
				if(q*q-2.0*mu*V0[ichannel][2]<0.0) q3=ci*q3;
				x1b=a[ichannel][0]*q1/HBARC;
				x2a=a[ichannel][0]*q2/HBARC;
				x2b=a[ichannel][1]*q2/HBARC;
				x3a=a[ichannel][1]*q3/HBARC;
				x3b=a[ichannel][2]*q3/HBARC;
				x=a[ichannel][2]*q/HBARC;
				eta1=q1q2*mu_coulomb*ALPHA/q1;
				eta2=q1q2*mu_coulomb*ALPHA/q2;
				eta3=q1q2*mu_coulomb*ALPHA/q3;
				eta0=q1q2*mu_coulomb*ALPHA/q;
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x1b,eta1,&F1b,&G1b,&F1bprime,&G1bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2a,eta2,&F2a,&G2a,&F2aprime,&G2aprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x2b,eta2,&F2b,&G2b,&F2bprime,&G2bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x3a,eta3,&F3a,&G3a,&F3aprime,&G3aprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x3b,eta3,&F3b,&G3b,&F3bprime,&G3bprime);
				CoulWave::GetFGprime_ComplexQ(ell[ichannel],x,eta0,&F,&G,&Fprime,&Gprime);
				
				for(i=0;i<6;i++){
					Y[i]=0.0;
					for(j=0;j<6;j++) M[i][j]=0.0;
				}
				M[0][0]=F1b;              M[0][1]=-F2a;              M[0][2]=-G2a;
				M[1][0]=abs(q1)*F1bprime; M[1][1]=-abs(q2)*F2aprime; M[1][2]=-abs(q2)*G2aprime;
				M[2][1]=F2b;              M[2][2]=G2b;               M[2][3]=-F3a;               M[2][4]=-G3a;
				M[3][1]=abs(q2)*F2bprime; M[3][2]=abs(q2)*G2bprime;  M[3][3]=-abs(q3)*F3aprime;  M[3][4]=-abs(q3)*G3aprime;
				M[4][3]=F3b;              M[4][4]=G3b;               M[4][5]=-0.5*(F+ci*G);
				M[5][3]=abs(q3)*F3bprime; M[5][4]=abs(q3)*G3bprime;  M[5][5]=-0.5*q*(Fprime+ci*Gprime);
				
				Y[4]=0.5*(F-ci*G); Y[5]=0.5*q*(Fprime-ci*Gprime);
				cmatrix->SolveLinearEqs(Y,M,A[ichannel][iq]);	
				delta[ichannel][iq]=-0.5*atan2(imag(A[ichannel][iq][5]),real(A[ichannel][iq][5]));
				if(delta[ichannel][iq]<0.0) delta[ichannel][iq]=delta[ichannel][iq]+PI;
			}
			delete(cmatrix);
			delete [] Y;
			for(i=0;i<6;i++) delete [] M[i];
			delete [] M;
		}
		else{
			sprintf(message,"nwells[%d] not equal to 1, 2 or 3??? =%d\n",ichannel,nwells[ichannel]);
			CLog::Fatal(message);
		}
	}
	
	for(iq=0;iq<nqmax;iq++){
		for(ichannel=0;ichannel<nchannels;ichannel++) DelPhiArray[iq][0][ichannel]=0.0;
		for(ir=1;ir<=DelPhiArray_NRMAX;ir++){
			r=ir*DelPhiArray_DELR;
			SquareWell_CalcDelPhi(iq,r,DelPhiArray[iq][ir]);
		}
	}
	sprintf(message,"FINISHED INITIALIZATION OF WAVEFUNCTIONS FOR PARTIAL WAVES\n");
	CLog::Info(message);
	
}

