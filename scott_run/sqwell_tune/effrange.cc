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

