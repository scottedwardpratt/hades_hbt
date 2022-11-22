#include "hades_hbt/hades_hbt.h"
using namespace std;

Chades_hbt_acceptance::Chades_hbt_acceptance(){
	thetamin=18.0;
	thetamax=85.0;
	pTmin=0.0;
	pTmax=2500.0;
}

bool Chades_hbt_acceptance::Acceptance(int pid, Chades_hbt_part *part, double &efficiency){
//// === HADES acceptance === ////
  efficiency=0.0;
  
  //momentum:
  int charge = pid / abs(pid);
	double pmag = sqrt(part->psmear[1]*part->psmear[1] + part->psmear[2]*part->psmear[2] + part->psmear[3]*part->psmear[3]);
	if(charge > 0 && pmag < 50.0)
		return false;
	if(charge < 0 && pmag < 105.0)
		return false;

	//pT:
	double pT = sqrt(part->psmear[1]*part->psmear[1] + part->psmear[2]*part->psmear[2]);
	if(pT<pTmin && pT>pTmax)
		return false;

	//theta:
	double theta = atan(pT/part->psmear[3])*(180.0/PI);
	if(theta < thetamin || theta > thetamax)
		return false;

	//rapidity:
	double rapidity = 0.5 * log((part->psmear[0]+part->psmear[3])/(part->psmear[0]-part->psmear[3]));
	if(rapidity > 2.1)
		return false;

	//// === Efficiency === ////
	efficiency=1.0;
	return true;
}

void Chades_hbt_acceptance::Smear(Chades_hbt_part *part){
	double sigmap1,sigmap2,sigmap3;
	sigmap1=30.0+fabs(part->p[1]*0.05);
	sigmap2=30.0+fabs(part->p[2]*0.05);
	sigmap3=30.0+fabs(part->p[3]*0.05);
	part->psmear[1]=part->p[1]+sigmap1*randy->ran_gauss();
	part->psmear[2]=part->p[2]+sigmap2*randy->ran_gauss();
	part->psmear[3]=part->p[3]+sigmap3*randy->ran_gauss();
	part->Setp0();
}
