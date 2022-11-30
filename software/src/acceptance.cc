#include "hades_hbt/hades_hbt.h"
using namespace std;

void Chades_hbt_acceptance::Init(CparameterMap *parmap){
	thetamin=parmap->getD("HADES_ACCEPTANCE_THETAMIN",18.0);
	thetamax=parmap->getD("HADES_ACCEPTANCE_THETAMAX",85.0);
	pTmin=parmap->getD("HADES_ACCEPTANCE_PTMIN",0.0);
	pTmax=parmap->getD("HADES_ACCEPTANCE_PTMAX",2500.0);
	ymin=parmap->getD("HADES_ACCEPTANCE_YMIN",-5.0);
	ymax=parmap->getD("HADES_ACCEPTANCE_YMAX",5.0);
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
	if(rapidity > ymax || rapidity < ymin)
		return false;

	//// === Efficiency === ////
	efficiency=1.0;
	return true;
}

void Chades_hbt_acceptance_smear::Smear(Chades_hbt_part *part){
	double sigmap1,sigmap2,sigmap3;
	sigmap1=10.0+fabs(part->p[1]*0.02);
	sigmap2=10.0+fabs(part->p[2]*0.02);
	sigmap3=10.0+fabs(part->p[3]*0.02);
	part->psmear[1]=part->p[1]+sigmap1*randy->ran_gauss();
	part->psmear[2]=part->p[2]+sigmap2*randy->ran_gauss();
	part->psmear[3]=part->p[3]+sigmap3*randy->ran_gauss();
	part->Setp0();
}

void Chades_hbt_acceptance_nosmear::Smear(Chades_hbt_part *part){
	part->psmear[1]=part->p[1];
	part->psmear[2]=part->p[2];
	part->psmear[3]=part->p[3];
	part->Setp0();
}

