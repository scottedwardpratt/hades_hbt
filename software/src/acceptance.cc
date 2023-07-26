#include "hades_hbt/hades_hbt.h"
using namespace std;

Chades_hbt_master *Chades_hbt_acceptance::master=NULL;

void Chades_hbt_acceptance::Init(CparameterMap *parmap){
	thetamin=parmap->getD("HADES_ACCEPTANCE_THETAMIN",18.0);
	thetamax=parmap->getD("HADES_ACCEPTANCE_THETAMAX",85.0);
	pTmin=parmap->getD("HADES_ACCEPTANCE_PTMIN",0.0);
	pTmax=parmap->getD("HADES_ACCEPTANCE_PTMAX",2500.0);
	ymin=parmap->getD("HADES_ACCEPTANCE_YMIN",-5.0);
	ymax=parmap->getD("HADES_ACCEPTANCE_YMAX",5.0);
	
	double ebeam=master->parmap.getD("HADES_BEAM_ENERGY_GEV",1.26); // KE of beam in lab frame per nucleon
	int Abeam=master->parmap.getI("HADES_BEAM_A",197);
	int Atarget=master->parmap.getI("HADES_TARGET_A",197);
	double m=0.931,Mtarget,Mbeam,roots,Pbeam,Ebeam;   // roughly account for binding energy
	Mbeam=Abeam*m;
	Mtarget=Atarget*m;
	Ebeam=(m+ebeam)*Abeam;
	Pbeam=sqrt(Ebeam*Ebeam-Mbeam*Mbeam);
	roots=sqrt((Ebeam+Mtarget)*(Ebeam+Mtarget)-Pbeam*Pbeam);
	ucm[0]=(Ebeam+Mtarget)/roots;
	ucm[1]=ucm[2]=0.0;
	ucm[3]=sqrt(ucm[0]*ucm[0]-1.0);
}

bool Chades_hbt_acceptance_nosmear::OneParticleAcceptance(int pid, Chades_hbt_part *part, double &efficiency){
//// === HADES acceptance === ////
  efficiency=0.0;
	FourVector Plab;
	Misc::Boost(ucm,part->psmear,Plab);
  
  //momentum:
  int charge = pid / abs(pid);
	double pmag = sqrt(Plab[1]*Plab[1]+Plab[2]*Plab[2]+Plab[3]*Plab[3]);
	if(charge > 0 && pmag < 50.0)
		return false;
	if(charge < 0 && pmag < 105.0)
		return false;

	//pT:
	double pT = sqrt(Plab[1]*Plab[1]+Plab[2]*Plab[2]);
	if(pT<pTmin && pT>pTmax)
		return false;

	//theta:
	double theta = atan(pT/Plab[3])*(180.0/PI);
	if(theta < thetamin || theta > thetamax)
		return false;

	//rapidity:
	double rapidity = 0.5 * log((Plab[0]+Plab[3])/(Plab[0]-Plab[3]));
	if(rapidity > ymax || rapidity < ymin)
		return false;

	//// === Efficiency === ////
	efficiency=1.0;
	
	return true;
}

void Chades_hbt_acceptance_nosmear::Smear(Chades_hbt_part *part){
	part->psmear[1]=part->p[1];
	part->psmear[2]=part->p[2];
	part->psmear[3]=part->p[3];
	part->Setp0();
}

bool Chades_hbt_acceptance_nosmear::TwoParticleAcceptance(Chades_hbt_part *parta,Chades_hbt_part *partb,
double qout,double qlong,double qside,double deleta,double dely,double delphi,
double &efficiency){
//// === HADES acceptance === ////
	bool acc=true;
  efficiency=1.0;
	if(fabs(dely)<0.001 && fabs(delphi)<0.001){
		efficiency=0.0;
		acc=false;
	}
	
	return acc;
}

bool Chades_hbt_acceptance_smear::OneParticleAcceptance(int pid, Chades_hbt_part *part, double &efficiency){
//// === HADES acceptance === ////
  efficiency=0.0;
	FourVector Plab;
	Misc::Boost(ucm,part->psmear,Plab);
  
  //momentum:
  int charge = pid / abs(pid);
	double pmag = sqrt(Plab[1]*Plab[1]+Plab[2]*Plab[2]+Plab[3]*Plab[3]);
	if(charge > 0 && pmag < 50.0)
		return false;
	if(charge < 0 && pmag < 105.0)
		return false;

	//pT:
	double pT = sqrt(Plab[1]*Plab[1]+Plab[2]*Plab[2]);
	if(pT<pTmin && pT>pTmax)
		return false;

	//theta:
	double theta = atan(pT/Plab[3])*(180.0/PI);
	if(theta < thetamin || theta > thetamax)
		return false;

	//rapidity:
	double rapidity = 0.5 * log((Plab[0]+Plab[3])/(Plab[0]-Plab[3]));
	if(rapidity > ymax || rapidity < ymin)
		return false;

	//// === Efficiency === ////
	efficiency=1.0;
	
	return true;
}

void Chades_hbt_acceptance_smear::Smear(Chades_hbt_part *part){
	double sigmap1,sigmap2,sigmap3;
	FourVector plab;
	part->Setp0();
	Misc::Boost(ucm,part->p,plab);
	
	sigmap1=10.0+fabs(plab[1]*0.02);
	sigmap2=10.0+fabs(plab[2]*0.02);
	sigmap3=10.0+fabs(plab[3]*0.02);
	
	plab[1]=plab[1]+sigmap1*randy->ran_gauss();
	plab[2]=plab[2]+sigmap2*randy->ran_gauss();
	plab[3]=plab[3]+sigmap3*randy->ran_gauss();
	plab[0]=sqrt(part->mass*part->mass+plab[1]*plab[1]+plab[2]*plab[2]+plab[3]*plab[3]);
	
	Misc::BoostToCM(ucm,plab,part->psmear);
	part->Setp0();
}

bool Chades_hbt_acceptance_smear::TwoParticleAcceptance(Chades_hbt_part *parta,Chades_hbt_part *partb,
double qout,double qlong,double qside,double deleta,double dely,double delphi,
double &efficiency){
//// === HADES acceptance === ////
	bool acc=true;
  efficiency=1.0;
	if(fabs(dely)<0.001 && fabs(delphi)<0.001){
		efficiency=0.0;
		acc=false;
	}
	
	return acc;
}

