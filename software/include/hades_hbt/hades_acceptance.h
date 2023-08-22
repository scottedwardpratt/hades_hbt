#ifndef __HADES_HBT_ACCEPTANCE_H__#ifndef __HADES_HBT_ACCEPTANCE_H__
#define __HADES_HBT_ACCEPTANCE_H__
#include "hades_hbt.h"

class Chades_hbt_part;
class Chades_hbt_master;

class Chades_hbt_acceptance{
public:
	double pTmax,pTmin,thetamin,thetamax,ymax,ymin;
	char message[CLog::CHARLENGTH];
	static Chades_hbt_master *master;
	CparameterMap *parmap;
	Crandy *randy;
	FourVector ucm;
	Chades_hbt_acceptance(){
		//
	}
	void Init(CparameterMap *parmap);
	
	Chades_hbt_acceptance(CparameterMap *parmap){
		Init(parmap);
	}
	
	virtual bool OneParticleAcceptance(int pid,Chades_hbt_part *part,double &efficiency){
		// this is a dummy function
		return true;
	}
	
	virtual bool TwoParticleAcceptance(Chades_hbt_part *parta,Chades_hbt_part *partb,
	double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
	double &efficiency){
		// this is a dummy function
		return true;
	}
	
	
	virtual void Smear(Chades_hbt_part *part){
		// this is a dummy function
	}
	
	
};

class Chades_hbt_acceptance_nosmear : public Chades_hbt_acceptance{
public:
	Chades_hbt_acceptance_nosmear(CparameterMap *parmap){
		Init(parmap);
	}
	bool OneParticleAcceptance(int pid,Chades_hbt_part *part,double &efficiency);
	bool TwoParticleAcceptance(Chades_hbt_part *parta,Chades_hbt_part *partb,
	double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
	double &efficiency);
	void Smear(Chades_hbt_part *part);
};

class Chades_hbt_acceptance_smear : public Chades_hbt_acceptance{
public:
	Chades_hbt_acceptance_smear(CparameterMap *parmap){
		Init(parmap);
	}
	bool OneParticleAcceptance(int pid,Chades_hbt_part *part,double &efficiency);
	bool TwoParticleAcceptance(Chades_hbt_part *parta,Chades_hbt_part *partb,
	double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
	double &efficiency);
	void Smear(Chades_hbt_part *part);
};

class Chades_hbt_acceptance_smear_maria : public Chades_hbt_acceptance{
public:
	Chades_hbt_acceptance_smear_maria(CparameterMap *parmap){
		Init(parmap);
	}
	bool OneParticleAcceptance(int pid,Chades_hbt_part *part,double &efficiency);
	bool TwoParticleAcceptance(Chades_hbt_part *parta,Chades_hbt_part *partb,
	double qinv,double qout,double qlong,double qside,double deleta,double dely,double delphi,
	double &efficiency);
	void Smear(Chades_hbt_part *part);
};

#endif
