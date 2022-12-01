#ifndef __HADES_HBT_H__
#define __HADES_HBT_H__

#include <list>
#include "msu_commonutils/commondefs.h"
#include "msu_commonutils/parametermap.h"
#include "msu_coral/coral.h"
#include "hades_acceptance.h"

using namespace std;

class Chades_hbt_part;
class Chades_hbt_cell;
class Chades_hbt_cell_list;
class Chades_hbt_part;
class Chades_hbt_resinfo;
class Chades_hbt_CFs;
class Chades_hbt_acceptance;

class Chades_hbt_master{
public:
	string parsfilename_prefix;
	Chades_hbt_master(string parmapfilename_prefix);
	int PIDA,PIDB;
	bool HADES_GAUSS;
	CparameterMap parmap;
	void ReadOSCAR_1997();
	void ReadOSCAR_2003();
	CWaveFunction *wf;
	double GetCorrelationWeight(Chades_hbt_part *parta,Chades_hbt_part *partb);
	Chades_hbt_acceptance *acceptance;
	void IncrementCFs(Chades_hbt_part *parta,Chades_hbt_part *partb);
	void CalcCFs();
	void CalcCFs_Gaussian();
	Chades_hbt_cell_list *cell_list;
	Chades_hbt_CFs *cfs;
	int nincrement,nsuccess;
	Crandy *randy;
	char message[200];
};

class Chades_hbt_cell_list{
public:
	Chades_hbt_cell_list(CparameterMap *parmap);
	int NRAPX,NRAPY,NRAPZ;
	double DRAPX,DRAPY,DRAPZ;
	double rapxmax,rapymax,rapzmax;
	double QMAX;
	vector<vector<vector<Chades_hbt_cell *> >> cell;
	void FindCell(Chades_hbt_part *part,Chades_hbt_cell *&cell);
	void Add2List(Chades_hbt_part &parta);
	static Chades_hbt_master *master;
	char message[200];
};

class Chades_hbt_cell{
public:
	Chades_hbt_cell();
	vector<Chades_hbt_part *> partlist_a;
	vector<Chades_hbt_part *> partlist_b;
	vector<vector<vector<Chades_hbt_cell *>>> neighbor;
	//vector<vector<int>> nlist[2][2][2];
};

class Chades_hbt_resinfo{
public:
	Chades_hbt_resinfo();
	int pida,pidb;
	double massa,massb;
};

class Chades_hbt_CFs{
public:
	Chades_hbt_CFs(CparameterMap *parmap);
	int NQINV;
	double DQINV;
	vector<double> C_of_qinv;
	vector<int> denom_of_qinv;
	void PrintC_of_qinv();
	void WriteC_of_qinv();
	void WriteC3D();
	bool XSYM,YSYM,ZSYM;
	
	C3DArray *threed_num,*threed_den;
	double Q3DMAX,DELQ3D;
	int NQ3D;
	static Chades_hbt_master *master;
	
	char message[200];
};

class Chades_hbt_part{
public:
	void Print(){
		double mass=sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]);
		printf("x=(%g,%g,%g,%g)\n",x[0],x[1],x[2],x[3]);
		printf("mass=%g, p=(%g,%g,%g,%g)\n",mass,p[0],p[1],p[2],p[3]);
		printf("psmear=(%g,%g,%g,%g)\n",psmear[0],psmear[1],psmear[2],psmear[3]);
	}
	int pid;
	FourVector x;
	FourVector p; // true momentum
	FourVector psmear; // smeared momentum, due to resolution
	double mass;
	void Setp0(){
		if(abs(pid)==2112)
			mass=NeutronMass;
		else if(abs(pid)==2212)
			mass=ProtonMass;
		else if(abs(pid)==211)
			mass=PionMass;
		p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
		psmear[0]=sqrt(mass*mass+psmear[1]*psmear[1]+psmear[2]*psmear[2]+psmear[3]*psmear[3]);
	}
	char message[200];
};

#endif
