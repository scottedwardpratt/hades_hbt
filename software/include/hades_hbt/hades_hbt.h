#ifndef __HADES_HBT_H__
#define __HADES_HBT_H__

#include <list>
#include "msu_commonutils/commondefs.h"
#include "msu_commonutils/parametermap.h"
#include "msu_coral/coral.h"

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
	Chades_hbt_master(string parsfilename);
	int PIDA,PIDB;
	CparameterMap parmap;
	void ReadOSCAR();
	CWaveFunction *wf;
	double GetCorrelationWeight(Chades_hbt_part *parta,Chades_hbt_part *partb);
	Chades_hbt_acceptance *acceptance;
	void IncrementCFs(Chades_hbt_part *parta,Chades_hbt_part *partb);
	void CalcCFs();
	void CalcCFs_Gaussian(double Rx,double Ry,double Rz);
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

class Chades_hbt_acceptance{
public:
	Chades_hbt_acceptance();
	bool Acceptance(int pid,Chades_hbt_part *part,double &efficiency);
	double pTmax,pTmin,thetamin,thetamax;
	char message[200];
};

class Chades_hbt_CFs{
public:
	Chades_hbt_CFs(CparameterMap *parmap);
	int NQINV;
	double DQINV;
	vector<double> C_of_qinv;
	vector<int> denom_of_qinv;
	void PrintC_of_qinv();
	void WriteC_of_qinv(string filename);
	void WriteC3D(string dirname);
	bool XSYM,YSYM,ZSYM;
	
	C3DArray *threed_num,*threed_den;
	double Q3DMAX,DELQ3D;
	int NQ3D;
	
	char message[200];
};

class Chades_hbt_part{
public:
	void Print(){
		double mass=sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]);
		printf("x=(%g,%g,%g,%g)\n",x[0],x[1],x[2],x[3]);
		printf("mass=%g, p=(%g,%g,%g,%g)\n",p[0],mass,p[1],p[2],p[3]);
	}
	FourVector x;
	FourVector p;
	char message[200];
};

#endif