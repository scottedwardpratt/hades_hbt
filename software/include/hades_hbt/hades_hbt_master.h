#ifndef __HADES_HBT_H__
#define __HADES_HBT_H__

#include "msu_commonutils/commondefs.h"
#include "msu_commonutils/parametermap.h"
#include "msu_coral/coral.h"

using namespace std;

class Chades_hbt_part;
class Chades_hbt_cell;
class Chades_hbt_part;
class Chades_hbt_resinfo;

class Chades_hades_hbt_master{
public:
	Chades_hhbt_master();
	CparameterMap *parmap;
	void ReadOSCAR();
	CWaveFunction *wf;
	double GetCorrelationWeight(Chades_hbt_part *parta,Chades_hbt_part *partb);
	void Acceptance(int pid,Chades_hbt_part *part);
	void IncrementCFs(Chades_hbt_part *parta,Chades_hbt_part *partb,double weight,double  efficiency);
	Chades_hbt_cell_list *cell_list;
	Chades_hbt_CFs *cfs;
};

class Chades_hbt_cell_list{
public:
	Chades_hbt_cell_list();
	int NX,NY,NZ;
	double DPX,DPY,DPZ;
	vector<vector<vector<Chades_hbt_cell *>>> *cell;
	void FindCell(Chades_hbt_part *part &*cell);
	void Add2List(Chades_hbt_part *parta);
	void IncrementCFs(Chades_hbt_part *parta,Chades_hbt_*partb); // Calulates relative momentum, phi^2, etc, then increments CFs accordingly.
};

class Chades_hbt_cell{
public:
	vector<Chades_part *> partlist;
	vector<vector<vector<Chades_hbt_cell>>> *neighbor;
};

class Chades_hbt_resinfo{
public:
	Chades_hbt_resinfo();
	int pida,pidb;
	double massa,massb;
};

class Chades_hbt_acceptance{
public:
	void acceptance(Chades_hbt_resinfo *resinfo,Chades_hbt_part,bool &accept,double &efficiency);
};

class Chades_hbt_CFs{
	int NQinv;
	double DQINV;
	vector<double> C_of_qinv;
};
