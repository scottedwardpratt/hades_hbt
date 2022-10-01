#ifndef __HADES_HBT_H__
#define __HADES_HBT_H__

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

class Chades_hbt_master{
public:
	Chades_hbt_master(string parsfilename);
	int PIDA,PIDB;
	CparameterMap parmap;
	void ReadOSCAR();
	CWaveFunction *wf;
	double GetCorrelationWeight(Chades_hbt_part *parta,Chades_hbt_part *partb);
	void Acceptance(int pid,Chades_hbt_part *part,bool &accept,double &efficiency);
	void IncrementCFs(Chades_hbt_part *parta,Chades_hbt_part *partb,double weight,double  efficiency);
	Chades_hbt_cell_list *cell_list;
	Chades_hbt_CFs *cfs;
};

class Chades_hbt_cell_list{
public:
	Chades_hbt_cell_list(CparameterMap *parmap);
	int NPX,NPY,NPZ;
	double DPX,DPY,DPZ;
	vector<vector<vector<Chades_hbt_cell *> >> cell;
	void FindCell(int pid,Chades_hbt_part *part,Chades_hbt_cell *&cell);
	void Add2List(Chades_hbt_part &parta);
	void IncrementCFs(Chades_hbt_part &parta,Chades_hbt_part &partb); // Calulates relative momentum, phi^2, etc, then increments CFs accordingly.
};

class Chades_hbt_cell{
public:
	Chades_hbt_cell();
	vector<Chades_hbt_part *> partlist_a;
	vector<Chades_hbt_part *> partlist_b;
	vector<vector<vector<Chades_hbt_cell *>>> neighbor;
};

class Chades_hbt_resinfo{
public:
	Chades_hbt_resinfo();
	int pida,pidb;
	double massa,massb;
};

class Chades_hbt_acceptance{
public:
	void acceptance(int pid,Chades_hbt_part *part,bool &accept,double &efficiency);
};

class Chades_hbt_CFs{
public:
	Chades_hbt_CFs(CparameterMap *parmap);
	int NQINV;
	double DQINV;
	vector<double> C_of_qinv;
};

class Chades_hbt_part{
public:
	double x[4];
	double p[4];
};

#endif