#include "msu_commonutils/commonutils.h"
#include "msu_coral/coral.h"
#include "hades_hbt/hades_hbt.h"

using namespace std;

int main(){
	string parsfilename="parameters/coralpars.dat";
	Chades_hbt_master *hades_hbt_master=new Chades_hbt_master(parsfilename);
	hades_hbt_master->ReadOSCAR_1997();
	//hades_hbt_master->ReadOSCAR_2003();
	hades_hbt_master->CalcCFs();
	hades_hbt_master->cfs->PrintC_of_qinv();
	hades_hbt_master->cfs->WriteC_of_qinv("qinv_output/qinv.txt");
	hades_hbt_master->cfs->WriteC3D("threed_output");
	return 0;
}

