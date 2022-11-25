#include "msu_commonutils/commonutils.h"
#include "msu_coral/coral.h"
#include "hades_hbt/hades_hbt.h"

using namespace std;

int main(int argc,char *argv[]){
	string parsfilename;
	if(argc!=2){
		printf("Usage hades_hbt_fromGAUSS parameter_file_name\n");
		exit(1);
	}
	else{
		parsfilename=argv[1];
	}
	Chades_hbt_master *hades_hbt_master=new Chades_hbt_master(parsfilename);
	hades_hbt_master->CalcCFs_Gaussian(2.5,3.5,5.0);
	hades_hbt_master->cfs->PrintC_of_qinv();
	hades_hbt_master->cfs->WriteC_of_qinv("qinv_output_gauss/qinv.txt");
	//hades_hbt_master->cfs->WriteC3D("threed_output_gauss");
	return 0;
}
