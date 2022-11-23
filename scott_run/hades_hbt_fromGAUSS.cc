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
	printf("check a\n");
	hades_hbt_master->CalcCFs_Gaussian(2.5,3.5,5.0);
	printf("check b\n");
	hades_hbt_master->cfs->PrintC_of_qinv();
	printf("check c\n");
	hades_hbt_master->cfs->WriteC_of_qinv("qinv_output_gauss/qinv.txt");
	printf("check d\n");
	//hades_hbt_master->cfs->WriteC3D("threed_output_gauss");
	return 0;
}
