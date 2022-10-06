#include "hades_hbt/hades_hbt.h"

Chades_hbt_CFs::Chades_hbt_CFs(CparameterMap *parmap){
	NQINV=parmap->getI("NQMAX",80);
	DQINV=parmap->getI("DQINV",1.0);
	C_of_qinv.resize(NQINV);
	denom_of_qinv.resize(NQINV);
}

void Chades_hbt_CFs::PrintC_of_qinv(){
	int iq;
	double q;
	printf(" qinv      CF\n");
	for(iq=0;iq<NQINV;iq++){
		q=(0.5+iq)*DQINV;
		printf("%5.1f %8.5g  %d\n",q,C_of_qinv[iq]/denom_of_qinv[iq],denom_of_qinv[iq]);
	}
}