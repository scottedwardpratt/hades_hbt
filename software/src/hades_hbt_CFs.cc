#include "hades_hbt/hades_hbt.h"

Chades_hbt_CFs::Chades_hbt_CFs(CparameterMap *parmap){
	NQinv=parmap->getI("NQinv",80);
	DQinv=parmap->getI("DQinv",1.0);
	C_of_qinv.resize(NQinv);
}