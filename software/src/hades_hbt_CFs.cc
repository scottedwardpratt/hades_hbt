#include "hades_hbt/hades_hbt.h"

Chades_hbt_CFs::Chades_hbt_CFs(CparameterMap *parmap){
	NQINV=parmap->getI("NQINV",80);
	DQINV=parmap->getI("DQINV",1.0);
	C_of_qinv.resize(NQINV);
}