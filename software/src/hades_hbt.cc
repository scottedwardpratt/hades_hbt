#include "hades_hbt.h"

Chades_hades_hbt_master::Chades_hades_hbt_master(CparamterMap *parmap){
	string message;
	PIDA=parmap->GetI("PIDA",211);
	PIDB=parmap->GetI("PIDB",211);
	
	log.init(parmap.getS("LOG_FILENAME","log.txt"));
	log.interactive=false;
	
	if((PIDA==2212 && PIDB==2212) || (PIDA==-2212 && PIDB==-2212)){
		wf=new CWaveFunction_pp_schrod(parmap);
	}
	else if((PIDA==211 && PIDB==211) || (PIDA==-211 && PIDB==-211)){
		wf=new CWaveFunction_pipluspiplus_sqwell(parmap);
	}
	else{
		log.fatal("Cannot recognize PIDA="+str(PIDA)+" or PIDB="+str(PIDB)+"\n");
	}
	
	cell_list=new Chades_hbt_cell_list(parmap);
	
	cfs=new Chades_hbt_CFs(parmap);
	
}