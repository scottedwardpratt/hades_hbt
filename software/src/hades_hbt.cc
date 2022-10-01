#include "hades_hbt/hades_hbt.h"

Chades_hbt_master::Chades_hbt_master(CparameterMap *parmapset){
	parmap=parmapset;	
	string message;
	PIDA=parmap->getI("PIDA",211);
	PIDB=parmap->getI("PIDB",211);
	
	string filename=parmap->getS("LOG_FILENAME","log.txt");
	log.init(parmap->getS(filename));
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
	
	cell_list::log=&log;
	Chades_hbt_acceptance::log=&log;
	Chades_hbt_cell::log=&log;
	
}