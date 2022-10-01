#include "hades_hbt/hades_hbt.h"

Chades_hbt_master::Chades_hbt_master(string parsfilename){
	parmap.ReadParsFromFile(parsfilename);	
	string message;
	PIDA=parmap.getI("PIDA",211);
	PIDB=parmap.getI("PIDB",211);
	
	string filename=parmap.getS("LOG_FILENAME","log.txt");
	CLog::Init(filename);
	CLog::INTERACTIVE=false;
	
	if((PIDA==2212 && PIDB==2212) || (PIDA==-2212 && PIDB==-2212)){
		wf=new CWaveFunction_pp_schrod(parsfilename);
	}
	else if((PIDA==211 && PIDB==211) || (PIDA==-211 && PIDB==-211)){
		wf=new CWaveFunction_pipluspiplus_sqwell(parsfilename);
	}
	else{
		CLog::Fatal("Cannot recognize PIDA="+to_string(PIDA)+" or PIDB="+to_string(PIDB)+"\n");
	}
	
	cell_list=new Chades_hbt_cell_list(&parmap);
	
	cfs=new Chades_hbt_CFs(&parmap);
	
}