#include "hades_hbt/hades_hbt.h"
using namespace std;
using namespace NMSUPratt;

Chades_hbt_master *Chades_hbt_CFs::master=NULL;

Chades_hbt_CFs::Chades_hbt_CFs(CparameterMap *parmap){
	NQINV=parmap->getI("HADES_NQINV",80);
	DQINV=parmap->getI("HADES_DELQINV",1.0);
	C_of_qinv.resize(NQINV);
	denom_of_qinv.resize(NQINV);
	
	// Instantiate 3D array
	NQ3D=parmap->getI("HADES_NQ3DARRAY",20);
	DELQ3D=parmap->getD("HADES_DELQ3D",4);
	Q3DMAX=NQ3D*DELQ3D;
	XSYM=YSYM=ZSYM=parmap->getB("XYZSYM",false);
	XSYM=parmap->getB("XSYM",XSYM);
	YSYM=parmap->getB("YSYM",YSYM);
	ZSYM=parmap->getB("ZSYM",ZSYM);
	
	
	threed_num=new C3DArray(NQ3D,DELQ3D,XSYM,YSYM,ZSYM);
	threed_den=new C3DArray(NQ3D,DELQ3D,XSYM,YSYM,ZSYM);
		
}

void Chades_hbt_CFs::PrintC_of_qinv(){
	int iq;
	double q;
	snprintf(message,CLog::CHARLENGTH," qinv      CF\n");
	CLog::Info(message);
	for(iq=0;iq<NQINV;iq++){
		q=(0.5+iq)*DQINV;
		snprintf(message,CLog::CHARLENGTH,"%5.1f %8.5g  %d\n",q,C_of_qinv[iq]/denom_of_qinv[iq],denom_of_qinv[iq]);
		CLog::Info(message);
	}
}

void Chades_hbt_CFs::WriteC_of_qinv(){
	string dirname="results/"+master->parsfilename_prefix;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/"+"qinv.txt";
	int iq;
	double q;
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr," qinv      CF\n");
	for(iq=0;iq<NQINV;iq++){
		q=(0.5+iq)*DQINV;
		fprintf(fptr,"%5.1f %8.5g  %d\n",q,C_of_qinv[iq]/denom_of_qinv[iq],denom_of_qinv[iq]);
	}
	fclose(fptr);
}

void Chades_hbt_CFs::WriteC3D(){
	string dirname="results/"+master->parsfilename_prefix;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	
	threed_num->DivideByArray(threed_den);
	threed_num->WriteArray(dirname);
	dirname=dirname+"_count";
	command="mkdir -p "+dirname;
	system(command.c_str());
	threed_den->WriteArray(dirname);

}