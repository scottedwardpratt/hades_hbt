#include "hades_hbt/hades_hbt.h"

Chades_hbt_master::Chades_hbt_master(string parsfilename){
	printf("parsfilename=%s\n",parsfilename.c_str());
	parmap.ReadParsFromFile(parsfilename);	
	string message;
	PIDA=parmap.getI("PIDA",211);
	PIDB=parmap.getI("PIDB",211);
	
	string filename=parmap.getS("LOG_FILENAME","log.txt");
	CLog::INTERACTIVE=false;
	CLog::Init(filename);
	
	if((PIDA==2212 && PIDB==2212) || (PIDA==-2212 && PIDB==-2212)){
		wf=new CWaveFunction_pp_schrod(parsfilename);
	}
	else if((PIDA==211 && PIDB==211) || (PIDA==-211 && PIDB==-211)){
		wf=new CWaveFunction_pipluspiplus_sqwell(parsfilename);
	}
	else{
		CLog::Fatal("Cannot recognize PIDA="+to_string(PIDA)+" or PIDB="+to_string(PIDB)+"\n");
	}
	
	Chades_hbt_cell_list::master=this;
	cell_list=new Chades_hbt_cell_list(&parmap);
	
	cfs=new Chades_hbt_CFs(&parmap);
	
	acceptance=new Chades_hbt_acceptance();
	
}

void Chades_hbt_master::CalcCFs(){
	int icx,icy,icz,jcx,jcy,jcz,inx,iny,inz,ia,ib,na,nb;
	int inxmin=1,inymin,inzmin,ibmin;
	int nincrement=0;
	if(PIDA!=PIDB)
		inxmin=0;
	Chades_hbt_cell *cella,*cellb;
	Chades_hbt_part *parta,*partb;
	for(icx=0;icx<cell_list->NRAPX;icx++){
		for(icy=0;icy<cell_list->NRAPY;icy++){
			for(icz=0;icz<cell_list->NRAPZ;icz++){
				cella=cell_list->cell[icx][icy][icz];
				na=cella->partlist_a.size();
				if(na>0)
				for(ia=0;ia<na;ia++){
					parta=cella->partlist_a[ia];
					
					for(inx=inxmin;inx<3;inx++){
						inymin=0;
						if(PIDA==PIDB && inx==1)
							inymin=1;
						for(iny=inymin;iny<3;iny++){
							inzmin=0;
							if(PIDA==PIDB && inx==1 && iny==1)
								inzmin=1;
							for(inz=inzmin;inz<3;inz++){
								if(cella->neighbor[inx][iny][inz]!=NULL){
									jcx=icx+inx;
									jcy=icy+iny;
									jcz=icz+inz;
									cellb=cell_list->cell[jcx][jcy][jcz];
									ibmin=0;
									if(PIDA==PIDB && inx==1 && iny==1 && inz==1)
										ibmin=ia+1;
									if(PIDA==PIDB){
										nb=cellb->partlist_a.size();
									}
									else
										nb=cellb->partlist_b.size();
									for(ib=ibmin;ib<nb;ib++){
										if(PIDA==PIDB)
											partb=cellb->partlist_a[ib];
										else
											partb=cellb->partlist_b[ib];
										nincrement+=1;
										IncrementCFs(parta,partb);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	printf("nincrement=%d\n",nincrement);
}

void Chades_hbt_master::IncrementCFs(Chades_hbt_part *parta,Chades_hbt_part *partb){
	double q,r,ctheta,weight;
	int iq;
	wf->getqrctheta(parta->p,parta->x,partb->p,partb->x,&q,&r,&ctheta);
	weight=wf->GetPsiSquared(q,r,ctheta);
	
	if(q<cfs->DQINV*cfs->NQINV){
		iq=lrint(floor(q/cfs->DQINV));
		cfs->C_of_qinv[iq]+=weight;
		cfs->denom_of_qinv[iq]+=1;		
	}
}