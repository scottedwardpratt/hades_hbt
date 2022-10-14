#include "hades_hbt/hades_hbt.h"

Chades_hbt_master::Chades_hbt_master(string parsfilename){
	parmap.ReadParsFromFile(parsfilename);	
	string filename=parmap.getS("LOG_FILENAME","log.txt");
	CLog::Init(filename);
	
	CLog::Info("parsfilename="+parsfilename+"\n");
	PIDA=parmap.getI("PIDA",211);
	PIDB=parmap.getI("PIDB",211);
	if(PIDA==PIDB)
		parmap.set("XYZSYM","true");
	
	
	
	
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
	randy=new Crandy(-12345);
	
}

void Chades_hbt_master::CalcCFs(){
	int icx,icy,icz,inx,iny,inz,ia,ib,na,nb;
	int inxmin,inymin,inzmin,ibmin;
	int natot=0;
	nincrement=nsuccess=0;
	inxmin=0;
	if(PIDA==PIDB)
		inxmin=1;
	Chades_hbt_cell *cella,*cellb;
	Chades_hbt_part *parta,*partb;
	printf("calculating CFs\n");
	for(icx=0;icx<cell_list->NRAPX;icx++){
		for(icy=0;icy<cell_list->NRAPY;icy++){
			for(icz=0;icz<cell_list->NRAPZ;icz++){
				cella=cell_list->cell[icx][icy][icz];
				na=cella->partlist_a.size();
				natot+=na;
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
								cellb=cella->neighbor[inx][iny][inz];
								if(cellb!=NULL){
									printf("asize=%lu bsize=%lu\n",cella->partlist_a.size(),cellb->partlist_a.size());
									ibmin=0;
									if(PIDA==PIDB && inx==1 && iny==1 && inz==1){
										ibmin=ia+1;
									}
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
	CLog::Info("nincrement="+to_string(nincrement)+", nincrement/npairs_tot="
		+to_string(2.0*double(nincrement)/(double(natot)*double(natot-1)))+"\n");//, nwf/nincrement="+to_string(double(nsuccess)/double(nincrement));
}

void Chades_hbt_master::CalcCFs_Gaussian(double Rx,double Ry,double Rz){
	double x,y,z,qx,qy,qz,q,r,ctheta,weight;
	int imc,NMC=parmap.getI("NMC_GAUSSIAN",1000);
	int iq,iqx,iqy,iqz,isx,isy,isz,nsx=2,nsy=2,nsz=2;
	if(cfs->XSYM)
		nsx=1;
	if(cfs->YSYM)
		nsy=1;
	if(cfs->ZSYM)
		nsz=1;
	for(iqx=0;iqx<cfs->NQ3D;iqx++){
		for(isx=0;isx<nsx;isx++){
			for(iqy=0;iqy<cfs->NQ3D;iqy++){
				for(isy=0;isy<nsy;isy++){
					for(iqz=0;iqz<cfs->NQ3D;iqz++){
						for(isz=0;isz<nsz;isz++){
							for(imc=0;imc<NMC;imc++){
								qx=cfs->DELQ3D*(iqx+randy->ran());
								qy=cfs->DELQ3D*(iqy+randy->ran());
								qz=cfs->DELQ3D*(iqz+randy->ran());
								if(isx>0)
									qx=-qx;
								if(isy>0)
									qy=-qy;
								if(isz>0)
									qz=-qz;
								q=sqrt(qx*qx+qy*qy+qz*qz);
								weight=1.0;
								if(q<cell_list->QMAX){
									if(q<cfs->DQINV*cfs->NQINV){
										x=Rx*randy->ran_gauss();
										y=Ry*randy->ran_gauss();
										z=Rz*randy->ran_gauss();
										r=sqrt(x*x+y*y+z*z);
										ctheta=(qx*x+qy*y+qz*z)/(q*r);
										weight=wf->GetPsiSquared(q,r,ctheta);
										iq=lrint(floor(q/cfs->DQINV));
										cfs->C_of_qinv[iq]+=weight;
										cfs->denom_of_qinv[iq]+=1;
										if(fabs(qx)<cfs->Q3DMAX && fabs(qy)<cfs->Q3DMAX && fabs(qz)<cfs->Q3DMAX){
											nsuccess+=1;
											cfs->threed_num->IncrementElement(qx,qy,qz,weight);
											cfs->threed_den->IncrementElement(qx,qy,qz,1.0);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void Chades_hbt_master::IncrementCFs(Chades_hbt_part *parta,Chades_hbt_part *partb){
	double q,r,ctheta,weight=1.0;
	int iq;
	wf->getqrctheta(parta->p,parta->x,partb->p,partb->x,q,r,ctheta);
	
	if(r>1.0E-8){
	
		if(q<cell_list->QMAX){
			nsuccess+=1;		
			if(q<cfs->DQINV*cfs->NQINV){
				weight=wf->GetPsiSquared(q,r,ctheta);
				if(weight!=weight){
					parta->Print();
					partb->Print();
					CLog::Fatal("weight=Nan\n");
				}
				iq=lrint(floor(q/cfs->DQINV));
				if(iq<int(cfs->C_of_qinv.size())){
					cfs->C_of_qinv[iq]+=weight;
					cfs->denom_of_qinv[iq]+=1;
				}
			}
		}
	}
	
	double qout,qlong,qside,deleta,dely,delphi;
	Misc::outsidelong(parta->p,partb->p,q,qout,qside,qlong,deleta,dely,delphi);
	if(fabs(qout)<cfs->Q3DMAX && fabs(qside)<cfs->Q3DMAX && fabs(qlong)<cfs->Q3DMAX){
		cfs->threed_num->IncrementElement(qout,qlong,qside,weight);
		cfs->threed_den->IncrementElement(qout,qlong,qside,1.0);
	}
	
}