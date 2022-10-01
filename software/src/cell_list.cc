#include "hades_hbt/hades_hbt.h"

using namespace std;

Chades_hbt_cell_list::Chades_hbt_cell_list(){
	int ix,iy,iz;
	int inx,iny,inz;
	NX=parmap->GetI("NX",10);
	NY=parmap->GetI("NY",10);
	NZ=parmap->GetI("NZ",10);
	lista.resize(NX);
	listb.resize(NX);
	for(ix=0;ix<NX;ix++){
		cell_list[ix].resize[NY];
		for(iy=0;iy<NY;iy++){
			cell_list[ix][iy].resize(NZ);
			for(iz=0;iz<NZ;iz++){
				cell_list[ix][iy][iz]=new Chades_hbt_cell();
			}
		}
	}
	for(ix=0;ix<NX;ix++){
		for(iy=0;iy<NY;iy++){
			for(iz=0;iz<NZ;iz++){
				for(inx=0;inx<3;inx++){
					neighbor[inx].resize(3);
					for(iny=0;iny<3;iny++){
						neighbor[inx][iny].resize(3);
						for(inz=0;inz<3;inz++){
							if((ix+inx-1<NX && ix+inx-1>=0) 
								&&(iy+iny-1><NY && iy+iny-1>=0)
									&&(iz+inz-1<NA && iz+in-1>=0){
										neighbor[inx][iny][inz]=&cell_list[ix+inx-1][iy+iny-1][iz+inz-1];
							}
						}
					}
				}
			}
		}
	}
}

void Chades_hbt_cell_list::FindCell(int pid,Chades_hbt_part *part,Chades_hbt_cell &*cell){
	double px=part->p[1],py=part->p[2],pz=part->p[3];
	double E=part->p[0];
	double mass=sqrt(E*E-px*px-py*py-pz*pz);
	double dpx=px-NX*DPX*0.5;
	int ipx=lrint(floor(dpx/DPX));
	cell=NULL;
	if(ipx>=0 && ipx<NPX){
		double dpy=py-NY*DPY*0.5;
		int ipy=lrint(floor(dpy/DPY));
		if(ipy>=0 && ipy<NPY){
			double dpz=pz-NZ*DPZ*0.5;
			int ipz=lrint(floor(dpz/DPZ));
			if(ipz>=0 && ipz<NPZ){
				cell=cell_list[ipx][ipy][ipz];
			}
		}
	}
	
}
