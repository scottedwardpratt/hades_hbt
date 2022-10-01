#include "hades_hbt/hades_hbt.h"

using namespace std;

Chades_hbt_cell_list::Chades_hbt_cell_list(CparameterMap *parmap){
	int ix,iy,iz;
	int inx,iny,inz;
	NPX=parmap->getI("NPX",10);
	NPY=parmap->getI("NPY",10);
	NPZ=parmap->getI("NPZ",10);
	cell.resize(NPX);
	for(ix=0;ix<NPX;ix++){
		cell[ix].resize(NPY);
		for(iy=0;iy<NPY;iy++){
			cell[ix][iy].resize(NPZ);
			for(iz=0;iz<NPZ;iz++){
				cell[ix][iy][iz]=new Chades_hbt_cell();
			}
		}
	}
	for(ix=0;ix<NPX;ix++){
		for(iy=0;iy<NPY;iy++){
			for(iz=0;iz<NPZ;iz++){
				cell[ix][iy][iz]->neighbor.resize(3);
				for(inx=0;inx<3;inx++){
					cell[ix][iy][iz]->neighbor[inx].resize(3);
					for(iny=0;iny<3;iny++){
						cell[ix][iy][iz]->neighbor[inx][iny].resize(3);
						for(inz=0;inz<3;inz++){
							if((ix+inx-1<NPX && ix+inx-1>=0) 
								&&(iy+iny-1<NPY && iy+iny-1>=0)
									&&(iz+inz-1<NPZ && iz+inz-1>=0)){
										cell[ix][iy][iz]->neighbor[inx][iny][inz]=cell[ix+inx-1][iy+iny-1][iz+inz-1];
							}
						}
					}
				}
			}
		}
	}
}

void Chades_hbt_cell_list::FindCell(int pid,Chades_hbt_part *part,Chades_hbt_cell *&cellptr){
	double px=part->p[1],py=part->p[2],pz=part->p[3];
	double E=part->p[0];
	double mass=sqrt(E*E-px*px-py*py-pz*pz);
	double dpx=px-NPX*DPX*0.5;
	int ipx=lrint(floor(dpx/DPX));
	cellptr=NULL;
	if(ipx>=0 && ipx<NPX){
		double dpy=py-NPY*DPY*0.5;
		int ipy=lrint(floor(dpy/DPY));
		if(ipy>=0 && ipy<NPY){
			double dpz=pz-NPZ*DPZ*0.5;
			int ipz=lrint(floor(dpz/DPZ));
			if(ipz>=0 && ipz<NPZ){
				cellptr=cell[ipx][ipy][ipz];
			}
		}
	}
	
}
