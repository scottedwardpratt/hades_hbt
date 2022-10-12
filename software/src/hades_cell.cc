#include "hades_hbt/hades_hbt.h"

Chades_hbt_cell::Chades_hbt_cell(){
	int ix,iy,iz;
	neighbor.resize(3);
	for(ix=0;ix<3;ix++){
		neighbor[ix].resize(3);
		for(iy=0;iy<3;iy++){
			neighbor[ix][iy].resize(3);
			for(iz=0;iz<3;iz++){
				neighbor[ix][iy][iz]=NULL;
			}
		}
	}
}