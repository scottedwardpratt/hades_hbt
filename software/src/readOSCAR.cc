#include "hades_hbt/hades_hbt.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

void Chades_hbt_master::ReadOSCAR(){
	int nr_ev=1;
	int firstFile=parmap.getI("FIRST_OSCARFILE",0);
	int lastFile=parmap.getI("LAST_OSCARFILE",10);
	Chades_hbt_cell *cell;
	Chades_hbt_part *tmp_particle;
	string directory=parmap.getS("OSCAR_DIRNAME","oscar_crap");
	
    // opening the files from first -> last:
    ifstream f_in;

    string dust, line; // dummy variable
    double t, x, y, z, mass, p0, px, py, pz, pdg, pid, charge, bim;
    int nrParticlesInEvent;
    int event = 0;
    for(int iFile = firstFile; iFile <= lastFile; iFile++){
      string infile_name = directory +"/"+ to_string(iFile); //create the name of the files to be read
      // open input file
      f_in.open( infile_name );
      if( f_in.fail() )
			CLog::Fatal( "cannot open input file" );
      //getting rid of some lines
      for(int i = 0; i < 3; i++) getline( f_in, line );
        while (!f_in.eof()){
          event = nr_ev;
          for(int i = 0; i < 4; i++) f_in >> dust;
          f_in >> nrParticlesInEvent;
          for(int i = 0; i < nrParticlesInEvent; i++){//reading particles in event loop
            f_in >> t >> x >> y >> z >> mass >> p0 >> px >> py >> pz >> pdg >> pid >> charge;
				
            bool accept; double eff;
            if(pid == PIDA){
              Acceptance(pid,tmp_particle, accept, eff);
				  if(accept){
					  cell_list->FindCell(pid,tmp_particle,cell);
					  tmp_particle=new Chades_hbt_part;
					  tmp_particle->p[0]=p0;
					  tmp_particle->p[1]=px; tmp_particle->p[2]=py; tmp_particle->p[3]=pz;
					  tmp_particle->x[0]=t;
					  tmp_particle->x[1]=x; tmp_particle->x[2]=y; tmp_particle->x[3]=z;
					  cell->partlist_a.push_back(tmp_particle);
				  }
            }
            else if(pid == PIDB ){
              Acceptance(pid, tmp_particle, accept, eff);
              if(accept){
					  cell_list->FindCell(pid,tmp_particle,cell);
					  tmp_particle=new Chades_hbt_part;
					  tmp_particle->p[0]=p0;
					  tmp_particle->p[1]=px; tmp_particle->p[2]=py; tmp_particle->p[3]=pz;
					  tmp_particle->x[0]=t;
					  tmp_particle->x[1]=x; tmp_particle->x[2]=y; tmp_particle->x[3]=z;
					  cell->partlist_b.push_back(tmp_particle);
				  }
            }
          } // end particles loop
          for(int i = 0; i < 6; i++) f_in >> dust;
          f_in >> bim;
          for(int i = 0; i < 2; i++) f_in >> dust;
          event++;
      }
    }//end of loop over files
}
