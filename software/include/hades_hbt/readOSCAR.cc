#include "hades_hbt_master.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

void ReadOSCAR(int nr_ev = 1, string directory, int firstFile, int lastFile, varctor<Chades_part*> partList_a, varctor<Chades_part*> partList_b){

    // opening the files from first -> last:
    ifstream f_in;

    string dust, line; // dummy variable
    double t, x, y, z, mass, p0, px, py, pz, pdg, id, charge, bim;
    int nrParticlesInEvent;
    int event = 0;
    for(int iFile = firstFile; iFile <= last; FileiFile++){
      string infile_name = directory + to_string(iFile); //create the name of the files to be read
      // open input file
      f_in.open( infile_name );
      if( f_in.fail() ) error_message( "cannot open input file" );
      //getting rid of some lines
      for(int i = 0; i < 3; i++) getline( f_in, line );
        while (!iFile.eof()){
          event == nr_ev;
          for(int i = 0; i < 4; i++) f_in >> dust;
          f_in >> nrParticlesInEvent;
          for(int i = 0; i < nrParticlesInEvent; i++){//reading particles in event loop
            f_in >> t >> x >> y >> z >> mass >> p0 >> px >> py >> pz >> pdg >> id >> charge;
            Chades_part *tmp_particle(p0, px, py, pz, t, x, y, z);
            bool accept; double eff;
            if(pid == parmap.PIDA){
              acceptance(pdg, tmp_particle, accept, eff);
              if(accept)partList_a.push_back(tmp_particle);
            }
            if(pid == parmap.PIDB){
              acceptance(pdg, tmp_particle, accept, eff);
              if(accept)partList_b.push_back(tmp_particle);
            }
          } // end particles loop
          for(int i = 0; i < 6; i++) f_in >> dust;
          f_in >> bim;
          for(int i = 0; i < 2; i++) f_in >> dust;
          event++;
      }
    }//end of loop over files
}
