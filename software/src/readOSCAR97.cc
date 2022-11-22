#include "hades_hbt/hades_hbt.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>

using namespace std;

void Chades_hbt_master::ReadOSCAR_1997(){
	Chades_hbt_cell *cell;
	Chades_hbt_part *tmp_particle=new Chades_hbt_part();
	//double taucompare=parmap.getD("OSCAR_TAUCOMPRE",25.0);
	string filenamefilenames=parmap.getS("OSCAR_FILENAME_FILENAMES","oscarfilenames_URQMD.txt");

	double BMIN=parmap.getD("OSCAR_BMIN",4.0);
	double BMAX=parmap.getD("OSCAR_BMAX",4.6);

	// opening the files from first -> last:
	FILE *fptr_in;

	double t, x, y, z, mass, p0, px, py, pz, bim, dumbo;
	int pid,pdg,nparts=0;
	int nrParticlesInEvent,tracknumber=0;
	int nr_event;

	list<string> oscar_filenames;
	list<string>::iterator fiter;
	string filename;
	char dummy[200];
	FILE  *fptr_filenames=fopen(filenamefilenames.c_str(),"r");
	do{
		fscanf(fptr_filenames,"%s",dummy);
		filename=dummy;
		oscar_filenames.push_back(filename);
	}while(!feof(fptr_filenames));
	fclose(fptr_filenames);

	for(fiter=oscar_filenames.begin();fiter!=oscar_filenames.end();++fiter){
		filename=*fiter;
		CLog::Info("Reading "+filename+"\n");
		fptr_in=fopen(filename.c_str(),"r");
		//getting rid of some lines
		for(int i = 0; i < 3; i++){
			fgets(dummy,200,fptr_in);
		}
		do{
			fscanf(fptr_in,"%d %d %lf %lf",&nr_event,&nrParticlesInEvent,&bim,&dumbo);
			fgets(dummy,200,fptr_in);
			if(!feof(fptr_in)){
				for(int i = 1; i <= nrParticlesInEvent; i++){//reading particles in event loop
					fscanf(fptr_in,"%d %d  %lf %lf %lf %lf  %lf  %lf %lf %lf %lf",&tracknumber,&pid,&px,&py,&pz,&p0,&mass,&x,&y,&z,&t);
					fgets(dummy,200,fptr_in);
					if(i!=tracknumber){
						CLog::Info("Warning for file "+filename+", tracknumber suspicious\n");
						CLog::Info("nr_event="+to_string(nr_event)+",tracknumber="+to_string(tracknumber)+"\n");
						Misc::Pause();
					}
					if(bim>=BMIN && bim<=BMAX){
						mass*=1000.0; p0*=1000.0; px*=1000.0; py*=1000.0; pz*=1000.0;
						p0=sqrt(mass*mass+px*px+py*py+pz*pz);
						pdg = pid;
						bool accept; double eff;
						if(pdg == PIDA){
							tmp_particle->p[1]=px;
							tmp_particle->p[2]=py;
							tmp_particle->p[3]=pz;
							tmp_particle->pid=pdg;
							acceptance->Smear(tmp_particle);
							tmp_particle->Setp0();
							tmp_particle->x[1]=x;//-(px/p0)*(t-taucompare);
							tmp_particle->x[2]=y;//-(py/p0)*(t-taucompare);
							tmp_particle->x[3]=z;//-(pz/p0)*(t-taucompare);
							tmp_particle->x[0]=t;//taucompare;

							accept=acceptance->Acceptance(pdg,tmp_particle, eff);
							if(accept){
								cell_list->FindCell(tmp_particle,cell);
								if(cell!=NULL){
									cell->partlist_a.push_back(tmp_particle);
									tmp_particle=new Chades_hbt_part;
									nparts+=1;
								}
							}
						}
						else if(pdg == PIDB ){
							accept=acceptance->Acceptance(pdg, tmp_particle, eff);
							tmp_particle->p[1]=px;
							tmp_particle->p[2]=py;
							tmp_particle->p[3]=pz;
							tmp_particle->pid=pdg;
							acceptance->Smear(tmp_particle);
							tmp_particle->Setp0();
							tmp_particle->x[0]=t;
							tmp_particle->x[1]=x;
							tmp_particle->x[2]=y;
							tmp_particle->x[3]=z;
							if(accept){
								cell_list->FindCell(tmp_particle,cell);
								if(cell!=NULL){
									cell->partlist_b.push_back(tmp_particle);
									tmp_particle=new Chades_hbt_part;
									nparts+=1;
								}
							}
						}
					}
				} // end particles loop
				nrParticlesInEvent = 0;
			}
		}while(!feof(fptr_in));
		fclose(fptr_in);
		CLog::Info("readOSCAR: read in "+to_string(nparts)+" parts\n");
	}//end of loop over files
	delete tmp_particle;
	cout << "Im here, in the end of read OSCAR97" << endl;
}
