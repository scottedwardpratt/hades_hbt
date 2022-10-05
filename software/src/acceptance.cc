#include "hades_hbt/hades_hbt.h"
using namespace std;

Chades_hbt_acceptance::Chades_hbt_acceptance(){
	thetamin=18.0;
	thetamax=85.0;
	pTmin=0.0;
	pTmax=2500.0;
}

bool Chades_hbt_acceptance::Acceptance(int pid, Chades_hbt_part *part, double &efficiency){
//// === HADES acceptance === ////
  efficiency=0.0;
  
  //momentum:
  int charge = pid / abs(pid);
  double pmag = sqrt(part->p[1]*part->p[1] + part->p[2]*part->p[2] + part->p[3]*part->p[3]);
  if(charge > 0 && pmag < 50.0)
	  return false;
  if(charge < 0 && pmag < 105.0)
	  return false;

  //pT:
  double pT = sqrt(part->p[1]*part->p[1] + part->p[2]*part->p[2]);
  if(pT<pTmin && pT>pTmax)
	  return false;

  //theta:
  double theta = atan(pT/part->p[3])*(180.0/PI);
  if(theta < thetamin || theta > thetamax)
	  return false;

  //rapidity:
  double rapidity = 0.5 * log((part->p[0]+part->p[3])/(part->p[0]-part->p[3]));
  if(rapidity > 2.1)
	  return false;
    

//// === Efficiency === ////
    efficiency=1.0;
	 return true;
}
