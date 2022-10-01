#include "hades_hbt/hades_hbt.h"
using namespace std;

void acceptance(int pid, Chades_hbt_part *part, bool &accept, double &efficiency){
//// === HADES acceptance === ////
  accept = true;
  //momentum:
  int charge = pid / abs(pid);
  double momentum = sqrt(part->p[1]*part->p[1] + part->p[2]*part->p[2] + part->p[3]*part->p[3]);
  if(charge > 0 && momentum < 50.0)  accept = false;
  if(charge < 0 && momentum < 105.0) accept = false;

  //pT:
  double pT = sqrt(part->p[1]*part->p[1] + part->p[2]*part->p[2]);
  if(pT > 2500) accept = false;

  //theta:
  double theta = atan(sqrt(part->p[1]*part->p[1] + part->p[2]*part->p[2])/part->p[3]) *(PI/180);
  if(theta < 18 || theta > 85) accept = false;

  //rapidity:
  double rapidity = 0.5 * log((part->p[0]+part->p[3])/(part->p[0]-part->p[3]));
  if(rapidity > 2.1) accept = false;
    
    
//// === Efficiency === ////
    efficiency=1.0;
}
