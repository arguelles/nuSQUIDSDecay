#include <vector>
#include <iostream>
#include <string>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>
#include <nuSQuIDS/tools.h>
#include "nusquids_decay.h"

using namespace nusquids;

int main(int argc, char* argv[])
{

  squids::Const units;
  const unsigned int e = 0;
  const unsigned int mu = 1;
  const unsigned int tau = 2;
  const unsigned int numneu = 4;

  nusquids::marray<double,1> e_nodes = logspace(1.0e0*units.GeV,1.0e4*units.GeV,100);

  nuSQUIDSDecay nusqdec(e_nodes,numneu);

  std::shared_ptr<Vacuum> vacuum = std::make_shared<Vacuum>();
  std::shared_ptr<Vacuum::Track> vacuum_track = std::make_shared<Vacuum::Track>(12000.*units.km);
  nusqdec.Set_Body(vacuum);
  nusqdec.Set_Track(vacuum_track);

  marray<double,3> neutrino_state({e_nodes.size(),2,nusqdec.GetNumNeu()});
  std::fill(neutrino_state.begin(),neutrino_state.end(),0);

  //fill this with a real neutrino flux later
  for(size_t ie=0; ie<neutrino_state.extent(0); ie++){
    for(size_t ir=0; ir<neutrino_state.extent(1); ir++){
      for(size_t iflv=0; iflv<neutrino_state.extent(2); iflv++){
        neutrino_state[ie][ie][iflv] = (iflv == mu) ? 1.: 0;
        neutrino_state[ie][ie][iflv] = (iflv == mu) ? 1.: 0;
      }
    }
  }

  nusqdec.Set_initial_state(neutrino_state,flavor);
  //nusqdec.Set_ProgressBar(true);
  nusqdec.Set_IncoherentInteractions(false);

  squids::Const decay_angles;
  std::vector<double> decay_strength {0.,0.,0.,0.};
  nusqdec.Set_Decay_Matrix(decay_angles,decay_strength);

  nusqdec.EvolveState();

  for(size_t ie=0; ie<e_nodes.size(); ie++){
    std::cout << e_nodes[ie]/units.GeV << " ";
    std::cout << nusqdec.EvalFlavorAtNode(e,ie,0) << " " << nusqdec.EvalFlavorAtNode(e,ie,1) << " ";
    std::cout << nusqdec.EvalFlavorAtNode(mu,ie,0) << " " << nusqdec.EvalFlavorAtNode(mu,ie,1) << " ";
    std::cout << nusqdec.EvalFlavorAtNode(tau,ie,0) << " " << nusqdec.EvalFlavorAtNode(tau,ie,1) << " ";
    std::cout << std::endl;
  }

  return 0;
}
