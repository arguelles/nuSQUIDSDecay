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
  const unsigned int numneu = 3;

  nusquids::marray<double,1> e_nodes = logspace(1.0e0*units.GeV,1.0e4*units.GeV,200);

  nuSQUIDSDecay nusqdec(e_nodes,numneu);

  nusqdec.Set_MixingParametersToDefault();

  std::shared_ptr<EarthAtm> body = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track = std::make_shared<EarthAtm::Track>(acos(-1.));

  //std::shared_ptr<Vacuum> body = std::make_shared<Vacuum>();
  //std::shared_ptr<Vacuum::Track> track = std::make_shared<Vacuum::Track>(12000.*units.km);

  nusqdec.Set_Body(body);
  nusqdec.Set_Track(track);

  marray<double,3> neutrino_state({e_nodes.size(),2,nusqdec.GetNumNeu()});
  std::fill(neutrino_state.begin(),neutrino_state.end(),0);

  //fill this with a real neutrino flux later
  for(size_t ie=0; ie<neutrino_state.extent(0); ie++){
    for(size_t ir=0; ir<neutrino_state.extent(1); ir++){
      for(size_t iflv=0; iflv<neutrino_state.extent(2); iflv++){

        neutrino_state[ie][ir][iflv] = (iflv == mu) ? 1.: 0;
		//Dummy spectrum w/ soft tail.
		neutrino_state[ie][ir][iflv]*= (1.0/e_nodes[ie]);
        //neutrino_state[ie][ir][iflv] = (iflv == mu) ? 1.: 0;
      }
    }
  }


	
  nusqdec.Set_m_phi(1.0e-4);
  nusqdec.Set_m_nu(1.0e-3, 0);
  nusqdec.Set_m_nu(1.0e-2, 1);
  nusqdec.Set_m_nu(1.0e-1, 2);
//  nusqdec.Set_m_nu(1.0, 3);

  nusqdec.Set_initial_state(neutrino_state,flavor);
  //nusqdec.Set_ProgressBar(true);
  //nusqdec.Set_IncoherentInteractions(false);
  nusqdec.Set_IncoherentInteractions(true);

  squids::Const decay_angles;
  std::vector<double> decay_strength(numneu);
  std::fill(decay_strength.begin(),decay_strength.end(),0.);

  decay_angles.SetMixingAngle(0,2,3.1415/4.0);
  decay_angles.SetMixingAngle(1,2,3.1415/4.0);
  decay_angles.SetMixingAngle(0,1,3.1415/4.0);

  decay_strength[1] = 1.0e-14;

  nusqdec.Set_Decay_Matrix(decay_angles,decay_strength);

	std::cout << "About to Evolve!" << std::endl;

  nusqdec.EvolveState();
	std::cout << "Finished Evolution" << std::endl;

  for(size_t ie=0; ie<e_nodes.size(); ie++){
    std::cout << e_nodes[ie]/units.GeV << " ";
    for(size_t flv=0; flv< numneu; flv++){
      std::cout << nusqdec.EvalFlavorAtNode(flv,ie,0) << " " << nusqdec.EvalFlavorAtNode(flv,ie,1) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
