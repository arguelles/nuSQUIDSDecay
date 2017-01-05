#include <vector>
#include <iostream>
#include <string>
#include <limits>
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

  nusquids::marray<double,1> e_nodes = logspace(1.0e2*units.GeV,1.0e5*units.GeV,100);

  nuSQUIDSDecay nusqdec(e_nodes,numneu);

	// setup integration settings
	double tolerance=1.0e-12;
	nusqdec.Set_rel_error(tolerance);
	nusqdec.Set_abs_error(tolerance);

  std::shared_ptr<EarthAtm> body = std::make_shared<EarthAtm>();
  std::shared_ptr<EarthAtm::Track> track = std::make_shared<EarthAtm::Track>(acos(-1.));

  //std::shared_ptr<Vacuum> body = std::make_shared<Vacuum>();
  //std::shared_ptr<Vacuum::Track> track = std::make_shared<Vacuum::Track>(2.0*6371.0*units.km);

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
		//neutrino_state[ie][ir][iflv]*= (1.0/e_nodes[ie]);
        //neutrino_state[ie][ir][iflv] = (iflv == mu) ? 1.: 0;
      }
    }
  }


  nusqdec.Set_SquareMassDifference(3,1.0);  //dm^2_41
	nusqdec.SetPhiMass(0.0);
	nusqdec.SetNeutrinoMasses(0.0);
	
  nusqdec.Set_MixingParametersToDefault();

  // mixing angles
  //nusqdec.Set_MixingAngle(0,1,0.563942);
  //nusqdec.Set_MixingAngle(0,2,0.154085);
  //nusqdec.Set_MixingAngle(1,2,0.785398);
  nusqdec.Set_MixingAngle(0,3,0.785398);
  nusqdec.Set_MixingAngle(1,3,0.785398);
  nusqdec.Set_MixingAngle(2,3,0.785398);

  nusqdec.Set_initial_state(neutrino_state,flavor);

	//------------------------//
	//   Physics Switches	    //
	//------------------------//

  //nusqdec.Set_ProgressBar(true);
  nusqdec.SetIncoherentInteractions(false);
	nusqdec.SetDecayRegeneration(true);
	
	//------------------------//


  gsl_matrix* TauMatrix = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_all(TauMatrix, 1e60); // Set lifetimes to effective stability.

	//Setting for 4 neutrino case with stable nu_1. 
  double lifetime = 1.0e2;

  gsl_matrix_set(TauMatrix,0,3,lifetime); //tau_41
  gsl_matrix_set(TauMatrix,1,3,lifetime); //tau_42
  gsl_matrix_set(TauMatrix,2,3,lifetime); //tau_43

  nusqdec.SetDecayMatrix(TauMatrix);

  nusqdec.EvolveState();

  for(size_t ie=0; ie<e_nodes.size(); ie++){
    std::cout << e_nodes[ie]/units.GeV << " ";
    for(size_t flv=0; flv< numneu; flv++){
      std::cout << nusqdec.EvalFlavorAtNode(flv,ie,0) << " ";
    }
    for(size_t flv=0; flv< numneu; flv++){
	  if (flv==numneu-1)
	  {
        std::cout << nusqdec.EvalFlavorAtNode(flv,ie,1) << std::endl;
	  }
	  else
	  {
        std::cout << nusqdec.EvalFlavorAtNode(flv,ie,1) << " ";
	  }
    }
  }

  gsl_matrix_free(TauMatrix);  
  return 0;
}
