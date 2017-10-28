/*========================="Couplings" Example==========================//
Below is an example implementation of the NuSQuIDSDecay class.
An example neutrino flux (kaon/pion) is read in over a specified
range and binning in cos(zenith angle) and energy. NuSQuIDSDecay
is then used to evolve this flux through the earth to the South 
Pole. Both the initial and final fluxes are written to text files
which can be used to produce oscillograms. 
	The neutrinos here are majorana, and incoherent interactions,
tau regeneration, and decay regeneration effects are all being
simulated. We consider a simplified decay scenario where 
all mass states except m_4 are stable, the only
decay channel is from m_4 to m_3, and the only non-zero mixing
angle between the light mass states and m_4 is theta_24. The phi
mass is always assumed to be zero. All m_4->m_3 decay processes 
({CPP,CVP}x{SCALAR,PSEUDOSCALAR}) are allowed, but they are 
computed internally by NuSQuIDSDecay as functions of the lagrangian 
coupling matrix g_ij which we supply to the constructor.
//==========================================================================*/

#include <vector>
#include <iostream>
#include <string>
#include <limits>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>
#include <nuSQuIDS/tools.h>
#include "nusquids_decay.h"

using namespace nusquids;

//Convenience function to write initial and final fluxes to a text file.
void WriteFlux(std::shared_ptr<nuSQUIDSAtm<nuSQUIDSDecay>> nusquids, std::string fname){
	enum {NU_E,NU_MU,NU_TAU};
	enum {NEUTRINO,ANTINEUTRINO};
	std::ofstream ioutput("../output/" + fname + ".dat");
	for(double costh : nusquids->GetCosthRange()){
		for(double enu : nusquids->GetERange()){
			ioutput << costh << " ";
			ioutput << enu << " ";
			ioutput << nusquids->EvalFlavor(NU_MU,costh,enu,NEUTRINO) << " ";
			ioutput << nusquids->EvalFlavor(NU_MU,costh,enu,ANTINEUTRINO) << " ";
			ioutput << nusquids->EvalFlavor(NU_MU,costh,enu,NEUTRINO) << " ";
			ioutput << nusquids->EvalFlavor(NU_MU,costh,enu,ANTINEUTRINO) << " ";
			ioutput << std::endl;
		}
	}
	ioutput.close();
}

//Convenience function to read a flux file into an marray input for nuSQuIDS.
void ReadFlux(std::shared_ptr<nuSQUIDSAtm<nuSQUIDSDecay>> nusquids, marray<double,4>& inistate, std::string type, 
				std::string input_flux_path, std::string modelname, double GeV){

	std::fill(inistate.begin(),inistate.end(),0);
	// read file
	marray<double,2> input_flux = quickread(input_flux_path + "/" + "initial_"+ type + "_atmopheric_" + modelname + ".dat");

	marray<double,1> cos_range = nusquids->GetCosthRange();
	marray<double,1> e_range = nusquids->GetERange();
	for ( int ci = 0 ; ci < nusquids->GetNumCos(); ci++){
		for ( int ei = 0 ; ei < nusquids->GetNumE(); ei++){
			double enu = e_range[ei]/GeV;
			double cth = cos_range[ci];

			inistate[ci][ei][0][0] = 0.;
			inistate[ci][ei][0][1] = input_flux[ci*e_range.size() + ei][2];
			inistate[ci][ei][0][2] = 0.;
			inistate[ci][ei][0][3] = 0.;

			inistate[ci][ei][1][0] = 0.;
			inistate[ci][ei][1][1] = input_flux[ci*e_range.size() + ei][3];
			inistate[ci][ei][1][2] = 0.;
			inistate[ci][ei][1][3] = 0.;
		}
	}
}

//====================================MAIN=========================================//

int main(int argc, char** argv){
	bool oscillogram = true;
	bool quiet = false;
	// getting input parameters
	double nu3mass, theta24;
	nu3mass = 1.0; //Set the mass of the sterile neutrino (eV)
	theta24 = 1.0; //Set the mixing angle between sterile and tau flavors.

	//Set scalar and pseudoscalar couplings.
	double scalar_coupling=1.0;
	double pseudoscalar_coupling=2.0;

	//Toggle majorana/dirac, incoherent interactions, and decay regeneration.
	bool iinteraction=true;
	bool decay_regen=true;
	bool majorana=true;

	//Path for input fluxes
	std::string input_flux_path = "../fluxes";
	
	//Flux model
	const std::string modelname = "PolyGonato_QGSJET-II-04";

	// oscillation physics parameters and nusquids setup
	double dm41sq = nu3mass*nu3mass; // assume m_0 is massless
	const unsigned int numneu = 4;
	const squids::Const units;
	double m1 = 0.0;
	double m2 = sqrt(7.65e-05);
	double m3 = sqrt(0.0024);
	double m4 = nu3mass;
	std::vector<double> nu_mass{m1,m2,m3,m4};

	enum{SCALAR, PSEUDOSCALAR};
	
	//Allocate memory for rate matrices 
	gsl_matrix* couplings[2];
	for (size_t s=0; s<2; s++){				 
		couplings[s] = gsl_matrix_alloc(numneu,numneu);
		gsl_matrix_set_zero(couplings[s]);
	}

	//Set chirality-preserving scalar coupling (43 nonzero).
	gsl_matrix_set(couplings[SCALAR],3,2,scalar_coupling); //g_43
	//Set chirality-preserving pseudoscalar coupling (43 nonzero).
	gsl_matrix_set(couplings[PSEUDOSCALAR],3,2,pseudoscalar_coupling); //g_43

	//Declare NuSQuIDSDecay objects. They are declared within a NuSQuIDSAtm wrapper to incorporate atmospheric simulation.
	//Here, we use the partial rate constructor of NuSQuIDSDecay. One object is created for the kaon flux component, and
	//one for the pion flux component.
	//The first two arguments (linspaces) define ranges of cos(zenith) and energy over which to simulate, respectively.
	if(!quiet){std::cout << "Declaring nuSQuIDSDecay atmospheric objects" << std::endl;}
	std::shared_ptr<nuSQUIDSAtm<nuSQUIDSDecay>> nusquids_pion = std::make_shared<nuSQUIDSAtm<nuSQUIDSDecay>>(linspace(-1.,0.2,40),
																logspace(1.e2*units.GeV,1.e6*units.GeV,150),numneu,both,iinteraction,
																decay_regen,majorana,nu_mass,couplings);

	std::shared_ptr<nuSQUIDSAtm<nuSQUIDSDecay>> nusquids_kaon = std::make_shared<nuSQUIDSAtm<nuSQUIDSDecay>>(linspace(-1.,0.2,40),
																logspace(1.e2*units.GeV,1.e6*units.GeV,150),numneu,both,iinteraction,
																decay_regen,majorana,nu_mass,couplings);

	//Include tau regeneration in simulation.
	nusquids_kaon->Set_TauRegeneration(true);
	nusquids_pion->Set_TauRegeneration(true);

	//Set mixing angles and masses.
	nusquids_kaon->Set_MixingAngle(0,1,0.563942);
	nusquids_kaon->Set_MixingAngle(0,2,0.154085);
	nusquids_kaon->Set_MixingAngle(1,2,0.785398); 
	nusquids_kaon->Set_MixingAngle(0,3,0.0);
	nusquids_kaon->Set_MixingAngle(1,3,theta24);
	nusquids_kaon->Set_MixingAngle(2,3,0.0);

	nusquids_kaon->Set_SquareMassDifference(1,7.65e-05);
	nusquids_kaon->Set_SquareMassDifference(2,0.00247);
	nusquids_kaon->Set_SquareMassDifference(3,dm41sq);
	nusquids_kaon->Set_CPPhase(0,2,0.0);
	nusquids_kaon->Set_CPPhase(0,3,0.0);
	nusquids_kaon->Set_CPPhase(1,3,0.0);

	nusquids_pion->Set_MixingAngle(0,1,0.563942);
	nusquids_pion->Set_MixingAngle(0,2,0.154085);
	nusquids_pion->Set_MixingAngle(1,2,0.785398);
	nusquids_pion->Set_MixingAngle(0,3,0.0);
	nusquids_pion->Set_MixingAngle(1,3,theta24);
	nusquids_pion->Set_MixingAngle(2,3,0.0);

	nusquids_pion->Set_SquareMassDifference(1,7.65e-05);
	nusquids_pion->Set_SquareMassDifference(2,0.00247);
	nusquids_pion->Set_SquareMassDifference(3,dm41sq);
	nusquids_pion->Set_CPPhase(0,2,0.0);
	nusquids_pion->Set_CPPhase(0,3,0.0);
	nusquids_pion->Set_CPPhase(1,3,0.0);

	//Setup integration settings
	double error = 1.0e-15;
	nusquids_pion->Set_GSL_step(gsl_odeiv2_step_rkf45);
	nusquids_pion->Set_rel_error(error);
	nusquids_pion->Set_abs_error(error);

	nusquids_kaon->Set_GSL_step(gsl_odeiv2_step_rkf45);
	nusquids_kaon->Set_rel_error(error);
	nusquids_kaon->Set_abs_error(error);

	if(!quiet)
		std::cout << "Setting up the initial fluxes for the nuSQuIDSDecay objects." << std::endl;

	//Read kaon flux and initialize nusquids object with it.
	marray<double,4> inistate_kaon {nusquids_kaon->GetNumCos(),nusquids_kaon->GetNumE(),2,numneu};
	ReadFlux(nusquids_kaon,inistate_kaon,std::string("kaon"),input_flux_path,modelname,units.GeV);
	nusquids_kaon->Set_initial_state(inistate_kaon,flavor);
	//Write initial flux to text file.
	if(oscillogram){WriteFlux(nusquids_kaon, std::string("kaon_initial"));}
	//Evolve flux through the earth. 
	if(!quiet){std::cout << "Evolving the kaon fluxes." << std::endl;}
	nusquids_kaon->EvolveState();
	//Write final flux to text file.
	if(oscillogram){WriteFlux(nusquids_kaon, std::string("kaon_final"));}

	//Read pion flux and initialize nusquids object with it.
	marray<double,4> inistate_pion {nusquids_pion->GetNumCos(),nusquids_pion->GetNumE(),2,numneu};
	ReadFlux(nusquids_pion,inistate_pion,std::string("pion"),input_flux_path,modelname,units.GeV);
	nusquids_pion->Set_initial_state(inistate_pion,flavor);
	//Write initial flux to text file.
	if(oscillogram){WriteFlux(nusquids_pion, std::string("pion_initial"));}
	//Evolve flux through the earth. 
	if(!quiet){std::cout << "Evolving the pion fluxes." << std::endl;}
	nusquids_pion->EvolveState();
	//Write final flux to text file.
	if(oscillogram){WriteFlux(nusquids_pion, std::string("pion_final"));}

	//Free memory for couplings 
	for (size_t s=0; s<2; s++){				 
		gsl_matrix_free(couplings[s]);
	}
	return 0;
}
