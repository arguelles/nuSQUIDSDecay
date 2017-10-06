#include <vector>
#include <iostream>
#include <string>
#include <limits>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>
#include <nuSQuIDS/tools.h>
#include "nusquids_decay.h"
//#include "Verosimilitud.h"

//===================MAIN======================================================//
//===================MAIN======================================================//

int main(int argc, char** argv){
  const unsigned int neutrino=0;
  const unsigned int antineutrino=1;
  const unsigned int nu_mu=1;
  bool quiet = false;
  // getting input parameters
  double nu3mass, theta24, lifetime;

  if(argc != 4){
      printf("ERROR:USAGE: program nu3mass [eV] th24 [rad] lifetime [eV^-1]\n");
      exit(0);
  } else {
      nu3mass = atof(argv[1]);
      theta24 = atof(argv[2]);
      lifetime = atof(argv[3]);
  }

  // important paths
  std::string data_path = "/home/carguelles/work/NeutrinoDecay/verosimilitud/data/observed_events.dat";
  std::string flux_path = "/home/carguelles/work/NeutrinoDecay/verosimilitud/data/PolyGonato_QGSJET-II-04.h5";
  std::string effective_area_path = "/home/carguelles/work/NeutrinoDecay/verosimilitud/data/";
  std::string input_flux_path = "/home/carguelles/work/TheSterileSearch/flux_calculation/flux_models/";

  // oscillation physics parameters and nusquids setup
  double dm41sq = nu3mass*nu3mass; // asume m_0 is massless
  const unsigned int numneu = 4;
  const squids::Const units;
  const std::string modelname = "PolyGonato_QGSJET-II-04";
  double m1 = 0.0;
  double m2 = sqrt(7.65e-05);
  double m3 = sqrt(0.0024);
  double m4 = nu3mass;
  double mphi = 0.0;
  std::vector<double> nu_mass{m1,m2,m3,m4};

  // RATE MATRICES

//------------------------//

  double x = m4/m3;
  double f = x/2.+2.+2.*log(x)/x - 2/x/x - 1/x/x/x/2;
  double k = x/2.-2.*log(x)/x - 1/x/x/x/2;
  double r = f/k;

  //std::cout<<"r: "<<r<<std::endl;

  double cpp_scalar_lifetime = lifetime*(1+1/r);
  double cvp_scalar_lifetime = r*cpp_scalar_lifetime;
  double cpp_pseudoscalar_lifetime = 1.0e90;
  double cvp_pseudoscalar_lifetime = 1.0e90;

  //Set chirality-preserving scalar process lifetimes.
  gsl_matrix* cpp_scalar_tau_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_all(cpp_scalar_tau_mat, 1e90); // Set lifetimes to effective stability.
  //Setting for 4 neutrino case with stable nu_1. 
  //  gsl_matrix_set(cpp_scalar_tau_mat,0,1,cpp_scalar_lifetime); //tau_21
  //  gsl_matrix_set(cpp_scalar_tau_mat,0,2,cpp_scalar_lifetime); //tau_31
  //  gsl_matrix_set(cpp_scalar_tau_mat,1,2,cpp_scalar_lifetime); //tau_32
  //  gsl_matrix_set(cpp_scalar_tau_mat,0,3,0.); //tau_41
  //  gsl_matrix_set(cpp_scalar_tau_mat,1,3,0.); //tau_42
  gsl_matrix_set(cpp_scalar_tau_mat,2,3,cpp_scalar_lifetime); //tau_43

  //Set chirality-violating scalar process lifetimes.
  gsl_matrix* cvp_scalar_tau_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_all(cvp_scalar_tau_mat, 1e90); // Set lifetimes to effective stability.
  //Setting for 4 neutrino case with stable nu_1. 
  //  gsl_matrix_set(cvp_scalar_tau_mat,0,1,cvp_scalar_lifetime); //tau_21
  //  gsl_matrix_set(cvp_scalar_tau_mat,0,2,cvp_scalar_lifetime); //tau_31
  //  gsl_matrix_set(cvp_scalar_tau_mat,1,2,cvp_scalar_lifetime); //tau_32
  //  gsl_matrix_set(cvp_scalar_tau_mat,0,3,0.); //tau_41
  //  gsl_matrix_set(cvp_scalar_tau_mat,1,3,0.); //tau_42
  gsl_matrix_set(cvp_scalar_tau_mat,2,3,cvp_scalar_lifetime); //tau_43

  //Set chirality-preserving pseudoscalar process lifetimes.
  gsl_matrix* cpp_pseudoscalar_tau_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_all(cpp_pseudoscalar_tau_mat, 1e90); // Set cpp_pseudoscalar_lifetimes to effective stability.
  //Setting for 4 neutrino case with stable nu_1. 
  //double cpp_pseudoscalar_lifetime = 1.0e60;
  //  gsl_matrix_set(cpp_pseudoscalar_tau_mat,0,1,cpp_pseudoscalar_lifetime); //tau_21
  //  gsl_matrix_set(cpp_pseudoscalar_tau_mat,0,2,cpp_pseudoscalar_lifetime); //tau_31
  //  gsl_matrix_set(cpp_pseudoscalar_tau_mat,1,2,cpp_pseudoscalar_lifetime); //tau_32
  //  gsl_matrix_set(cpp_pseudoscalar_tau_mat,0,3,0.); //tau_41
  //  gsl_matrix_set(cpp_pseudoscalar_tau_mat,1,3,0.); //tau_42
  gsl_matrix_set(cpp_pseudoscalar_tau_mat,2,3,cpp_pseudoscalar_lifetime); //tau_43

  //Set chirality-violating pseudoscalar process lifetimes.
  gsl_matrix* cvp_pseudoscalar_tau_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_all(cvp_pseudoscalar_tau_mat, 1e90); // Set cvp_pseudoscalar_lifetimes to effective stability.
  //Setting for 4 neutrino case with stable nu_1. 
  //double cvp_pseudoscalar_lifetime = 1.0e60;
  //  gsl_matrix_set(cvp_pseudoscalar_tau_mat,0,1,cvp_pseudoscalar_lifetime); //tau_21
  //  gsl_matrix_set(cvp_pseudoscalar_tau_mat,0,2,cvp_pseudoscalar_lifetime); //tau_31
  //  gsl_matrix_set(cvp_pseudoscalar_tau_mat,1,2,cvp_pseudoscalar_lifetime); //tau_32
  //  gsl_matrix_set(cvp_pseudoscalar_tau_mat,0,3,0.); //tau_41
  //  gsl_matrix_set(cvp_pseudoscalar_tau_mat,1,3,0.); //tau_42
  gsl_matrix_set(cvp_pseudoscalar_tau_mat,2,3,cvp_pseudoscalar_lifetime); //tau_43

  gsl_matrix* cpp_scalar_decay_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_zero(cpp_scalar_decay_mat);

	double rate;
	double colrate;
	for (size_t col=0; col<numneu; col++)
	{
		colrate=0;

		for (size_t row=0; row<col; row++)
		{	
			rate = 1.0/gsl_matrix_get(cpp_scalar_tau_mat,row,col);
			gsl_matrix_set(cpp_scalar_decay_mat,row,col,rate);
			colrate+=rate*nu_mass[col];
		}

//		std::cout<< "cpp_scalar: "<< colrate << std::endl;
		gsl_matrix_set(cpp_scalar_decay_mat,col,col,colrate);
	}	


  gsl_matrix* cvp_scalar_decay_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_zero(cvp_scalar_decay_mat);

	for (size_t col=0; col<numneu; col++)
	{
		colrate=0;

		for (size_t row=0; row<col; row++)
		{	
			rate = 1.0/gsl_matrix_get(cvp_scalar_tau_mat,row,col);
			gsl_matrix_set(cvp_scalar_decay_mat,row,col,rate);
			colrate+=rate*nu_mass[col];
		}

//		std::cout<< "cvp_scalar: "<< colrate << std::endl;
		gsl_matrix_set(cvp_scalar_decay_mat,col,col,colrate);
	}	

  gsl_matrix* cpp_pseudoscalar_decay_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_zero(cpp_pseudoscalar_decay_mat);

	for (size_t col=0; col<numneu; col++)
	{
		colrate=0;

		for (size_t row=0; row<col; row++)
		{	
			rate = 1.0/gsl_matrix_get(cpp_pseudoscalar_tau_mat,row,col);
			gsl_matrix_set(cpp_pseudoscalar_decay_mat,row,col,rate);
			colrate+=rate*nu_mass[col];
		}

//		std::cout<< "cpp_pseudo: "<< colrate << std::endl;
		gsl_matrix_set(cpp_pseudoscalar_decay_mat,col,col,colrate);
	}	


  gsl_matrix* cvp_pseudoscalar_decay_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_zero(cvp_pseudoscalar_decay_mat);

	for (size_t col=0; col<numneu; col++)
	{
		colrate=0;

		for (size_t row=0; row<col; row++)
		{	
			rate = 1.0/gsl_matrix_get(cvp_pseudoscalar_tau_mat,row,col);
			gsl_matrix_set(cvp_pseudoscalar_decay_mat,row,col,rate);
			colrate+=rate*nu_mass[col];
		}
//		std::cout<< "cvp_pseudo: "<< colrate << std::endl;
		gsl_matrix_set(cvp_pseudoscalar_decay_mat,col,col,colrate);
	}	



/* OLD RATE CALCULATION
  gsl_matrix* tau_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_all(tau_mat, 1e60); // Set lifetimes to effective stability.
//Setting for 4 neutrino case with stable nu_1.
  gsl_matrix_set(tau_mat,2,3,lifetime); //tau_43

  gsl_matrix* rate_mat = gsl_matrix_alloc(numneu,numneu);
  gsl_matrix_set_zero(rate_mat);

	double rate;
	double colrate;
	for (size_t col=0; col<numneu; col++){
		colrate=0;
		for (size_t row=0; row<col; row++){	
			rate = 1.0/gsl_matrix_get(tau_mat,row,col);
			gsl_matrix_set(rate_mat,row,col,rate);
			colrate+=rate*nu_mass[col];
		}
		gsl_matrix_set(rate_mat,col,col,colrate);
	}
*/



  ////////////////////////////////
  // NUSQUIDS DARK ARTS START HERE
  ////////////////////////////////
  using namespace nusquids;

  if(!quiet)
    std::cout << "Declaring nuSQuIDSDecay atmospheric objects" << std::endl;
  std::shared_ptr<nuSQUIDSAtm<nuSQUIDSDecay>> nusquids_pion = std::make_shared<nuSQUIDSAtm<nuSQUIDSDecay>>(linspace(-1.,0.2,40),logspace(1.e2*units.GeV,1.e6*units.GeV,150),
                                                                                                           numneu,both,true,
                                                                                                           cpp_scalar_decay_mat,cvp_scalar_decay_mat,cpp_pseudoscalar_decay_mat,cvp_pseudoscalar_decay_mat,nu_mass,mphi);

  std::shared_ptr<nuSQUIDSAtm<nuSQUIDSDecay>> nusquids_kaon = std::make_shared<nuSQUIDSAtm<nuSQUIDSDecay>>(linspace(-1.,0.2,40),logspace(1.e2*units.GeV,1.e6*units.GeV,150),
                                                                                                           numneu,both,true,
                                                                                                           cpp_scalar_decay_mat,cvp_scalar_decay_mat,cpp_pseudoscalar_decay_mat,cvp_pseudoscalar_decay_mat,nu_mass,mphi);



  nusquids_kaon->Set_TauRegeneration(true);
  nusquids_pion->Set_TauRegeneration(true);

  // set mixing angles and masses
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

  // setup integration settings
  double error = 1.0e-15;
  //double error = 1.0e-19;
  nusquids_pion->Set_GSL_step(gsl_odeiv2_step_rkf45);
  nusquids_pion->Set_rel_error(error);
  nusquids_pion->Set_abs_error(error);

  nusquids_kaon->Set_GSL_step(gsl_odeiv2_step_rkf45);
  nusquids_kaon->Set_rel_error(error);
  nusquids_kaon->Set_abs_error(error);

  if(!quiet)
    std::cout << "Setting up the initial fluxes for the nuSQuIDSDecay objects." << std::endl;

  marray<double,4> inistate_kaon {nusquids_kaon->GetNumCos(),nusquids_kaon->GetNumE(),2,numneu};
  std::fill(inistate_kaon.begin(),inistate_kaon.end(),0);

  // read file
  marray<double,2> input_kaon_flux = quickread(input_flux_path + "/" + "initial_kaon_atmopheric_" + modelname + ".dat");

  marray<double,1> cos_range = nusquids_kaon->GetCosthRange();
  marray<double,1> e_range = nusquids_kaon->GetERange();
  for ( int ci = 0 ; ci < nusquids_kaon->GetNumCos(); ci++){
    for ( int ei = 0 ; ei < nusquids_kaon->GetNumE(); ei++){
      double enu = e_range[ei]/units.GeV;
      double cth = cos_range[ci];

      inistate_kaon[ci][ei][0][0] = 0.;
      inistate_kaon[ci][ei][0][1] = input_kaon_flux[ci*e_range.size() + ei][2];
      inistate_kaon[ci][ei][0][2] = 0.;
      inistate_kaon[ci][ei][0][3] = 0.;

      inistate_kaon[ci][ei][1][0] = 0.;
      inistate_kaon[ci][ei][1][1] = input_kaon_flux[ci*e_range.size() + ei][3];
      inistate_kaon[ci][ei][1][2] = 0.;
      inistate_kaon[ci][ei][1][3] = 0.;
    }
  }

  nusquids_kaon->Set_initial_state(inistate_kaon,flavor);
//  nusquids_kaon->WriteStateHDF5("/data/user/mmoulai/osc/initial_kaon_"+ std::to_string(nu3mass) + "_" + std::to_string(theta24) + "_" + std::to_string(lifetime) +".hdf5");

  std::ofstream ioutput("/home/mmoulai/osc_rdo/initial_" + std::to_string(nu3mass) + "_" + std::to_string(theta24) + "_" + std::to_string(lifetime) + ".dat");
  for(double costh : nusquids_kaon->GetCosthRange()){
    for(double enu : nusquids_kaon->GetERange()){
    //for(double log10E = 3.0; log10E < 6.0; log10E +=  0.001 ){
        // const double enu = pow(10.0,log10E);
        ioutput << costh << " ";
        ioutput << enu << " ";
        ioutput << nusquids_kaon->EvalFlavor(nu_mu,costh,enu,neutrino) << " ";
        ioutput << nusquids_kaon->EvalFlavor(nu_mu,costh,enu,antineutrino) << " ";
        ioutput << nusquids_kaon->EvalFlavor(nu_mu,costh,enu,neutrino) << " ";
        ioutput << nusquids_kaon->EvalFlavor(nu_mu,costh,enu,antineutrino) << " ";
        ioutput << std::endl;
    }
  }
  ioutput.close();

  if(!quiet)
    std::cout << "Evolving the kaon fluxes." << std::endl;
  nusquids_kaon->EvolveState();
//  nusquids_kaon->WriteStateHDF5("/data/user/mmoulai/osc/final_kaon_" + std::to_string(nu3mass) + "_" + std::to_string(theta24) + "_" + std::to_string(lifetime) + ".hdf5");

 // std::ofstream output("final_" + std::to_string(nu3mass) + "_" + std::to_string(theta24) + "_" + std::to_string(lifetime) + ".dat");
  std::ofstream output("/home/mmoulai/osc_rdo/final_" + std::to_string(nu3mass) + "_" + std::to_string(theta24) + "_" + std::to_string(lifetime) + ".dat");
  
  for(double costh : nusquids_kaon->GetCosthRange()){
    for(double enu : nusquids_kaon->GetERange()){
    //for(double log10E = 3.0; log10E < 6.0; log10E +=  0.001 ){
        // const double enu = pow(10.0,log10E);
        output << costh << " ";
        output << enu << " ";
        output << nusquids_kaon->EvalFlavor(nu_mu,costh,enu,neutrino) << " ";
        output << nusquids_kaon->EvalFlavor(nu_mu,costh,enu,antineutrino) << " ";
        output << nusquids_kaon->EvalFlavor(nu_mu,costh,enu,neutrino) << " ";
        output << nusquids_kaon->EvalFlavor(nu_mu,costh,enu,antineutrino) << " ";
        output << std::endl;
    }
  }
  output.close();

//  marray<double,4> inistate_pion {nusquids_pion->GetNumCos(),nusquids_pion->GetNumE(),2,numneu};
//  std::fill(inistate_pion.begin(),inistate_pion.end(),0);
//
//  // read file
//  marray<double,2> input_pion_flux = quickread(input_flux_path + "/" + "initial_pion_atmopheric_" + modelname + ".dat");
//
//  for ( int ci = 0 ; ci < nusquids_pion->GetNumCos(); ci++){
//    for ( int ei = 0 ; ei < nusquids_pion->GetNumE(); ei++){
//      double enu = e_range[ei]/units.GeV;
//      double cth = cos_range[ci];
//
//      inistate_pion[ci][ei][0][0] = 0.;
//      inistate_pion[ci][ei][0][1] = input_pion_flux[ci*e_range.size() + ei][2];
//      inistate_pion[ci][ei][0][2] = 0.;
//      inistate_pion[ci][ei][0][3] = 0.;
//
//      inistate_pion[ci][ei][1][0] = 0.;
//      inistate_pion[ci][ei][1][1] = input_pion_flux[ci*e_range.size() + ei][3];
//      inistate_pion[ci][ei][1][2] = 0.;
//      inistate_pion[ci][ei][1][3] = 0.;
//    }
//  }
//
//  nusquids_pion->Set_initial_state(inistate_pion,flavor);
//  nusquids_pion->WriteStateHDF5("/data/user/mmoulai/osc/initial_pion"+ std::to_string(nu3mass) + "_" + std::to_string(theta24) + "_" + std::to_string(lifetime) +".hdf5");
//  if(!quiet)
//    std::cout << "Evolving the pion fluxes." << std::endl;
//  nusquids_pion->EvolveState();
//  nusquids_pion->WriteStateHDF5("/data/user/mmoulai/osc/final_pion_" + std::to_string(nu3mass) + "_" + std::to_string(theta24) + "_" + std::to_string(lifetime) + ".hdf5");
//
  return 0;
}
