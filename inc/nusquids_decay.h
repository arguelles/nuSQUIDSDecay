#ifndef nusquidslv_H
#define nusquidslv_H

/*
Copyright © 2016 Alexander Moss, Marjon Moulai, Carlos Arguelles, and Janet
Conrad.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the “Software”), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

For more information please email:

Alexander Moss (zander@mit.edu)
Marjon Moulai (marjon@mit.edu)
Carlos Arguelles (caad@mit.edu)
Janet Conrad (conrad@mit.edu)

*/

#include <vector>
#include <iostream>
#include <nuSQuIDS/nuSQuIDS.h>
#include "exCross.h"


namespace nusquids {

//! Implements neutrino decay using the nusquids package.
class nuSQUIDSDecay : public nuSQUIDS {
private:
	//! A switch for incoherent interactions.
	/*!
	This does not control decay regeneration, although it is
	an incoherent process. That is controlled by SetDecayRegeneration().
	*/
  bool iincoherent_int = false;

	//! A flag for setting the decay rate matrix.
  bool DecayParametersSet = false;

	//! A flag for setting neutrino masses.
  bool NeutrinoMassesSet = false;
	
	//! A flag for setting neutrino masses.
  bool PhiMassSet = false;

	//decay parameters: FIXME clean up! Remove from HDF5?
  squids::Const decay_parameters;
  std::vector<double> decay_strength;

	//! The Gamma matrix from the paper.
  squids::SU_vector DT;
	//! The Gamma matrix, in evolving basis.
  std::vector<squids::SU_vector> DT_evol;

	//! The decay rate matrix.
	/*!
	Entries are calculated in SetDecayMatrix.
	The matrix is upper-triangular. The strictly upper triangular 
	elements are the decay rates from the mass state indexed by the
	matrix column to the mass state indexed by the matrix row.
	The diagonal elements are total decay rates from the mass state
	indexed by the column to all daughters. 	
	*/
  gsl_matrix* rate_mat;

	//! A matrix containing parent-frame momenta of daughter neutrinos.
	/*!
	Parent-daughter pair indexed by column-row pair, as in the case of
	the rate matrix.
	*/ 
	gsl_matrix* pstar_mat;

	//! A vector of neutrino masses.
  std::vector<double> NeutrinoMasses;
	
	//! Mass of the phi particle.
  double PhiMass;


	//--------------------------------------------------------//
	//! Calculates daughter neutrino momentum from neutrino and phi masses. 
	/*! 
	\param parent is the index of the parent neutrino mass state. This is a column index in the rate matrix rate_mat.
	\param daughter is the index of the daughter neutrino mass state. This is a row index in the rate matrix rate_mat.
	\return retval is the desired momentum.
	*/

  double pstar(unsigned int daughter, unsigned int parent) const {
    if (NeutrinoMasses[parent] < NeutrinoMasses[daughter] + PhiMass) {
      throw std::runtime_error("non physical case");
    }
    double retval = (1.0 / (2.0 * NeutrinoMasses[parent])) *
                    sqrt((pow(NeutrinoMasses[parent], 2) - pow(NeutrinoMasses[daughter] + PhiMass, 2)) *
                         (pow(NeutrinoMasses[parent], 2) - pow(NeutrinoMasses[daughter] - PhiMass, 2)));
    return retval;
  }

	//--------------------------------------------------------//
	//! 	
	/*! Precomputes entries of pstar_mat
	Note: called within SetNeutrinoMasses(). Requires PhiMass to
	be previously set.
	*/

	void PreComputePstarMat(void){
    if (!NeutrinoMassesSet)
      throw std::runtime_error("neutrino masses not set");
    if (!PhiMassSet)
      throw std::runtime_error("phi mass not set");

	  for (unsigned int col=0; col<numneu; col++)
	  {
	    for (unsigned int row=0; row<col; row++)
	    {
				double val = pstar(row,col); 
	      gsl_matrix_set(pstar_mat,row,col,val);
	    }
		}
	}

	//--------------------------------------------------------//
	//! Given a double, finds the nearest double in the E_range array.	
	/*! 
	Given a double argument, searches E_range for the element of smallest
	absolute value difference from the argument. 
	\param value is the argument whose closest match in E_range we're looking for.
	\return the index of the closest element in E_range to value.
	*/

  unsigned int nearest_element(double value) const {
		//Generate a vector of absolute differences between value
		//and each element in E_range. 
    std::vector<double> diffs(E_range.size());
    for (size_t i = 0; i < E_range.size(); i++) {
      diffs[i] = fabs(value - E_range[i]);
    }

		//Return the index corresponding to the element of minimum difference
		//in E_range.
    return std::distance(diffs.begin(),
                         std::min_element(diffs.begin(), diffs.end()));
  }


	//--------------------------------------------------------//
	//! Constructs the Gamma decay matrix using rate_mat.	
	/*! 
		Constructing the initial Gamma decay matrix. The matrix is the sum of the 
		mass-basis projectors times the corresponding neutrino mass. The mass
		weighting factor is divided by the neutrino energy in GammaRho() yield a time
		dilation factor.
	*/

	void ConstructGammaMatrix(void){
		DT = squids::SU_vector(numneu);
		for(size_t i = 0; i < numneu; i++){
			double entry = NeutrinoMasses[i]*gsl_matrix_get(rate_mat,i,i);
			DT += entry*squids::SU_vector::Projector(numneu, i);
		}
	}


protected:

	//--------------------------------------------------------//
	//! FIXME: Carlos, can you elaborate on the implementation?	
	/*! 
	\param
	\return 
	*/

  void AddToPreDerive(double x) {
    if (!DecayParametersSet)
      throw std::runtime_error("decay parameters not set");
    if (!NeutrinoMassesSet)
      throw std::runtime_error("neutrino masses not set");
    if (!PhiMassSet)
      throw std::runtime_error("phi mass not set");

    for (int ei = 0; ei < ne; ei++) {
      // asumming same mass hamiltonian for neutrinos/antineutrinos
      squids::SU_vector h0 = H0(E_range[ei], 0);
      DT_evol[ei] = DT.Evolve(h0, (x - Get_t_initial()));
    }
  }

	//--------------------------------------------------------//
	//! FIXME: Carlos, can you elaborate on the implementation?	
	/*! 
	\param
	\return 
	*/

  void AddToWriteHDF5(hid_t hdf5_loc_id) const {
    // here we write the new parameters to be saved in the HDF5 file
    hsize_t dim[1]{1};

    // creating arrays to save stuff
    H5LTmake_dataset(hdf5_loc_id, "rate_matrix", 1, dim, H5T_NATIVE_DOUBLE,
                     0);
    H5LTmake_dataset(hdf5_loc_id, "phi_mass", 1, dim, H5T_NATIVE_DOUBLE,
                     0);
    H5LTmake_dataset(hdf5_loc_id, "lightest_nu_mass", 1, dim, H5T_NATIVE_DOUBLE,
                     0);

    // save decay rate matrix 
    for (unsigned int i = 0; i < numneu; i++) {
      for (unsigned int j = i; j < numneu; j++) {
        std::string rate_label =
            "rate" + std::to_string(i + 1) + std::to_string(j + 1);
        double rate_value = gsl_matrix_get(rate_mat,i,j);
        H5LTset_attribute_double(hdf5_loc_id, "rate_matrix", rate_label.c_str(),
                                 &rate_value, 1);
      }
    }

		// save phi mass	
		std::string phi_mass_label = "phimass";
		H5LTset_attribute_double(hdf5_loc_id, "phi_mass",
														 phi_mass_label.c_str(),
														 &PhiMass, 1);

		// save lightest neutrino mass	
		std::string light_mass_label = "lightmass";
		H5LTset_attribute_double(hdf5_loc_id, "lightest_nu_mass",
														 light_mass_label.c_str(),
														 &(NeutrinoMasses[0]), 1);
  }

	//--------------------------------------------------------//
	//! FIXME: Carlos, can you elaborate on the implementation?	
	/*! 
	\param
	\return 
	*/

  void AddToReadHDF5(hid_t hdf5_loc_id) {
    // read and set decay rates 
    for (unsigned int i = 0; i < numneu; i++) {
      for (unsigned int j = i; j < numneu; j++) {
        double rate_value;
        std::string rate_label =
            "rate" + std::to_string(i + 1) + std::to_string(j + 1);
        H5LTget_attribute_double(hdf5_loc_id, "rate_matrix", rate_label.c_str(),
                                 &rate_value);
				gsl_matrix_set(rate_mat,i,j,rate_value);
      }
    }

		//Construct the gamma matrix from rate_mat
		ConstructGammaMatrix();

    DecayParametersSet=true;

		// read and set phi mass
		std::string phi_mass_label = "phimass";
		H5LTget_attribute_double(hdf5_loc_id, "phi_mass",
														 phi_mass_label.c_str(),
														 &PhiMass);
		// read and set lightest neutrino mass
		double lightmass;
		std::string light_mass_label = "lightmass";
		H5LTget_attribute_double(hdf5_loc_id, "lightest_nu_mass",
														 light_mass_label.c_str(),
														 &lightmass);
		// calculate neutrino masses using the lightest mass and splittings.
		SetNeutrinoMasses(lightmass);
  }

	//--------------------------------------------------------//
	//! Returns the hamiltonian term corresponding to the Gamma matrix with additional interaction terms (can be switched on or off) described in the nuSQUIDS documentation under "GammaRho".
	/*! 
	If iincoherent_int switch is set to false, then this returns only the Gamma matrix
	corresponding to neutrino decay. If the switch is set to true, this function
	returns the sum of the decay Gamma matrix and additional neutrino-matter
	interactions described in the nuSQUIDS documentation.
	\param ie is the energy index of the desired matrix.
	\param irho is the neutrino/antineutrino index of the desired matrix, though the decay term is the same in each case.
	\return the decay/interaction matrix.
	*/

  squids::SU_vector GammaRho(unsigned int ie, unsigned int irho) const {
		//Toggling additional interaction terms.
		/*
		In each case, the 1/E divides the neutrino mass weighting applied in
		SetDecayMatrix() to DT. m/E=gamma, the time dilation factor, which 
		suppresses lab-frame decay rates.
		*/

    if (iincoherent_int)
      return nuSQUIDS::GammaRho(ie, irho) + DT_evol[ie] * (0.5 / E_range[ie]);
    else
      return DT_evol[ie] * (0.5 / E_range[ie]);
  }

	//--------------------------------------------------------//
	//! Returns the hamiltonian term corresponding to neutrino regeneration from interactions. This includes decay regeneration, as well as additional interaction regeneration terms which can be toggled on or off.
	/*! 
	If iincoherent_int switch is set to false, then this returns only the decay regeneration term. If set to true, this function returns the sum of additional interaction regeneration terms and the decay regeneration. The form of the decay regeneration term is given in the paper as the "R" matrix, while the additional interaction regeneration terms are described in the nuSQUIDS documentation under "InteractionsRho".
	\param ie is the energy index of the desired matrix.
	\param irho is the neutrino/antineutrino index of the desired matrix, though the decay term is the same in each case.
	\return the incoherent interaction matrix. 
	*/

  squids::SU_vector InteractionsRho(unsigned int ie, unsigned int irho) const {
    squids::SU_vector decay_regeneration(numneu);
		//ie is the index of the daughter neutrino energy. Fetch the energy.
    double Ef = E_range[ie];
		//Scanning through daughter neutrino mass states
    // i-daughter index
    for (size_t i = 0; i < numneu; i++) {
			//Scanning through parent neutrino mass states
      // j-parent index
      for (size_t j = i + 1; j < numneu; j++) {
				//Calculate parent-rest -- lab frame boost factor from daughter energy.
        double gamma = Ef/gsl_matrix_get(pstar_mat,i,j);
				//Calculate parent neutrino energy 
        double E0 = NeutrinoMasses[j]*gamma;
				//Finds the E_range energy closest to the calculated parent neutrino energy.
        size_t E0_index = nearest_element(E0);
				/*
				Determines expectation of parent mass state in E0_index energy bin given
				the density matrix rho at the current time step. Multiplies this by the
				parent-frame decay rate in the corresponding channel, modified by the time
				dialation factor from the parent-lab frame gamma factor. Finally, this scalar
				weight is multiplied by the daughter mass state projector in the daughter
				energy bin. This is the "R" matrix in the hamiltonian specified in the paper.
				*/
        decay_regeneration += (state[E0_index].rho[irho]
					*evol_b0_proj[irho][j][E0_index])
					*(gsl_matrix_get(rate_mat,i,j)/gamma)
					*(evol_b0_proj[irho][i][ie]);
      }
    }

		//Toggling additional regeneration terms.
    if (iincoherent_int)
      return nuSQUIDS::InteractionsRho(ie, irho) + decay_regeneration;
    else
      return decay_regeneration;
  }

public:

	//--------------------------------------------------------//
	//! nuSQUIDSDecay constructor.  	
	/*! 
	Calls nuSQUIDS parent constructor and allocates memory for member
	matrices.
	\param e_nodes is the array of neutrino propagation energies.
	\param numneu_ number of neutrino states in the system. Defaults to 3.
	\param NT_ is the neutrino type from (neutrino/antineutrino). Defaults to both.	
	\param iinteraction_ FIXME: Carlos: does this just get fed to the parent constructor? 
	*/

  nuSQUIDSDecay(marray<double, 1> e_nodes, unsigned int numneu_ = 3,
                NeutrinoType NT_ = NeutrinoType::both,
                bool iinteraction_ = true)
      : nuSQUIDS(e_nodes, numneu_, NT_, iinteraction_,
                 std::make_shared<
                     nusquids::NeutrinoDISCrossSectionsFromTablesExtended>()) {

    // just allocate some matrices
    DT_evol.resize(ne);
    for (int ei = 0; ei < ne; ei++) {
      DT_evol[ei] = squids::SU_vector(nsun);
    }

    // allocating space for neutrino masses
    NeutrinoMasses.resize(numneu);

    // allocating memory for rate matrix
    rate_mat = gsl_matrix_alloc(numneu,numneu);

    // allocating memory for pstar matrix
    pstar_mat = gsl_matrix_alloc(numneu,numneu);
  }

	//--------------------------------------------------------//
	//! nuSQUIDSDecay destructor. 	
	/*! 
	Freeing memory allocated to gsl_matrices.
	*/

  ~nuSQUIDSDecay(){
    gsl_matrix_free(rate_mat);
    gsl_matrix_free(pstar_mat);
  }

	//--------------------------------------------------------//
	//! Fills the decay rate matrix rate_mat given a matrix tau of lifetimes. 	
	/*! 
	\param tau is a gsl_matrix containing parent neutrino lifetimes in eV^-1. This matrix is strictly upper triangular. The column indicates the parent neutrino, and the row the daughter.
	*/

  void SetDecayMatrix(gsl_matrix* tau){
		//Perform some matrix dimension checks.
    if (rate_mat->size1 != tau->size1){
      throw std::runtime_error("size1 mismatch while constructing decay matrix.");
    }
    if (rate_mat->size2 != tau->size2){
      throw std::runtime_error("size2 mismatch while constructing decay matrix.");
    }

		//Zero out the elements of the rate matrix.
	  gsl_matrix_set_zero(rate_mat);

		/*
		Setting each strictly upper-triangular decay rate to 1/lt, where lt is
		the corresponding lifetime. The diagonal elements are used to construct
		the decay Gamma matrix, so they contain the sum of the rates in each column.
		That is to say, each diagonal initially corresponds to the total rate of 
		depletion of the corresponding parent mass state.
		*/
	  for (size_t col=0; col<numneu; col++)
	  {
	    double colrate=0;
	    for (size_t row=0; row<col; row++)
	    {
	      double rate = 1.0/gsl_matrix_get(tau,row,col);
	      gsl_matrix_set(rate_mat,row,col,rate);
	      colrate+=rate;
	    }
	    gsl_matrix_set(rate_mat,col,col,colrate);
		}

		//Construct the gamma matrix from rate_mat
		ConstructGammaMatrix();

    DecayParametersSet=true;
  }

	//--------------------------------------------------------//
	//! Returns the matrix of lifetimes, rate_mat. 	
	/*! 
	\return rate_mat: the lifetime matrix. 
	*/

	gsl_matrix* GetDecayMatrix(void){
		return rate_mat;
	}

	//--------------------------------------------------------//
	//! Toggles additional neutrino interactions.	
	/*! 
			The switch, iincoherent_int, toggles additional interactions
			in both GammaRho() and InteractionsRho(). These interactions are
			described in the nuSQUIDS documentation under the corresponding
			function names. If the switch is set to true, these terms are 
			added to the decay Gamma matrix and decay "R" matrix. If set 
			to false, then the only interaction in play is neutrino decay.
	\param opt: the boolean value to toggle interactions.
	*/

  void SetIncoherentInteractions(bool opt) { iincoherent_int = opt; }

	//--------------------------------------------------------//
	//! Toggles decay regeneration. 	
	/*!
			The switch is internal to SQUIDS/nuSQUIDS. If set to true, the 
			terms returned by InteractionsRho() are present in the evolution
			equation. If not, there is no regeneration. It is perhaps helpful
			to include a truth table below.
			iincoherent_int | DecayRegeneration | Physics
			--------------- | ----------------- | -------------------------------
			True            | True              | All terms included. 
			True            | False             | All except regeneration terms.
			False           | True              | Gamma + R (decay with regen only).
			False           | False             | Gamma only (decay only without regen).
	\param opt: the boolean value to toggle regeneration.
	*/

  void SetDecayRegeneration(bool opt) { Set_OtherRhoTerms(opt); }

	//--------------------------------------------------------//
	//! Sets neutrino masses given splittings and the lightest mass. 	
	/*! 
			One must call SetPhiMass before calling this function! PreComputePstarMat()
			is called before the function terminates, and it requires a phi mass.
			The mass splittings must be set before calling this function.
			Default splittings will be used if they are not set.
	\param lightmass is the lightest neutrino mass.
	*/

	void SetNeutrinoMasses(double lightmass){
    if (!PhiMassSet){
      throw std::runtime_error("Phi mass not set. Set phi mass before setting neutrino masses.");
    }
		//Set the first mass to lightmass.
		NeutrinoMasses[0] = lightmass;
		/*
		Calculate the ith mass using the (i-1)th mass and the corresponding
		mass-squared splitting.
		*/
		for (unsigned int i=1; i<numneu; i++){
			NeutrinoMasses[i] = sqrt(Get_SquareMassDifference(i) + NeutrinoMasses[0]*NeutrinoMasses[0]);
		}
		NeutrinoMassesSet=true;
		PreComputePstarMat();
	}
		
	//--------------------------------------------------------//
	//! Gets the neutrino mass corresponding to the given state index. 	
	/*! 
	\param index is the index of the desired mass state.
	\return NeutrinoMasses[index] is the corresponding neutrino mass.
	*/

	double GetNeutrinoMass(unsigned int index){
		//Check that the requested index is not out of range.
    if (index>=numneu){
      throw std::runtime_error("Index larger than numneu-1");
    }
		return NeutrinoMasses[index];
	}

	//--------------------------------------------------------//
	//! Sets the mass of the phi particle.	
	/*! 
			This function must be called before SetNeutrinoMasses().
	\param PhiMass_ is the provided mass of the phi particle.
	*/

	void SetPhiMass(double PhiMass_){
		PhiMass = PhiMass_;
		PhiMassSet=true;
	}

	//--------------------------------------------------------//
	//! Gets the mass of the phi particle.
	/*! 
	\return PhiMass is the mass of the phi particle. 
	*/

	double GetPhiMass(void){
		return PhiMass;
	}

}; // close class
} // close nusquids namespace

#endif // nusquidslv_h
