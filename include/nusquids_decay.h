#ifndef nusquids_decay_H
#define nusquids_decay_H

/*
Header implementing NuSQUIDSDecay class. See documentation for technical details.
See arXiv:1711.05921 (referred to here as [1]) for details about the decay model.

For more information please email:

Alexander Moss (zander@mit.edu)
Marjon Moulai (marjon@mit.edu)
Carlos Arguelles (caad@mit.edu)
Janet Conrad (conrad@mit.edu)
*/

//Get mathematical constants.
#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>
#include <iostream>
#include <nuSQuIDS/nuSQuIDS.h>
#include "exCross.h"

namespace nusquids {

class nuSQUIDSDecay : public nuSQUIDS {
private:
	//-----------------------------Variables----------------------------//
	//! Toggles Majorana or Dirac neutrinos
	/*!
	If neutrinos are Dirac, the chirality violating process
	induces decay into right-handed neutrinos, which do not
	contribute to regeneration.
	If neutrinos are Majorana, the chirality violating process
	induces decay into left-handed antineutrinos, which do contribute
	to regeneration. This switch alters the regeneration calculation
	(see InteractionsRho())
	Default: false.
	*/
	bool majorana=false;

	//! Toggles additional incoherent interactions
	/*!
	See Set_DecayRegeneration() for details.
	*/
	bool ihard_interactions=false;

	//Chirality Preserving Process or Chirality Violating Process
	enum{CPP,CVP};

	//! Toggles scalar or pseudoscalar couplings.
	/*!
	For the time being, we are considering either purely scalar or
	purely pseudoscalar couplings. This switch takes the value false in the 
	scalar case, and true in the pseudoscalar case. Its value will be set in 
	the constructor.
	*/
	bool pscalar;

	//Parent<->Row, Daughter<->Column (lower triangular)
	//One for each of {CPP,CVP}

	//! The coupling matrices, g_ij.
	/*!
	The row index corresponds to the parent mass state, and the column to
	the daughter. The matrix will then be strictly lower triangular.
	*/
	gsl_matrix* couplings;

	//! The decay rate matrices, Gamma_ij.
	/*!
	There are two rate matrices, one for each element of
	{CPP,CVP}. The row index of each
	corresponds to the parent mass state, and the column to
	the daughter. The matrix is then strictly lower-triangular.
	These decay rates are computed in the *rest frame* of the
	parent. Eqns. (2) and (3) in [1] are lab-frame, and differ
	by a factor of 1/gamma.
	*/
	gsl_matrix* rate_matrices[2];

	//! Vector of neutrino masses.
	/*!
	The lightest mass may be zero, but all other neutrino masses
	must be non-zero.
	*/
	std::vector<double> m_nu;

	//! The "Gamma" matrix appearing in the full Hamiltonian.
	/*!
	Encodes the loss of probability current from a given energy bin,
	due to decay.
	*/
	squids::SU_vector DT;

	//! The "Gamma" matrix in evolving basis.
	std::vector<squids::SU_vector> DT_evol;

	//----------------------------------Functions---------------------------------//
	//Trying to keep the "model-specific" functions in private,
	//and have moved all functions related to more general decay models
	//to protected.

	//! Kinematic funtion m_j*f(x_i/x_j), related to equation (4a) in [1].
	/*!
	In this and the other auxiliary kinematic function definitions, the modification
	is to make the m_j->0 limit numerically well defined. Luckily, the cmath pow()
	function implicitly takes the limit: pow(0,0)=1, so we're happy.
	\param m_i is the parent neutrino mass.
	\param m_j is the daughter neutrino mass.
	*/
	double f(double m_i, double m_j) const {
		double result = m_i/2.0 + 2.0*m_j + (2.0*m_j/m_i)*(m_j*log(m_i) - log(pow(m_j,m_j))) - (2.0*m_j*m_j*m_j/(m_i*m_i)) - (m_j*m_j*m_j*m_j/(2.0*m_i*m_i*m_i));
		return result;
	}
	//! Kinematic funtion m_j*g(x_i/x_j), related to equation (4b) in [1].
	/*!
	\param m_i is the parent neutrino mass.
	\param m_j is the daughter neutrino mass.
	*/
	double g(double m_i, double m_j) const {
		double result = m_i/2.0 - 2.0*m_j + (2.0*m_j/m_i)*(m_j*log(m_i) - log(pow(m_j,m_j))) + (2.0*m_j*m_j*m_j/(m_i*m_i)) - (m_j*m_j*m_j*m_j/(2.0*m_i*m_i*m_i));
		return result;
	}
	//! Kinematic funtion m_j*k(x_i/x_j), related to equation (4c) in [1].
	/*!
	\param m_i is the parent neutrino mass.
	\param m_j is the daughter neutrino mass.
	*/
	double k(double m_i, double m_j) const {
		double result = m_i/2.0 - (2.0*m_j/m_i)*(m_j*log(m_i) - log(pow(m_j,m_j))) - (m_j*m_j*m_j*m_j/(2.0*m_i*m_i*m_i));
		return result;
	}

	//! Computes the four decay rate matrices using the coupling matrix couplings.
	/*!
	This function implements equations (2) and (3) of [1] to generate the two partial
	rate matrices corresponding to each decay channel in {CPP,CVP}.
	Decay rates are in the rest frame of the parent neutrino.
	*/
	void Compute_Rate_Matrices(){
		//Compute *rest frame* decay rate matrices.
		if (!pscalar){
			//CPP,SCALAR
			rate_matrices[CPP] = gsl_matrix_alloc(numneu,numneu);
			for (unsigned int i=0; i<numneu; i++){
				for (unsigned int j=0; j<i; j++){
					double g_ij = gsl_matrix_get(couplings,i,j);
					double rate = (1.0/(16.0*M_PI))*g_ij*g_ij*f(m_nu[i],m_nu[j]);
					gsl_matrix_set(rate_matrices[CPP],i,j,rate);
				}
			}
			//CVP,SCALAR
			rate_matrices[CVP] = gsl_matrix_alloc(numneu,numneu);
			for (unsigned int i=0; i<numneu; i++){
				for (unsigned int j=0; j<i; j++){
					double g_ij = gsl_matrix_get(couplings,i,j);
					double rate = (1.0/(16.0*M_PI))*g_ij*g_ij*k(m_nu[i],m_nu[j]);
					gsl_matrix_set(rate_matrices[CVP],i,j,rate);
				}
			}
		}
		if (pscalar){
			//CPP,PSEUDOSCALAR
			rate_matrices[CPP] = gsl_matrix_alloc(numneu,numneu);
			for (unsigned int i=0; i<numneu; i++){
				for (unsigned int j=0; j<i; j++){
					double g_ij = gsl_matrix_get(couplings,i,j);
					double rate = (1.0/(16.0*M_PI))*g_ij*g_ij*g(m_nu[i],m_nu[j]);
					gsl_matrix_set(rate_matrices[CPP],i,j,rate);
				}
			}
			//CVP,PSEUDOSCALAR
			rate_matrices[CVP] = gsl_matrix_alloc(numneu,numneu);
			for (unsigned int i=0; i<numneu; i++){
				for (unsigned int j=0; j<i; j++){
					double g_ij = gsl_matrix_get(couplings,i,j);
					double rate = (1.0/(16.0*M_PI))*g_ij*g_ij*k(m_nu[i],m_nu[j]);
					gsl_matrix_set(rate_matrices[CVP],i,j,rate);
				}
			}
		}
	}

	//! Returns a sum of all Hamiltonian interaction terms except DT, the "Gamma" matrix.
	/*!
	The contribution from decay regeneration (the "R" terms from eqns. (18) and (19) in [1]) is computed,
	and then additional incoherent interactions are added internally by nuSQuIDS if iinteractions is true.
	These additional interaction regeneration terms are described in the nuSQUIDS documentation under "InteractionsRho".
	If majorana is true, there are regeneration contributions from both CVP and CPP. Otherwise, there is no contribution 
	because right-handed neutrinos are sterile.
	Note that the contribution from CVP in the majorana case converts neutrinos to antineutrinos.
	\param iedaughter is the energy index of the daughter density matrix.
	\param irho is the neutrino/antineutrino index of the desired matrix.
	\return the incoherent interaction matrix.
	*/

	squids::SU_vector InteractionsRho(unsigned int iedaughter, unsigned int irho) const {
		squids::SU_vector decay_regeneration(numneu);
		// Get the daughter neutrino energy.
		double edaughter = E_range[iedaughter];
		// j-daughter index
		for (size_t j = 0; j < numneu; j++) {
			// i-parent index: sum over contributions from all
			// states heavier than m_i
			for (size_t i = j+1; i < numneu; i++) {
				/*
				Decay kinematics dictate an integral of the regeneration
				contribution over parent momenta in the range [edaughter,edaughter*x_ij^2]
				See (18) and (19) in [1]. Here, we approximate the integral with a
				left-rectangular sum over energy bins in this range.
				We change from x_ij to y_ij=1/x_ij to clean up behavior in the m_j-> case.
				This is just an algebraic manipulation to keep everything stable.
				*/

				//parent-to-daughter mass ratio
				double xij = m_nu[i]/m_nu[j];
				double ieparent_high = nearest_element(edaughter*(xij*xij));
				//Check that m_nu[j] is not too close to zero.
				//if it isn't, we can use the formulae directly from the paper.	
				if (fabs(m_nu[j]-0.0)>1e-6){
					// i-energy (parent energy) index
					//left-rectangular integral approximation
					for (size_t ieparent = iedaughter; ieparent < ieparent_high-1; ieparent++) {
						//get parent neutrino energy
						double eparent = E_range[ieparent];
						//boost factor to lab frame
						double gamma = eparent/m_nu[i];
						double delta_eparent = E_range[ieparent+1]-E_range[ieparent];
						if (!pscalar){
							decay_regeneration += (delta_eparent)*(state[ieparent].rho[irho]
												*evol_b0_proj[irho][i][ieparent])*
												(xij*xij/(xij*xij-1))*
												(1/(eparent*eparent*edaughter))*
												((gsl_matrix_get(rate_matrices[CPP],i,j)/gamma)*
												pow(eparent+xij*edaughter,2)/pow(xij+1,2))*
												(evol_b0_proj[irho][j][iedaughter]);
						}
						if (pscalar){
							decay_regeneration += (delta_eparent)*(state[ieparent].rho[irho]
												*evol_b0_proj[irho][i][ieparent])*
												(xij*xij/(xij*xij-1))*
												(1/(eparent*eparent*edaughter))*
												((gsl_matrix_get(rate_matrices[CPP],i,j)/gamma)*
												pow(eparent-xij*edaughter,2)/pow(xij-1,2))*
												(evol_b0_proj[irho][j][iedaughter]);
						}
					}
	
					//Include chirality-violating term. Majorana CVP sends nu to nubar.
					//The procedure is the same, but the parent irho index is inverted.
					unsigned int parent_irho;
					if (irho==0) parent_irho=1;
					if (irho==1) parent_irho=0;
					// i-energy (parent energy) index
					//left-rectangular integral approximation
					for (size_t ieparent = iedaughter; ieparent < ieparent_high-1; ieparent++) {
						double eparent = E_range[ieparent];
						double gamma = eparent/m_nu[i];
						double delta_eparent = E_range[ieparent+1]-E_range[ieparent];
						if (!pscalar){
							decay_regeneration += (delta_eparent)*(state[ieparent].rho[parent_irho]
												*evol_b0_proj[parent_irho][i][ieparent])*
												(xij*xij/(xij*xij-1))*
												((eparent-edaughter)/(eparent*eparent*edaughter))*
												((gsl_matrix_get(rate_matrices[CVP],i,j)/gamma)*
												(edaughter*pow(xij,2)-eparent)/pow(xij+1,2))*
												(evol_b0_proj[irho][j][iedaughter]);
						}
						if (pscalar){
							decay_regeneration += (delta_eparent)*(state[ieparent].rho[parent_irho]
												*evol_b0_proj[parent_irho][i][ieparent])*
												(xij*xij/(xij*xij-1))*
												((eparent-edaughter)/(eparent*eparent*edaughter))*
												((gsl_matrix_get(rate_matrices[CVP],i,j)/gamma)*
												(edaughter*pow(xij,2)-eparent)/pow(xij-1,2))*
												(evol_b0_proj[irho][j][iedaughter]);
						}
					}
				}
				
				//If m_nu[j] is too close to zero, xij diverges, and we switch to an alternative
				//form for the differential decay rates, in terms of yij=1/xij.
				else{
					double yij = m_nu[j]/m_nu[i];
	
					for (size_t ieparent = iedaughter; ieparent < ieparent_high-1; ieparent++) {
						//get parent neutrino energy
						double eparent = E_range[ieparent];
						//boost factor to lab frame
						double gamma = eparent/m_nu[i];
						double delta_eparent = E_range[ieparent+1]-E_range[ieparent];
						if (!pscalar){
							decay_regeneration += (delta_eparent)*(state[ieparent].rho[irho]
												*evol_b0_proj[irho][i][ieparent])*
												(1/(1-yij*yij))*
												(1/(eparent*eparent*edaughter))*
												((gsl_matrix_get(rate_matrices[CPP],i,j)/gamma)*
												pow(eparent*yij+edaughter,2)/pow(yij+1,2))*
												(evol_b0_proj[irho][j][iedaughter]);
						}
						if (pscalar){
							decay_regeneration += (delta_eparent)*(state[ieparent].rho[irho]
												*evol_b0_proj[irho][i][ieparent])*
												(1/(1-yij*yij))*
												(1/(eparent*eparent*edaughter))*
												((gsl_matrix_get(rate_matrices[CPP],i,j)/gamma)*
												pow(eparent*yij-edaughter,2)/pow(1-yij,2))*
												(evol_b0_proj[irho][j][iedaughter]);
						}
					}
	
	
					//Include chirality-violating term. Majorana CVP sends nu to nubar.
					//The procedure is the same, but the parent irho index is inverted.
					unsigned int parent_irho;
					if (irho==0) parent_irho=1;
					if (irho==1) parent_irho=0;
					// i-energy (parent energy) index
					//left-rectangular integral approximation
					for (size_t ieparent = iedaughter; ieparent < ieparent_high-1; ieparent++) {
						double eparent = E_range[ieparent];
						double gamma = eparent/m_nu[i];
						double delta_eparent = E_range[ieparent+1]-E_range[ieparent];
						if (!pscalar){
							decay_regeneration += (delta_eparent)*(state[ieparent].rho[parent_irho]
												*evol_b0_proj[parent_irho][i][ieparent])*
												(1/(1-yij*yij))*
												((eparent-edaughter)/(eparent*eparent*edaughter))*
												((gsl_matrix_get(rate_matrices[CVP],i,j)/gamma)*
												(edaughter-pow(yij,2)*eparent)/pow(yij+1,2))*
												(evol_b0_proj[irho][j][iedaughter]);
						}
						if (pscalar){
							decay_regeneration += (delta_eparent)*(state[ieparent].rho[parent_irho]
												*evol_b0_proj[parent_irho][i][ieparent])*
												(1/(1-yij*yij))*
												((eparent-edaughter)/(eparent*eparent*edaughter))*
												((gsl_matrix_get(rate_matrices[CVP],i,j)/gamma)*
												(edaughter-pow(yij,2)*eparent)/pow(1-yij,2))*
												(evol_b0_proj[irho][j][iedaughter]);
						}
					} 

				} //Close second m_nu[j] conditional.
			} //Close i loop
		} //Close j loop 

		//Toggling additional regeneration terms (from nuSQuIDS).
		if (ihard_interactions){
			return nuSQUIDS::InteractionsRho(iedaughter, irho) + decay_regeneration;
		}
		else{
			return decay_regeneration;
		}
	}

protected:
	//! Given a double, finds the nearest double in the E_range array.
	/*!
	Given a double argument, searches E_range for the element of smallest
	absolute value difference from the argument.
	\param value is the argument whose closest match in E_range we're looking for.
	\return the index of the closest element in E_range to value.
	*/
	unsigned int nearest_element(double value) const {
		//If value is larger than last element, return
		//last element.
		if (value>E_range[E_range.size()-1]){
			return E_range.size()-1;
		}
		//Otherwise, search in range.
		std::vector<double> diffs(E_range.size());
		for (size_t i = 0; i < E_range.size(); i++) {
			diffs[i] = fabs(value - E_range[i]);
		}
		return std::distance(diffs.begin(),std::min_element(diffs.begin(), diffs.end()));
	}

	//! Prints the contents of a squids::SU_vector. (Useful for debugging)
	/*!
	\param mat is the SU_vector.
	\param dim is its dimension (assume square).
	\param mname is the name of the matrix to print with it's entries.
	*/
	void printmat(const squids::SU_vector mat, unsigned int dim,std::string mname) const {
		std::cout << "Matrix: " << mname << std::endl;
		unsigned int i, j;
		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				std::cout << (mat)[j + dim * i] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	//! Prints the contents of a gsl_matrix. (Useful for debugging)
	/*!
	\param mat is the gsl_matrix.
	\param dim is its dimension (assume square).
	\param mname is the name of the matrix to print with it's entries.
	*/
	void gslprintmat(gsl_matrix* mat, unsigned int dim,std::string mname) const {
		std::cout << "Matrix: " << mname << std::endl;
		unsigned int i, j;
		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				std::cout << gsl_matrix_get(mat,i,j) << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	//! Checks a pair of gsl_matrices for matching row and column dimensions.
	/*!
	\param m1 the first matrix.
	\param m2 the second matrix.
	*/
	void Check_Matrix_Size(gsl_matrix* m1, gsl_matrix* m2) const {
		if (m1->size1 != m2->size1){
			throw std::runtime_error("size1 mismatch while copying matrix.");
		}
		if (m1->size2 != m2->size2){
			throw std::runtime_error("size2 mismatch while copying matrix.");
		}
	}

	//! Copies a rank-two array of gsl_matrices to the rate_matrices data member.
	/*!
	First allocates empty gsl_matrices to the rate_matrices[] pointer-array,
	then compares source and destination matrix dimensions, then copies.
	This is only used internally to set the rate_matrices data member. Except
	in the case of the special "partial rates" constructor, the user should provide
	a *coupling* matrix, and the rate matrices will be computed internally.
	\param rate_matrices_ the input array.
	*/
	void Set_Rate_Matrices(gsl_matrix* rate_matrices_[2]){
		for (unsigned int chi=0; chi<2; chi++){
			rate_matrices[chi] = gsl_matrix_alloc(numneu,numneu);
			Check_Matrix_Size(rate_matrices[chi],rate_matrices_[chi]);
			gsl_matrix_memcpy(rate_matrices[chi],rate_matrices_[chi]);
		}
	}

	//! Set the Lagrangian couplings between mass states.
	/*!
	First allocate memory gsl_matrix to the couplings pointer.
	Then, check source and target matrix dimensions match, and copy the matrix.
	See #couplings.
	The user specifies whether the couplings are scalar or pseudoscalar 
	in the constructor.
	\param couplings_ is a gsl_matrix pointer. 	
	*/
	void Set_Couplings(gsl_matrix* couplings_){
		couplings = gsl_matrix_alloc(numneu,numneu);
		Check_Matrix_Size(couplings,couplings_);
		gsl_matrix_memcpy(couplings,couplings_);
	}

	//! Computes DT, the "Gamma" term in the Hamiltonian, using rate_matrices.
	/*!
	In the mass basis, DT is a diagonal matrix such that the (i,i) entry
	is the sum over all rate matrices of their ith rows, weighted by m_i,
	the corresponding neutrino mass. Essentially, this matrix encodes the
	total rate of decay of mass state i to all lighter states, including
	all decay channels ({CPP,CVP}) if the neutrino
	is Majorana, and only the channel (CVP) if it 
	is Dirac. The rate is weighted by the mass m_i for convenience, so that, 
	in GammaRho(), one can simply divide DT by the energy of the parent 
	neutrino, and each rate will acquire the proper factor of 1/gamma 
	characterizing lab-frame decay retarded by time-dialation.
	*/
	void Compute_DT(){
		DT = squids::SU_vector(numneu);
		//Include chirality preserving processes only in majorana case.
		int chi_min;
		if(majorana){chi_min=0;}
		else{chi_min=1;}
		//Sum over parent mass states.
		for(size_t i = 0; i < numneu; i++){
			double rate=0;
			//Sum over daughter mass states.
			for(size_t j=0; j<i; j++){
				//Sum over all decay channels.
				for (size_t chi=chi_min; chi<2; chi++){
					rate+=gsl_matrix_get(rate_matrices[chi],i,j);	
				}
			}
			//Weight rate by m_i, and add a projector to the m_i state,
			//multiplied by the weighted rate, to DT.
			DT += m_nu[i]*rate*squids::SU_vector::Projector(numneu, i);
		}
	}

	//! Evolves the interaction picture DT Hamiltonian term.
	/*!
	\param x the target evolution time.
	*/
	void AddToPreDerive(double x) {
		for (int ei = 0; ei < ne; ei++) {
			// asumming same mass hamiltonian for neutrinos/antineutrinos
			squids::SU_vector h0 = H0(E_range[ei], 0);
			DT_evol[ei] = DT.Evolve(h0, (x - Get_t_initial()));
		}
	}

    //! Returns the hamiltonian term corresponding to the "Gamma" matrix with additional terms from nuSQuIDS.
    /*! 
    If interaction switch is set to false, then this returns only the Gamma matrix
    corresponding to neutrino decay (DT_evol properly weighted). If the switch is set 
	to true, this function returns the sum of the decay Gamma matrix and additional neutrino-matter
    interactions described in the nuSQUIDS documentation under "GammaRho".
    \param ie is the energy index of the desired term.
    \param irho is the neutrino/antineutrino index of the desired term, though the decay term is the same in each case.
    \return the modified "Gamma" matrix.
    */
	squids::SU_vector GammaRho(unsigned int ie, unsigned int irho) const {
		if (ihard_interactions)
			return nuSQUIDS::GammaRho(ie, irho) + DT_evol[ie] * (0.5 / E_range[ie]);
		else
			return DT_evol[ie] * (0.5 / E_range[ie]);
	}


public:
	//! Basic nuSQUIDSDecay constructor.
	/*!
	Calls nuSQUIDS parent constructor and allocates memory for member
	matrices and neutrino masses. This constuctor should not be called
	directly. It is just an encapsulation of basic funcitonality to be
	called by the two, more complete, overloaded constructors.
	\param e_nodes is the array of neutrino propagation energies.
	\param numneu_ is the number of neutrino states in the system. Defaults to 3.
	\param NT_ is the neutrino type from (neutrino/antineutrino/both). Defaults to both.	
	\param iinteraction_ is a switch for incoherent interactions. See Set_DecayRegeneration() .
	*/
	nuSQUIDSDecay(marray<double, 1> e_nodes, unsigned int numneu_ = 3,
					NeutrinoType NT_ = NeutrinoType::both,
					bool ihard_interactions = true)
					: nuSQUIDS(e_nodes, numneu_, NT_, ihard_interactions), ihard_interactions(ihard_interactions) {
    // activate gamma term always
    Set_NonCoherentRhoTerms(true);
    // do not allow the system to skip the numerics
    Set_AllowConstantDensityOscillationOnlyEvolution(false);
		// just allocate some matrices
		DT_evol.resize(ne);
		for (int ei = 0; ei < ne; ei++) {
			DT_evol[ei] = squids::SU_vector(nsun);
		}
		// allocating space for neutrino masses
		m_nu.resize(numneu);
	}

	//! nuSQUIDSDecay "majorana coupling" constructor.
	/*!
	Calls the basic constructor, and then sets neutrino masses, phi mass,
	the coupling matrices, as well as switches for incoherent interactions,
	decay regeneration (See SetIncoherentInteractions(),
	Set_DecayRegeneration()). The constructor then calls
	Compute_Rate_Matrices() to calculate decay rates as functions of masses and couplings,
	and then calls Compute_DT() to calculate the "Gamma" matrix decay term as a function
	of masses and decay rates.
	This constructor only applies to tha majorana case, where it is preferred because it 
	allows the user to input Lagrangian parameters only, as opposed to manually computing 
	partial decay rate matrices and passing them to nuSQuIDS decay. This both simplifies 
	the interface, and guarantees that only physical decay rate matrices are used in the 
	simulations (an arbitrary selection of rate matrices might not be possible to generate 
	with a Lagrangian of the form (1) in [1]). There is no majorana boolean parameter in this
	case. If the user wishes to implement a dirac neutrino, they should use the "partial rate"
	constructor, described below. Note: only m_1 may be set to zero!
	This calculation *does not apply* if more than one neutrino mass is zero.

	\param e_nodes is the array of neutrino propagation energies.
	\param numneu_ is the number of neutrino states in the system. Defaults to 3.
	\param NT_ is the neutrino type from (neutrino/antineutrino/both). Defaults to both.	
	\param iinteraction_ is a switch for incoherent interactions. See Set_DecayRegeneration() .
	\param pscalar_ is a switch for scalar/pseudoscalar couplings. See #pscalar.
	\param decay_regen_ is a switch for decay regeneration. See Set_DecayRegeneration()
	\param m_nu_ is a vector of neutrino masses. See #m_nu .
	\param couplings_ is a gsl_matrix* pointer. See #couplings .
	*/
	nuSQUIDSDecay(marray<double, 1> e_nodes, unsigned int numneu_,
					NeutrinoType NT_, bool ihard_interactions_,
					bool decay_regen_, bool pscalar_, 
					std::vector<double> m_nu_, gsl_matrix* couplings_ 
					):
					nuSQUIDSDecay(e_nodes,numneu_,NT_,ihard_interactions_){
		ihard_interactions=ihard_interactions_;
		pscalar=pscalar_;
		majorana=true;
		Set_DecayRegeneration(decay_regen_);
		m_nu=m_nu_;
		Set_Couplings(couplings_);
		Compute_Rate_Matrices();
		Compute_DT();
	}

	//! nuSQUIDSDecay "partial rate" constructor.
	/*!
	Calls the basic constructor, and then sets neutrino masses,
	the four partial rate matrices, as well as switches for incoherent interactions,
	decay regeneration, and majorana/dirac neutrinos (See SetIncoherentInteractions(),
	Set_DecayRegeneration(), and). The constructor then allocates memory
	to the unused #couplings and calls Compute_DT() to calculate the "Gamma" matrix
	decay term as a function of masses and decay rates.
	This constructor is used in the analysis in [1] because it allows the user to input
	partial decay rates directly, which is useful for characterizing the effect of
	neutrino lifetimes on the evolution of the neutrino system. However, one must be careful
	with the interpretation of these "partial lifetimes": these four matrices are functions of
	the scalar and pseudoscalar couplings in (1) from [1], so an arbitrary choice for rate_matrices_
	might not be physical.
	If the neutrino is majorana, then the "coupling constructor" is preferred, but this can be used
	if the user remembers to set the majorana bool to true. If the user wishes to implement a Dirac
	neutrino, they should use this constructor, set majorana to false, and compute the decay rate
	matrices appropriate to the corresponding lagrangian ahead of time, then feed them to this constructor. Note: only m_1 may be set to zero!
	This calculation *does not apply* if more than one neutrino mass is zero.


	\param e_nodes is the array of neutrino propagation energies.
	\param numneu_ is the number of neutrino states in the system. Defaults to 3.
	\param NT_ is the neutrino type from (neutrino/antineutrino/both). Defaults to both.	
	\param iinteraction_ is a switch for incoherent interactions. See Set_DecayRegeneration() .
	\param pscalar_ is a switch for scalar/pseudoscalar couplings. See #pscalar.
	\param decay_regen_ is a switch for decay regeneration. See Set_DecayRegeneration()
	\param majorana_ is a switch for Majorana/Dirac neutrinos. See #majorana .
	\param m_nu_ is a vector of neutrino masses. See #m_nu .
	\param rate_matrices_ is a two-element array of gsl_matrix* pointers. See #rate_matrices .
	*/
	nuSQUIDSDecay(marray<double, 1> e_nodes, unsigned int numneu_,
					NeutrinoType NT_, bool ihard_interactions,
					bool decay_regen_, bool pscalar_, 
					bool majorana_, std::vector<double> m_nu_,
					gsl_matrix* rate_matrices_[2]
					):
					nuSQUIDSDecay(e_nodes,numneu_,NT_,ihard_interactions){
		ihard_interactions=ihard_interactions;
		pscalar=pscalar_;
		majorana=majorana_;
		/*
		If the neutrino is majorana, both CPP and CVP processes contribute to regeneration,
		and CVP sends neutrinos to antineutrinos. If the neutrino is Dirac, there is no CPP,
		and the CVP sends left-handed neutrinos to right-handed neutrinos, which are sterile,
		and therefore do not contribute to regeneration of visible neutrino flux. That is to say,
		the regeneration terms are automatically zero if the neutrino is Dirac.
		*/
		if (majorana){
			Set_DecayRegeneration(decay_regen_);
		}
		else{
			Set_DecayRegeneration(false);
		}
		m_nu=m_nu_;
		Set_Rate_Matrices(rate_matrices_);
		couplings = gsl_matrix_alloc(numneu,numneu);
		Compute_DT();
	}

	//! nuSQUIDSDecay move constructor.
	/*!
	This constructor is of technical utility in wrapping a nuSQuIDSDecay object in
	a nuSQuIDS atmospheric object. See the example scripts.
	*/
	nuSQUIDSDecay(nuSQUIDSDecay&& other):
	nuSQUIDS(std::move(other)),
	ihard_interactions(other.ihard_interactions),
	pscalar(other.pscalar),	majorana(other.majorana), 
	DT(other.DT), DT_evol(other.DT_evol), m_nu(other.m_nu)
	{
		couplings = gsl_matrix_alloc(other.numneu,other.numneu);
		gsl_matrix_memcpy(couplings,other.couplings);
		for (unsigned int chi=0; chi<2; chi++){
			rate_matrices[chi] = gsl_matrix_alloc(other.numneu,other.numneu);
			gsl_matrix_memcpy(rate_matrices[chi],other.rate_matrices[chi]);
		}
	}

	//! nuSQUIDSDecay destructor.
	/*!
	Freeing memory allocated to gsl_matrices.
	*/
	~nuSQUIDSDecay(){
		gsl_matrix_free(couplings);
		for (size_t chi=0; chi<2; chi++){
			gsl_matrix_free(rate_matrices[chi]);
		}
	}

	//! Toggles decay regeneration.
	/*!
		The switch is internal to SQUIDS/nuSQUIDS. If set to true, the
		terms returned by InteractionsRho() are present in the evolution
		equation. If not, there is no regeneration.
		The switch, iinteractions (set in the constructor, also nuSQuIDS internal), 
		toggles additional interactions in both GammaRho() and InteractionsRho(). 
		These interactions are described in the nuSQUIDS documentation 
		under the corresponding function names. If the switch is set to true, 
		these terms are added to the decay Gamma matrix and decay "R" matrix. 
		If set to false, then the only interaction in play is neutrino decay.
		To clarify the behavior of NuSQuIDSDecay in the various 
		iinteractions/DecayRegeneration cases, consider the table below.
		iinteractions | DecayRegeneration | Physics
		--------------- | ----------------- | -------------------------------
		True						| True							| All terms included.
		True						| False						 | All except regeneration terms.
		False					 | True							| Gamma + R (decay with regen only).
		False					 | False						 | Gamma only (decay only without regen).
	\param opt is the boolean value to toggle regeneration.
	*/
	void Set_DecayRegeneration(bool opt) { Set_OtherRhoTerms(opt); }
}; // close nusquids class definition
} // close nusquids namespace
#endif // nusquids_decay_h
