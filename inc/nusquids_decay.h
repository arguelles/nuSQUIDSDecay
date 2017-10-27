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

FIXME add reference to paper??

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
	//------------------------Data------------------------//
	bool majorana;
	bool iincoherent_int;
	enum{SCALAR, PSEUDOSCALAR};
	//Chirality Preserving Process or Chirality Violating Process
	enum{CPP,CVP};
	//One for each of {CPP,CVP}x{SCALAR,PSEUDOSCALAR}
	gsl_matrix* rate_matrices[2][2];
	//One scalar g_ij, one pseudoscalar.
	gsl_matrix* couplings[2];
	std::vector<double> m_nu;
	double m_phi;
	squids::SU_vector DT;
	std::vector<squids::SU_vector> DT_evol;

	//----------------------Functions---------------------//
	//Trying to keep the "model-specific" functions in private,
	//and have moved all functions related to more general decay models
	//to protected.

	//Define auxiliary functions (4a),(4b),(4c) in [1]
	double f(double x){	
		double result =  x/2.0 + 2.0 + (2.0/x)*log(x) - 2/(x*x) - 1.0/(2.0*x*x*x);
		return result;
	double g(double x){	
		double result =  x/2.0 - 2.0 + (2.0/x)*log(x) + 2/(x*x) - 1.0/(2.0*x*x*x);
		return result;
	double k(double x){	
		double result =  x/2.0 - (2.0/x)*log(x) - 1.0/(2.0*x*x*x);
		return result;

	void Compute_Rate_Matrices(){
		//Compute *rest frame* decay rate matrices.
		
		//CPP,SCALAR	
		rate_matrices[CPP][SCALAR] = gsl_matrix_alloc(numneu,numneu); 
		for (unsigned int i=0; i<numneu; i++){
			for (unsigned int j=0; j<numneu; j++){
				double x_ij = m_nu[i]/m_nu[j]; 
				double g_ij = gsl_matrix_get(couplings[SCALAR],i,j);
				double rate = (m_nu[i]/(16.0*M_PI))*(1.0/x_ij)*g_ij*g_ij*f(x_ij);
				gsl_matrix_set(rate_matrices[CPP][SCALAR],i,j,rate);
			}
		}
		//CPP,PSEUDOSCALAR	
		rate_matrices[CPP][PSEUDOSCALAR] = gsl_matrix_alloc(numneu,numneu); 
		for (unsigned int i=0; i<numneu; i++){
			for (unsigned int j=0; j<numneu; j++){
				double x_ij = m_nu[i]/m_nu[j]; 
				double g_ij = gsl_matrix_get(couplings[PSEUDOSCALAR],i,j);
				double rate = (m_nu[i]/(16.0*M_PI))*(1.0/x_ij)*g_ij*g_ij*g(x_ij);
				gsl_matrix_set(rate_matrices[CPP][PSEUDOSCALAR],i,j,rate);
			}
		}
		//CVP,SCALAR	
		rate_matrices[CVP][SCALAR] = gsl_matrix_alloc(numneu,numneu); 
		for (unsigned int i=0; i<numneu; i++){
			for (unsigned int j=0; j<numneu; j++){
				double x_ij = m_nu[i]/m_nu[j]; 
				double g_ij = gsl_matrix_get(couplings[SCALAR],i,j);
				double rate = (m_nu[i]/(16.0*M_PI))*(1.0/x_ij)*g_ij*g_ij*k(x_ij);
				gsl_matrix_set(rate_matrices[CVP][SCALAR],i,j,rate);
			}
		}
		//CVP,PSEUDOSCALAR	
		rate_matrices[CVP][PSEUDOSCALAR] = gsl_matrix_alloc(numneu,numneu); 
		for (unsigned int i=0; i<numneu; i++){
			for (unsigned int j=0; j<numneu; j++){
				double x_ij = m_nu[i]/m_nu[j]; 
				double g_ij = gsl_matrix_get(couplings[PSEUDOSCALAR],i,j);
				double rate = (m_nu[i]/(16.0*M_PI))*(1.0/x_ij)*g_ij*g_ij*k(x_ij);
				gsl_matrix_set(rate_matrices[CVP][PSEUDOSCALAR],i,j,rate);
			}
		}
	}

	squids::SU_vector InteractionsRho(unsigned int iedaughter, unsigned int irho) const {
		squids::SU_vector decay_regeneration(numneu);
		// here one needs to fill in the extra decay regeneration terms
		double edaughter = E_range[iedaughter];
		// j-daughter index
		for (size_t j = 0; j < numneu; j++) {
			// i-parent index
			for (size_t i = j+1; i < numneu; i++) {
				//parent-to-daughter mass ratio
				double xij = m_nu[i]/m_nu[j];
				//boost factor to lab frame
				double ieparent_high = nearest_element(edaughter*xij);
				// i-energy (parent energy) index 
				//left-rectangular integral approximation
				for (size_t ieparent = iedaughter; ieparent < ieparent_high-1; ieparent++) {
					double eparent = E_range[ieparent];
					double gamma = eparent/m_nu[i];
					double delta_eparent = E_range[ieparent+1]-E_range[ieparent];
					//combining scalar and pseudoscalar contributions
					decay_regeneration += (delta_eparent)*(state[ieparent].rho[irho]
											*evol_b0_proj[irho][i][ieparent])*
											(1/(eparent*eparent*edaughter))*
											((gsl_matrix_get(rate_matrices[CPP][SCALAR],j,i)/gamma)*
											pow(eparent+xij*edaughter,2)/pow(xij+1,2)+
											(gsl_matrix_get(rate_matrices[CPP][PSEUDOSCALAR],j,i)/gamma)*
											pow(eparent-xij*edaughter,2)/pow(xij-1,2))*
											(evol_b0_proj[irho][j][iedaughter]);
				}

				//Include chirality-violating term if neutrino is majorana.
				if(majorana){
					unsigned int parent_irho;
					if (irho==0) parent_irho=1;
					if (irho==1) parent_irho=0;
					// i-energy (parent energy) index 
					//left-rectangular integral approximation
					for (size_t ieparent = iedaughter; ieparent < ieparent_high-1; ieparent++) {
						double eparent = E_range[ieparent];
						double gamma = eparent/m_nu[i];
						double delta_eparent = E_range[ieparent+1]-E_range[ieparent];
						//combining scalar and pseudoscalar contributions
						decay_regeneration += (delta_eparent)*(state[ieparent].rho[parent_irho]
												*evol_b0_proj[parent_irho][i][ieparent])*
												((eparent-edaughter)/(eparent*eparent*edaughter))*
												((gsl_matrix_get(rate_matrices[CVP][SCALAR],j,i)/gamma)*
												(edaughter*pow(xij,2)-eparent)/pow(xij+1,2)+
												(gsl_matrix_get(rate_matrices[CVP][PSEUDOSCALAR],j,i)/gamma)*
												(edaughter*pow(xij,2)-eparent)/pow(xij-1,2))*
												(evol_b0_proj[irho][j][iedaughter]);
					}
				}
			}
		}
		//	do not modify after this line
		if (iincoherent_int)
			return nuSQUIDS::InteractionsRho(iedaughter, irho) + decay_regeneration;
		else
			return decay_regeneration;
	}

protected:
	unsigned int nearest_element(double value) const {
		std::vector<double> diffs(E_range.size());
		for (size_t i = 0; i < E_range.size(); i++) {
			diffs[i] = fabs(value - E_range[i]);
		}
		return std::distance(diffs.begin(),std::min_element(diffs.begin(), diffs.end()));
	}
	
	//FIXME: REMOVE IN RELEASE
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

	void Check_Matrix_Size(gsl_matrix* m1, gsl_matrix* m2){
		if (m1->size1 != m2->size1){
			throw std::runtime_error("size1 mismatch while copying matrix.");
		}
		if (m1->size2 != m2->size2){
			throw std::runtime_error("size2 mismatch while copying matrix.");
		}
	}
	void Set_Rate_Matrices(gsl_matrix  *rate_matrices_[2][2]){
		for (unsigned int chi=0; chi<2; chi++){
			for (unsigned int s=0; s<2; s++){			
				rate_matrices[chi][s] = gsl_matrix_alloc(numneu,numneu);
				Check_Matrix_Size(rate_matrices[chi][s],rate_matrices_[chi][s]);
				gsl_matrix_memcpy(rate_matrices[chi][s],rate_matrices_[chi][s]);
			}
		}
	}

	void Compute_DT(){
		DT = squids::SU_vector(numneu);
		for(size_t i = 0; i < numneu; i++){
			double rate=0;
			for (size_t chi=0; chi<2; chi++){
				for (size_t s=0; s<2; s++){			
					rate+=gsl_matrix_get(rate_matrices[chi][s],i,i);	
				}
			}
			DT += m_nu[i]*rate*squids::SU_vector::Projector(numneu, i);
		}
	}

	void AddToPreDerive(double x) {
		for (int ei = 0; ei < ne; ei++) {
			// asumming same mass hamiltonian for neutrinos/antineutrinos
			squids::SU_vector h0 = H0(E_range[ei], 0);
			DT_evol[ei] = DT.Evolve(h0, (x - Get_t_initial()));
		}
	}

	squids::SU_vector GammaRho(unsigned int ie, unsigned int irho) const {
		if (iincoherent_int)
			return nuSQUIDS::GammaRho(ie, irho) + DT_evol[ie] * (0.5 / E_range[ie]);
		else
			return DT_evol[ie] * (0.5 / E_range[ie]);
	}


public:
	//Preferred constructor (foolproof)
	nuSQUIDSDecay(marray<double, 1> e_nodes, unsigned int numneu_,
					NeutrinoType NT_, bool iinteraction_,
					std::vector<double> m_nu_, double m_phi_,
					gsl_matrix *couplings_[], 
					bool iincoherent_int, bool majorana_):
					nuSQUIDSDecay(e_nodes,numneu_,NT_,iinteraction_){
		m_nu=m_nu_;
		m_phi=m_phi_;
		Set_Couplings(couplings_);
		Compute_Rate_Matrices();
		Compute_DT();

		iincoherent_int=iincoherent_int_;
		majorana=majorana_;
	}

	//Our constructor (allowing for NuSQuIDS Atmospheric Wrapping)
	nuSQUIDSDecay(marray<double, 1> e_nodes, unsigned int numneu_,
					NeutrinoType NT_, bool iinteraction_,
					std::vector<double> m_nu_, double m_phi_,
					gsl_matrix *rate_matrices_[][]):
					nuSQUIDSDecay(e_nodes,numneu_,NT_,iinteraction_){
		m_nu=m_nu_;
		m_phi=m_phi_;
		Set_Rate_Matrices(rate_matrices_);
		Compute_DT();

		//Default values
		iincoherent_int=false;
		majorana=false;
	}

	//Move constructor
	nuSQUIDSDecay(nuSQUIDSDecay&& other):
	nuSQUIDS(std::move(other)), numneu(other.numneu),
	iincoherent_int(other.iincoherent_int),
	majorana(other.majorana), DT(other.DT),
	DT_evol(other.DT_evol), m_nu(other.m_nu),
	m_phi(other.m_phi), rate_matrices(other.rate_matrices),
	couplings(other.couplings)
	{
		for (unsigned int s=0; s<2; s++){			
			couplings[s] = gsl_matrix_alloc(numneu,numneu);
			gsl_matrix_memcpy(couplings[s],other.couplings[s]);
			for (unsigned int chi=0; chi<2; chi++){
				rate_matrices[chi][s] = gsl_matrix_alloc(numneu,numneu);
				gsl_matrix_memcpy(rate_matrices[chi][s],other.rate_matrices[chi][s]);
			}
		}
	}

	//Destructor
	~nuSQUIDSDecay(){
		for (size_t s=0; s<2; s++){			
			gsl_matrix_free(couplings[s]);
			for (size_t chi=0; chi<2; chi++){
				gsl_matrix_free(rate_matrices[chi][s]);
			}
		}
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
		True						| True							| All terms included. 
		True						| False						 | All except regeneration terms.
		False					 | True							| Gamma + R (decay with regen only).
		False					 | False						 | Gamma only (decay only without regen).
	\param opt: the boolean value to toggle regeneration.
	*/
	void SetDecayRegeneration(bool opt) { Set_OtherRhoTerms(opt); }




	void Set_Majorana(bool opt) { majorana = opt; }
	void Set_m_nu(double mass, unsigned int state) { m_nu[state] = mass; }
	void Set_m_phi(double mass) { m_phi = mass; }

	void Set_Couplings(gsl_matrix  *couplings_[2]){
		for (unsigned int s=0; s<2; s++){
			couplings[s] = gsl_matrix_alloc(numneu,numneu);
			Check_Matrix_Size(couplings[s],couplings_[s]);
			gsl_matrix_memcpy(couplings[s],couplings_[s]);
		}
	}

}; // close nusquids class definition
} // close nusquids namespace
#endif // nusquidslv_h
