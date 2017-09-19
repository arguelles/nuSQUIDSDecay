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
#include <math.h>

bool close_enough(double x, double y)
{
    double epsilon = std::numeric_limits<double>::epsilon();
    return (std::abs(x - y) <= epsilon * std::abs(x)
        && std::abs(x - y) <= epsilon * std::abs(y));
}

namespace nusquids {

class nuSQUIDSDecay : public nuSQUIDS {
private:
	bool majorana = false;
  bool iincoherent_int = false;
  bool scalar_decay_parameters_set = false;
  bool pseudoscalar_decay_parameters_set = false;
  bool dt_set = false;
  squids::Const scalar_decay_parameters;
  squids::Const pseudoscalar_decay_parameters;
  std::vector<double> scalar_decay_strength;
  std::vector<double> pseudoscalar_decay_strength;
  squids::SU_vector DT;
  std::vector<squids::SU_vector> DT_evol;

  gsl_matrix* scalar_decay_mat;
  gsl_matrix* pseudoscalar_decay_mat;

  std::vector<double> m_nu;
  double m_phi;


  double pstar(unsigned int i, unsigned int j) const {
    if (m_nu[i] < m_nu[j] + m_phi) {
		//std::cout << "I,J: " << i << " " << j << std::endl;
		//std::cout << "MI,MJ: " << m_nu[i] << " " << m_nu[j] << std::endl;
		//std::cout << "MPHI: " << m_phi << std::endl;
      throw std::runtime_error("non physical case");
    }
    double retval = (1.0 / (2.0 * m_nu[i])) *
                    sqrt((pow(m_nu[i], 2) - pow(m_nu[j] + m_phi, 2)) *
                         (pow(m_nu[i], 2) - pow(m_nu[j] - m_phi, 2)));
    return retval;
  }

  unsigned int nearest_element(double value) const {

    std::vector<double> diffs(E_range.size());
    for (size_t i = 0; i < E_range.size(); i++) {
      diffs[i] = fabs(value - E_range[i]);
    }

    return std::distance(diffs.begin(),
                         std::min_element(diffs.begin(), diffs.end()));
  }

protected:
  void AddToPreDerive(double x) {
    if (!(scalar_decay_parameters_set&&pseudoscalar_decay_parameters_set&&dt_set))
      throw std::runtime_error("decay parameters not set");
    for (int ei = 0; ei < ne; ei++) {
      // asumming same mass hamiltonian for neutrinos/antineutrinos
      squids::SU_vector h0 = H0(E_range[ei], 0);
      DT_evol[ei] = DT.Evolve(h0, (x - Get_t_initial()));
    }
  }

  /*
  void AddToWriteHDF5(hid_t hdf5_loc_id) const {
    // here we write the new parameters to be saved in the HDF5 file
    hsize_t dim[1]{1};

    // creating arrays to save stuff
    H5LTmake_dataset(hdf5_loc_id, "decay_strength", 1, dim, H5T_NATIVE_DOUBLE,
                     0);
    H5LTmake_dataset(hdf5_loc_id, "mixing_angles", 1, dim, H5T_NATIVE_DOUBLE,
                     0);
    H5LTmake_dataset(hdf5_loc_id, "CP_phases", 1, dim, H5T_NATIVE_DOUBLE, 0);

    // save decay strength
    for (size_t i = 0; i < numneu; i++) {
      std::string decay_strength_label = "lam" + std::to_string(i + 1);
      H5LTset_attribute_double(hdf5_loc_id, "decay_strength",
                               decay_strength_label.c_str(),
                               &(decay_strength[i]), 1);
    }

    // save decay mixing angles
    for (unsigned int i = 0; i < numneu; i++) {
      for (unsigned int j = i + 1; j < numneu; j++) {
        std::string th_label =
            "th" + std::to_string(i + 1) + std::to_string(j + 1);
        double th_value = decay_parameters.GetMixingAngle(i, j);
        H5LTset_attribute_double(hdf5_loc_id, "mixing_angles", th_label.c_str(),
                                 &th_value, 1);

        std::string delta_label =
            "delta" + std::to_string(i + 1) + std::to_string(j + 1);
        double delta_value = decay_parameters.GetPhase(i, j);
        H5LTset_attribute_double(hdf5_loc_id, "CP_phases", delta_label.c_str(),
                                 &delta_value, 1);
      }
    }
  }

  void AddToReadHDF5(hid_t hdf5_loc_id) {
    // read and set mixing parameters
    for (unsigned int i = 0; i < numneu; i++) {
      for (unsigned int j = i + 1; j < numneu; j++) {
        double th_value;
        std::string th_label =
            "th" + std::to_string(i + 1) + std::to_string(j + 1);
        H5LTget_attribute_double(hdf5_loc_id, "mixing_angles", th_label.c_str(),
                                 &th_value);
        decay_parameters.SetMixingAngle(i, j, th_value);

        double delta_value;
        std::string delta_label =
            "delta" + std::to_string(i + 1) + std::to_string(j + 1);
        H5LTget_attribute_double(hdf5_loc_id, "CP_phases", delta_label.c_str(),
                                 &delta_value);
        decay_parameters.SetPhase(i, j, delta_value);
      }
    }

    // strength vector reset
    decay_strength = std::vector<double>(numneu);
    for (size_t i = 0; i < numneu; i++) {
      std::string decay_strength_label = "lam" + std::to_string(i + 1);
      H5LTget_attribute_double(hdf5_loc_id, "decay_strength",
                               decay_strength_label.c_str(),
                               &(decay_strength[i]));
    }
  }
  */

  squids::SU_vector GammaRho(unsigned int ie, unsigned int irho) const {
    if (iincoherent_int)
      return nuSQUIDS::GammaRho(ie, irho) + DT_evol[ie] * (0.5 / E_range[ie]);
    else
	  //std::cout << "IE: " << ie << "  E: " << E_range[ie] << std::endl;
	  //printmat(DT_evol[ie] * (0.5 / E_range[ie]),4,"GAMMA");
	  //std::cout << std::endl;
      return DT_evol[ie] * (0.5 / E_range[ie]);
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
											((gsl_matrix_get(scalar_decay_mat,j,i)/gamma)*
											pow(eparent+xij*edaughter,2)/pow(xij+1,2)+
											(gsl_matrix_get(pseudoscalar_decay_mat,j,i)/gamma)*
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
												((gsl_matrix_get(scalar_decay_mat,j,i)/gamma)*
												(edaughter*pow(xij,2)-eparent)/pow(xij+1,2)+
												(gsl_matrix_get(pseudoscalar_decay_mat,j,i)/gamma)*
												(edaughter*pow(xij,2)-eparent)/pow(xij-1,2))*
	                      (evol_b0_proj[irho][j][iedaughter]);
	    		}
				}
			}
		}
    //  do not modify after this line
    if (iincoherent_int)
      return nuSQUIDS::InteractionsRho(iedaughter, irho) + decay_regeneration;
    else
      return decay_regeneration;
  }

public:
  nuSQUIDSDecay(marray<double, 1> e_nodes, unsigned int numneu_ = 3,
                NeutrinoType NT_ = NeutrinoType::both,
                bool iinteraction_ = true)
      : nuSQUIDS(e_nodes, numneu_, NT_, iinteraction_){
    //             std::make_shared<
    //                 nusquids::NeutrinoDISCrossSectionsFromTablesExtended>()) {
    // just allocate some matrices
    DT_evol.resize(ne);
    for (int ei = 0; ei < ne; ei++) {
      DT_evol[ei] = squids::SU_vector(nsun);
    }

    // allocating space for neutrino masses
    m_nu.resize(numneu);

    // allocating memory for rate matrix
    scalar_decay_mat = gsl_matrix_alloc(numneu,numneu);
    pseudoscalar_decay_mat = gsl_matrix_alloc(numneu,numneu);
  }

  nuSQUIDSDecay(marray<double, 1> e_nodes, unsigned int numneu_,
                NeutrinoType NT_, bool iinteraction_,
                gsl_matrix * scalar_decay_matrix_, 
                gsl_matrix * pseudoscalar_decay_matrix_, 
								std::vector<double> m_nu_, double m_phi_):
                nuSQUIDSDecay(e_nodes,numneu_,NT_,iinteraction_){
      m_nu=m_nu_;
      m_phi=m_phi_;
    	scalar_decay_mat = gsl_matrix_alloc(numneu,numneu);
    	pseudoscalar_decay_mat = gsl_matrix_alloc(numneu,numneu);
      Set_Scalar_Matrix(scalar_decay_matrix_);
      Set_Pseudoscalar_Matrix(pseudoscalar_decay_matrix_);
	  Compute_DT();
      iincoherent_int=true;
  }

  ~nuSQUIDSDecay(){
    gsl_matrix_free(scalar_decay_mat);
    gsl_matrix_free(pseudoscalar_decay_mat);
  }

  // move constructor
  nuSQUIDSDecay(nuSQUIDSDecay&& other):
  nuSQUIDS(std::move(other)),
  iincoherent_int(other.iincoherent_int),
  majorana(other.majorana),
  scalar_decay_parameters_set(other.scalar_decay_parameters_set),
  scalar_decay_parameters(std::move(other.scalar_decay_parameters)),
  pseudoscalar_decay_parameters_set(other.pseudoscalar_decay_parameters_set),
  pseudoscalar_decay_parameters(std::move(other.pseudoscalar_decay_parameters)),
  scalar_decay_strength(other.scalar_decay_strength),
  pseudoscalar_decay_strength(other.pseudoscalar_decay_strength),
  DT(other.DT),
	dt_set(std::move(other.dt_set)),
  DT_evol(other.DT_evol),
  m_nu(other.m_nu),
  m_phi(other.m_phi)
  {
    if(other.scalar_decay_parameters_set&&other.pseudoscalar_decay_parameters_set){
      scalar_decay_mat = gsl_matrix_alloc(other.numneu,other.numneu);
      pseudoscalar_decay_mat = gsl_matrix_alloc(other.numneu,other.numneu);
      gsl_matrix_memcpy(scalar_decay_mat,other.scalar_decay_mat);
      gsl_matrix_memcpy(pseudoscalar_decay_mat,other.pseudoscalar_decay_mat);
    }
  }

  void Set_Scalar_Matrix(gsl_matrix* m){

    if (scalar_decay_mat->size1 != m->size1){
      throw std::runtime_error("size1 mismatch while constructing decay matrix.");
    }
    if (scalar_decay_mat->size2 != m->size2){
      throw std::runtime_error("size2 mismatch while constructing decay matrix.");
    }

    gsl_matrix_memcpy(scalar_decay_mat,m);

    scalar_decay_parameters_set=true;
  }

  void Set_Pseudoscalar_Matrix(gsl_matrix* m){

    if (pseudoscalar_decay_mat->size1 != m->size1){
      throw std::runtime_error("size1 mismatch while constructing decay matrix.");
    }
    if (pseudoscalar_decay_mat->size2 != m->size2){
      throw std::runtime_error("size2 mismatch while constructing decay matrix.");
    }

    gsl_matrix_memcpy(pseudoscalar_decay_mat,m);

    pseudoscalar_decay_parameters_set=true;
  }

	void Compute_DT(){
    DT = squids::SU_vector(numneu);
    for(size_t i = 0; i < numneu; i++){
      double entry = gsl_matrix_get(scalar_decay_mat,i,i)+gsl_matrix_get(pseudoscalar_decay_mat,i,i);
      DT += entry*squids::SU_vector::Projector(numneu, i);
    }
		dt_set=true;
	}


  void Set_IncoherentInteractions(bool opt) { iincoherent_int = opt; }

  void Set_Majorana(bool opt) { majorana = opt; }

  void Set_m_nu(double mass, unsigned int state) { m_nu[state] = mass; }

  void Set_m_phi(double mass) { m_phi = mass; }

  void printmat(const squids::SU_vector mat, unsigned int dim,
                std::string mname) const {
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
};

} // close nusquids namespace

#endif // nusquidslv_h
