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

class nuSQUIDSDecay : public nuSQUIDS {
private:
  bool iincoherent_int = false;
  bool DecayParametersSet = false;
  bool NeutrinoMassesSet = false;
  bool PhiMassSet = false;

	//decay parameters: clean up! Remove from HDF5?
  squids::Const decay_parameters;
  std::vector<double> decay_strength;
  squids::SU_vector DT;
  std::vector<squids::SU_vector> DT_evol;

  gsl_matrix* rate_mat;
	gsl_matrix* pstar_mat;

  std::vector<double> NeutrinoMasses;
  double PhiMass;

  double pstar(unsigned int daughter, unsigned int parent) const {
    if (NeutrinoMasses[parent] < NeutrinoMasses[daughter] + PhiMass) {
      throw std::runtime_error("non physical case");
    }
    double retval = (1.0 / (2.0 * NeutrinoMasses[daughter])) *
                    sqrt((pow(NeutrinoMasses[daughter], 2) - pow(NeutrinoMasses[parent] + PhiMass, 2)) *
                         (pow(NeutrinoMasses[daughter], 2) - pow(NeutrinoMasses[parent] - PhiMass, 2)));
    return retval;
  }

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

	//FIXME
    //Set_Decay_Matrix(decay_parameters, decay_strength);
  }

  squids::SU_vector GammaRho(unsigned int ie, unsigned int irho) const {
    if (iincoherent_int)
      return nuSQUIDS::GammaRho(ie, irho) + DT_evol[ie] * (0.5 / E_range[ie]);
    else
      return DT_evol[ie] * (0.5 / E_range[ie]);
  }

  squids::SU_vector InteractionsRho(unsigned int ie, unsigned int irho) const {
    squids::SU_vector decay_regeneration(numneu);
    double Ef = E_range[ie];
    // i-daughter index
    for (size_t i = 0; i < numneu; i++) {
      // j-parent index
      for (size_t j = i + 1; j < numneu; j++) {
				std::cout << "Matrx Pstar: " << gsl_matrix_get(pstar_mat,i,j) << std::endl;
        double gamma = Ef/gsl_matrix_get(pstar_mat,i,j);
        double E0 = NeutrinoMasses[j]*gamma;
        size_t E0_index = nearest_element(E0);
        decay_regeneration += (state[E0_index].rho[irho]
					*evol_b0_proj[irho][j][E0_index])
					*(gsl_matrix_get(rate_mat,i,j)/gamma)
					*(evol_b0_proj[irho][i][ie]);
      }
    }

    if (iincoherent_int)
      return nuSQUIDS::InteractionsRho(ie, irho) + decay_regeneration;
    else
      return decay_regeneration;
  }

public:
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

  ~nuSQUIDSDecay(){
    gsl_matrix_free(rate_mat);
    gsl_matrix_free(pstar_mat);
  }

  void SetDecayMatrix(gsl_matrix* tau){
    if (rate_mat->size1 != tau->size1){
      throw std::runtime_error("size1 mismatch while constructing decay matrix.");
    }
    if (rate_mat->size2 != tau->size2){
      throw std::runtime_error("size2 mismatch while constructing decay matrix.");
    }

	  gsl_matrix_set_zero(rate_mat);
	  for (size_t col=0; col<numneu; col++)
	  {
	    double colrate=0;
	    for (size_t row=0; row<col; row++)
	    {
	      double rate = 1.0/gsl_matrix_get(tau,row,col);
	      gsl_matrix_set(rate_mat,row,col,rate);
	      colrate+=rate*NeutrinoMasses[col];
	    }
	    gsl_matrix_set(rate_mat,col,col,colrate);
		}

    DT = squids::SU_vector(numneu);
    for(size_t i = 0; i < numneu; i++){
      double entry = gsl_matrix_get(rate_mat,i,i);
      DT += entry*squids::SU_vector::Projector(numneu, i);
    }

    DecayParametersSet=true;
  }


	gsl_matrix* GetDecayMatrix(void){
		return rate_mat;
	}

  void SetIncoherentInteractions(bool opt) { iincoherent_int = opt; }

	void SetNeutrinoMasses(double lightmass){
		NeutrinoMasses[0] = lightmass;
		for (unsigned int i=1; i<numneu; i++){
			NeutrinoMasses[i] = sqrt(Get_SquareMassDifference(i) + NeutrinoMasses[0]*NeutrinoMasses[0]);
		}
		NeutrinoMassesSet=true;
		PreComputePstarMat();
	}
		

	double GetNeutrinoMass(unsigned int index){
    if (index>=numneu){
      throw std::runtime_error("Index larger than numneu-1");
    }
		return NeutrinoMasses[index];
	}

	void SetPhiMass(double PhiMass_){
		PhiMass = PhiMass_;
		PhiMassSet=true;
	}

	double GetPhiMass(void){
		return PhiMass;
	}


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


  void printgslmat(const gsl_matrix* mat, unsigned int dim,
                std::string mname) const {
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
};

} // close nusquids namespace

#endif // nusquidslv_h
