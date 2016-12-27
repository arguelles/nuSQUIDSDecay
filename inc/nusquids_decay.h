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

bool close_enough(double x, double y)
{
    double epsilon = std::numeric_limits<double>::epsilon();
    return (std::abs(x - y) <= epsilon * std::abs(x)
        && std::abs(x - y) <= epsilon * std::abs(y));
}

namespace nusquids {

class nuSQUIDSDecay : public nuSQUIDS {
private:
  bool iincoherent_int = false;
  bool decay_parameters_set = false;
  squids::Const decay_parameters;
  std::vector<double> decay_strength;
  squids::SU_vector DT;
  std::vector<squids::SU_vector> DT_evol;

  gsl_matrix* rate_mat;

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
    if (!decay_parameters_set)
      throw std::runtime_error("decay parameters not set");
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
	  //std::cout << "IE: " << ie << "  E: " << E_range[ie] << std::endl; 
	  //printmat(DT_evol[ie] * (0.5 / E_range[ie]),4,"GAMMA");
	  //std::cout << std::endl;
      return DT_evol[ie] * (0.5 / E_range[ie]);
  }

  squids::SU_vector InteractionsRho(unsigned int ie, unsigned int irho) const {
    squids::SU_vector decay_regeneration(numneu);
    // here one needs to fill in the extra decay regeneration terms
    double Ef = E_range[ie];
    // auxiliary variables
    double E0;
    size_t E0_index;
    double my_pstar;
	double gamma;

    //printmat(decay_regeneration, numneu, "DCY_REGEN");
    // i-daughter index
    for (size_t i = 0; i < numneu; i++) {
      // j-parent index
      for (size_t j = i + 1; j < numneu; j++) {
        my_pstar = pstar(j,i);
		gamma = Ef/my_pstar;
        E0 = m_nu[j]*gamma;
        E0_index = nearest_element(E0);
        decay_regeneration +=
            (state[E0_index].rho[irho]*evol_b0_proj[irho][j][E0_index])*
            (gsl_matrix_get(rate_mat,i,j) / gamma) * evol_b0_proj[irho][i][ie];


        //---------Testing trace functionality---------//

		/*
        std::vector<double> rhocomps =
            (state[E0_index].rho[irho]).GetComponents();

        std::cout << "rhocomps: " << std::endl;
        unsigned int l;
        for (l = 0; l < rhocomps.size(); l++) {
          std::cout << rhocomps[l] << " ";
        }
        std::cout << std::endl;

        std::vector<double> projcomps =
            (evol_b0_proj[irho][i][E0_index]).GetComponents();

        std::cout << "projcomps: " << std::endl;
        for (l = 0; l < projcomps.size(); l++) {
          std::cout << projcomps[l] << " ";
        }
        std::cout << std::endl;

        unsigned int k;
        for (k = 0; k < numneu; k++) {
          std::cout << k << " " << m_nu[k] << " ";
        }
        std::cout << std::endl;
        std::cout << "i,j: "
                  << "(" << i << " , " << j << ")" << std::endl;
        std::cout << "PHIMASS : " << m_phi << std::endl;
        std::cout << "E0 : " << E0 << std::endl;
        std::cout << "E0closest : " << E_range[E0_index] << std::endl;
        std::cout << "E0closestp1 : " << E_range[E0_index + 1] << std::endl;
        std::cout << "E0closestm1 : " << E_range[E0_index - 1] << std::endl;
        std::cout << "Ef : " << E_range[ie] << std::endl;
        //std::cout << "Pstar : " << pstar(i, j) << std::endl;
        printmat(DT_evol[ie], numneu, "DTEVOL");
        printmat(state[E0_index].rho[irho], numneu, "RHOMATRIX");
        printmat(evol_b0_proj[irho][i][E0_index], numneu, "MASSPROJ:I");
        printmat(evol_b0_proj[irho][j][E0_index], numneu, "MASSPROJ:J");
        printmat(decay_regeneration, numneu, "DCY_REGEN");
		*/
      }
    }

    //  do not modify after this line
    if (iincoherent_int)
      return nuSQUIDS::GammaRho(ie, irho) + decay_regeneration;
    else
      return decay_regeneration;
  }

public:
  nuSQUIDSDecay(){};
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
    m_nu.resize(numneu);


	// allocating memory for rate matrix
	rate_mat = gsl_matrix_alloc(numneu,numneu);

  }

  ~nuSQUIDSDecay()
  {
    gsl_matrix_free(rate_mat);	
  }

  void Set_Decay_Matrix(gsl_matrix* m) 
  {
	if (rate_mat->size1 != m->size1)
	{
		throw std::runtime_error("size1 mismatch while constructing decay matrix.");
	}

	if (rate_mat->size2 != m->size2)
	{
		throw std::runtime_error("size2 mismatch while constructing decay matrix.");
	}

	/*
	double colsum;
	for (size_t col=0; col<m->size2; col++)
	{
	    colsum=0;

	    for (size_t row=0; row<col; row++)
	    {
	        colsum+=gsl_matrix_get(m,row,col);
	    }

		if (!close_enough(gsl_matrix_get(m,col,col),colsum))
		{
			throw std::runtime_error("Off-diagonal decay rates do not sum to diagonal rates.");
		}
	}
	*/

	gsl_matrix_memcpy(rate_mat,m);
    DT = squids::SU_vector(numneu);

	double entry;
    for (size_t i = 0; i < numneu; i++) 
	{
		entry = gsl_matrix_get(rate_mat,i,i);
    	DT += entry*squids::SU_vector::Projector(numneu, i);
    }

	decay_parameters_set=true;
  }

  void Set_IncoherentInteractions(bool opt) { iincoherent_int = opt; }

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
