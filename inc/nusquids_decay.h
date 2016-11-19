#ifndef nusquidslv_H
#define nusquidslv_H

#include <vector>
#include <iostream>
#include <nuSQuIDS/nuSQuIDS.h>
#include "exCross.h"


namespace nusquids {

class nuSQUIDSDecay: public nuSQUIDS {
  private:
    bool iincoherent_int = false;
    bool decay_parameters_set = false;
    squids::Const decay_parameters;
    std::vector<double> decay_strength;
    squids::SU_vector DT;
    std::vector<squids::SU_vector> DT_evol;

	std::vector<double> m_nu;
	double m_phi;

	
	double pstar(unsigned int i, unsigned int j) const {
		return (1.0/(2.0*m_nu[i]))*sqrt((pow(m_nu[i],2)-pow(m_nu[j]+m_phi,2))
			*(pow(m_nu[i],2)-pow(m_nu[j]-m_phi,2)));	
	}


	unsigned int nearest_element(double value) const {

		std::vector<double> diffs(E_range.size());
	  	unsigned int i=0;
	  	for (i=0; i<E_range.size(); i++)
		{
				diffs[i] = fabs(value-E_range[i]);
		}
	
	  	return std::distance(diffs.begin(), std::min_element(diffs.begin(),diffs.end()));
	}


  protected:


    void AddToPreDerive(double x) {
      if(!decay_parameters_set)
        throw std::runtime_error("decay parameters not set");
      for(int ei = 0; ei < ne; ei++){
        // asumming same mass hamiltonian for neutrinos/antineutrinos
        squids::SU_vector h0 = H0(E_range[ei],0);
        DT_evol[ei] = DT.Evolve(h0,(x-Get_t_initial()));
      }
    }

    void AddToWriteHDF5(hid_t hdf5_loc_id) const {
      // here we write the new parameters to be saved in the HDF5 file
      hsize_t dim[1]{1};

      // creating arrays to save stuff
      H5LTmake_dataset(hdf5_loc_id,"decay_strength",1,dim,H5T_NATIVE_DOUBLE,0);
      H5LTmake_dataset(hdf5_loc_id,"mixing_angles",1,dim,H5T_NATIVE_DOUBLE,0);
      H5LTmake_dataset(hdf5_loc_id,"CP_phases",1,dim,H5T_NATIVE_DOUBLE,0);

      // save decay strength
      for(size_t i =0; i< numneu; i++){
        std::string decay_strength_label = "lam"+std::to_string(i+1);
        H5LTset_attribute_double(hdf5_loc_id, "decay_strength",decay_strength_label.c_str(),&(decay_strength[i]), 1);
      }

      // save decay mixing angles
      for( unsigned int i = 0; i < numneu; i++ ){
        for( unsigned int j = i+1; j < numneu; j++ ){
          std::string th_label = "th"+std::to_string(i+1)+std::to_string(j+1);
          double th_value = decay_parameters.GetMixingAngle(i,j);
          H5LTset_attribute_double(hdf5_loc_id, "mixing_angles",th_label.c_str(),&th_value, 1);

          std::string delta_label = "delta"+std::to_string(i+1)+std::to_string(j+1);
          double delta_value = decay_parameters.GetPhase(i,j);
          H5LTset_attribute_double(hdf5_loc_id, "CP_phases",delta_label.c_str(),&delta_value, 1);
        }
      }
    }

    void AddToReadHDF5(hid_t hdf5_loc_id){
      // read and set mixing parameters
      for( unsigned int i = 0; i < numneu; i++ ){
        for( unsigned int j = i+1; j < numneu; j++ ){
        double th_value;
        std::string th_label = "th"+std::to_string(i+1)+std::to_string(j+1);
        H5LTget_attribute_double(hdf5_loc_id,"mixing_angles", th_label.c_str(), &th_value);
        decay_parameters.SetMixingAngle(i,j,th_value);

        double delta_value;
        std::string delta_label = "delta"+std::to_string(i+1)+std::to_string(j+1);
        H5LTget_attribute_double(hdf5_loc_id,"CP_phases", delta_label.c_str(), &delta_value);
        decay_parameters.SetPhase(i,j,delta_value);
        }
      }

      // strength vector reset
      decay_strength = std::vector<double>(numneu);
      for(size_t i =0; i< numneu; i++){
        std::string decay_strength_label = "lam"+std::to_string(i+1);
        H5LTget_attribute_double(hdf5_loc_id, "decay_strength",decay_strength_label.c_str(),&(decay_strength[i]));
      }

      Set_Decay_Matrix(decay_parameters,decay_strength);
    }

    squids::SU_vector GammaRho(unsigned int ie,unsigned int irho) const {
      if(iincoherent_int)
        return nuSQUIDS::GammaRho(ie,irho) + DT_evol[ie]*(1./E_range[ie]);
      else
        return DT_evol[ie]*(1./E_range[ie]);
    }

    squids::SU_vector InteractionsRho(unsigned int ie, unsigned int irho) const {
      squids::SU_vector decay_regeneration(numneu);
      // here one needs to fill in the extra decay regeneration terms
			double E0;			
			double E0_index;
			double Ef = E_range[ie];
			double my_pstar;

			unsigned int i,j;
			for(i=0; i<numneu; i++)
			{
				for(j=i+1; i<numneu; j++)
				{
					my_pstar = pstar(i,j);	
					E0 = Ef*m_nu[i]/my_pstar;
					E0_index = nearest_element(E0); 
					decay_regeneration+=(state[E0_index].rho[irho]*evol_b0_proj[irho][i][E0_index])*(((DT_evol[ie])[j+numneu*i])/Ef)*evol_b0_proj[irho][j][ie];
					printmat(DT_evol[ie],numneu,"DTEVOL");
					printmat(evol_b0_proj[irho][i][E0_index],numneu,"MASSPROJ");
					printmat(state[E0_index].rho[irho],numneu,"RHOMATRIX");
					printmat(evol_b0_proj[irho][i][E0_index],numneu,"MASSPROJ:I");
					printmat(evol_b0_proj[irho][j][E0_index],numneu,"MASSPROJ:J");
					printmat(decay_regeneration,numneu,"DCY_REGEN");
					unsigned int k;
					for(k=0; k<numneu; k++)
					{
						std::cout << k << " " << m_nu[k] << " ";
					}
					std::cout << std::endl;
					std::cout << "i,j: " << "(" << i << " , " << j << ")"  << std::endl;
					std::cout << "PHIMASS : " <<m_phi <<std::endl;
					std::cout << "E0 : " << E0 << std::endl;
					std::cout << "E0closest : " << E_range[E0_index] << std::endl;
					std::cout << "E0closestp1 : " << E_range[E0_index+1] << std::endl;
					std::cout << "E0closestm1 : " << E_range[E0_index-1] << std::endl;
					std::cout << "Ef : " << E_range[ie] << std::endl;
					std::cout << "Pstar : " << pstar(i,j) << std::endl;
				}
			}



      //  do not modify after this line
      if(iincoherent_int)
        return nuSQUIDS::GammaRho(ie,irho) + decay_regeneration;
      else
        return decay_regeneration;
    }

  public:
    nuSQUIDSDecay() {};
    nuSQUIDSDecay(marray<double,1> e_nodes,
                  unsigned int numneu_ = 3,NeutrinoType NT_ = NeutrinoType::both,bool iinteraction_ = true):
      nuSQUIDS(e_nodes,numneu_,NT_,iinteraction_,std::make_shared<nusquids::NeutrinoDISCrossSectionsFromTablesExtended>())
    {
      // just allocate some matrices
			DT_evol.resize(ne);
			for(int ei = 0; ei < ne; ei++){
				DT_evol[ei] = squids::SU_vector(nsun);
			}
		
			// allocating space for neutrino masses	
			m_nu.resize(numneu);	

    }

    void Set_Decay_Matrix(squids::Const & decay_parameters,std::vector<double> decay_strength){
      if(decay_strength.size() != numneu){
        throw std::runtime_error("Mismatch size while constructing decay matrix.");
      }
      DT = squids::SU_vector(numneu);
      for(size_t i=0; i<numneu; i++){
       DT += decay_strength[i]*squids::SU_vector::Projector(numneu,i);
      }
      // rotate from decay basis to mass basis
      DT.RotateToB1(decay_parameters);
      decay_parameters_set = true;
    }

    void Set_IncoherentInteractions(bool opt){
      iincoherent_int = opt;
    }

	void Set_m_nu(double mass, unsigned int state)
	{
		m_nu[state]=mass;
	}

	void Set_m_phi(double mass)
	{
		m_phi=mass;
	}


	

	void printmat(const squids::SU_vector mat, unsigned int dim, std::string mname) const
	{
		std::cout << std::endl;
		std::cout << "Matrix: " << mname << std::endl;

		unsigned int i,j;
	
		for(i=0; i<dim; i++)
		{
			for(j=0; j<dim; j++)
			{
				std::cout << (mat)[j+dim*i] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}



};

} // close nusquids namespace

#endif //nusquidslv_h

