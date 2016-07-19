#ifndef nusquidslv_H
#define nusquidslv_H

#include <vector>
#include <iostream>
#include <nuSQuIDS/nuSQuIDS.h>


namespace nusquids {

class nuSQUIDSDecay: public nuSQUIDS {
  private:
    bool decay_parameters_set = false;
    squids::Const decay_parameters;
    std::vector<double> decay_strength;
    squids::SU_vector DT;
    std::vector<squids::SU_vector> DT_evol;
  protected:
    void AddToPreDerive(double x) {
      if(!lv_parameters_set)
        throw std::runtime_error("LV parameters not set");
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
        decay_parameters.Set_MixingAngle(i,j,th_value);

        double delta_value;
        std::string delta_label = "delta"+std::to_string(i+1)+std::to_string(j+1);
        H5LTget_attribute_double(hdf5_loc_id,"CP_phases", delta_label.c_str(), &delta_value);
        decay_parameters.Set_CPPhase(i,j,delta_value);
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

    squids:: SU_vector GammaRho(unsigned int ie,unsigned int irho) const {
      return nuSQUIDS::GammaRho(ie,irho) + DT_evol[ie];
    }

  public:
    nuSQUIDSDecay() {};
    nuSQUIDSDecay(double Emin_,double Emax_,int Esize_,int numneu_,NeutrinoType NT_,
         bool elogscale_,bool iinteraction_):
          nuSQUIDS(Emin_,Emax_,Esize_,numneu_,NT_,elogscale_,iinteraction_)
    {
      // just allocate some matrices
       DT_evol.resize(ne);
       for(int ei = 0; ei < ne; ei++){
         DT_evol[ei] = squids::SU_vector(nsun);
       }
    }

    void Set_Decay_Matrix(squids::Const & decay_parameters,std::vector<double> decay_strength){
      if(decay_strength.size() != numneu){
        throw std::runtime_error("Mismatch size while constructing decay matrix.");
      }
      DT = squids::SU_vector(numneu);
      for(size_t i=0; i<numneu; i++){
       DT += decay_strength[i]*squids::Projector(numneu,i);
      }
      // rotate from decay basis to mass basis
      DT.RotateToB1(decay_parameters);
      decay_parameters_set = true;
    }
};

} // close nusquids namespace

#endif //nusquidslv_h

