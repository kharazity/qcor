#include "observable_provider.hpp"
#include "xacc_plugin.hpp"
#include "qcor_observable.hpp"

#include "ObservableTransform.hpp"
#include "xacc.hpp"
#include "xacc_quantum_gate_api.hpp"
#include "xacc_service.hpp"

namespace qcor{
    using Observable = xacc::Observable;
    using PauliOperator = xacc::quantum::PauliOperator;
    using FermionOperator = xacc::quantum::FermionOperator;

    //het map options:
    /*
        x-dimension :int
        y-dimension :int
        tunneling :double
        coulomb :double
        chemical-potential=0.
        magnetic-field=0.
        periodic=True
        spinless=False
        particle-hole-symmetry=False
    */

    class HubbardProvider : public ObservableProvider{
        public:
            std::shared_ptr<Observable> createOperator(HeterogeneousMap &options){
                return hubbardHelper(options);
            }
        
        std::shared_ptr<Observable> hubbardHelper(HeterogeneousMap &args){
            int x_dimension;
            int y_dimension;
            double tunneling = 1.;
            double coulomb = 1.;
            double chemical_potential = 0;
            double magnetic_field = 0;
            bool periodic = true;
            bool spinless = false;
            bool particle_hole_symmetry = false;
            //parsing hetmap for arguments
            if (args.keyExists<int>("x-dimension")){
                x_dimension = args.get<int>("x-dimension");
            }
            else{
                xacc::error("must pass x-dimension argument i.e. in Python, {\"x-dimension\": int}");
            }
            if (args.keyExists<int>("y-dimension")){
                y_dimension = args.get<int>("y-dimension");
            }
            else{
                xacc::error("must pass y-dimension argument i.e. in Python,{\"y-dimension\": int}");
            }
            if (args.keyExists<double>("tunneling")){
                tunneling = args.get<double>("tunneling");
            }
            if (args.keyExists<double>("coulomb")){
                coulomb = args.get<double>("coulomb");
            }
            if (args.keyExists<double>("chemical-potential")){
                chemical_potential = args.get<double>("chemical-potential");
            }
            if (args.keyExists<double>("magnetic-field")){
                magnetic_field = args.get<double>("magnetic-field");
            }
            if (args.keyExists<bool>("periodic")){
                periodic = args.get<bool>("periodic");
            }
            if (args.keyExists<bool>("spinless")){
                spinless = args.get<bool>("spinless");
            }
            if(args.keyExists<bool>("particle-hole-symmetry")){
                particle_hole_symmetry = args.get<bool>("particle-hole-symmetry");
            }
            if(spinless){
                return spinlessHubbard(x_dimension, y_dimension, tunneling,
                                        coulomb, chemical_potential, magnetic_field,
                                        periodic, particle_hole_symmetry);
            }
            else{
                return spinfullHubbard(x_dimension, y_dimension, tunneling,
                                        coulomb, chemical_potential, magnetic_field,
                                        periodic, particle_hole_symmetry);
            }
        }

        std::shared_ptr<Observable> spinlessHubbard(int x_dimension, int y_dimension, double tunneling,
                            double coulomb, double chemical_potential, double magnetic_field,
                            bool periodic, bool particle_hole_symmetry){
            int n_sites = x_dimension*y_dimension;
            auto hubbard_model = FermionOperator();
            for(auto &site : range(n_sites)){
                int right_neighbor = _right_neighbor(site, x_dimension, y_dimension, periodic);
                int bottom_neighbor = _bottom_neighbor(site, x_dimension, y_dimension, periodic);

                if((x_dimension == 2) && periodic && (site % 2 == 1)){
                    right_neighbor = -1;
                }
                if((y_dimension == 2) && periodic && (site >= x_dimension)){
                    bottom_neighbor = -1;
                }
                if(right_neighbor != -1){
                    hubbard_model += _hopping_term(site, right_neighbor, -tunneling);
                    hubbard_model += _coulomb_interaction_term(n_sites, site,
                                                                right_neighbor, coulomb,
                                                                particle_hole_symmetry);
                }
                if(bottom_neighbor != -1){
                    hubbard_model += _hopping_term(site, bottom_neighbor, -tunneling);
                    hubbard_model += _coulomb_interaction_term(n_sites, site, bottom_neighbor,
                                                                coulomb, particle_hole_symmetry);
                }

                hubbard_model += -chemical_potential*adag(site)*a(site);
            }//site in sites
            std::shared_ptr<Observable> obs = std::make_shared<FermionOperator>(hubbard_model);
            return obs;
        } //spinlessHubbard

        std::shared_ptr<Observable> spinfullHubbard(int x_dimension, int y_dimension, double tunneling,
                                    double coulomb, double chemical_potential, double magnetic_field,
                                    bool periodic, bool particle_hole_symmetry){
            int n_sites = x_dimension*y_dimension;
            int n_spin_orbitals = 2*n_sites; //two spins per site
            auto hubbard_model = FermionOperator();

            for(auto &site :range(n_sites)){
                int right_neighbor = _right_neighbor(site, x_dimension, y_dimension, periodic);
                int bottom_neighbor = _bottom_neighbor(site, x_dimension, y_dimension, periodic);

                //avoid double counting when x_dim or y_dim = 2, and periodic
                if((x_dimension == 2) && periodic && (site%2 == 1)){
                    right_neighbor = -1;
                }
                if((y_dimension == 2) && periodic && (site >= x_dimension)){
                    bottom_neighbor = -1;
                }
                if(right_neighbor != -1){
                    hubbard_model += _hopping_term(up_index(site), up_index(right_neighbor),
                                                    -tunneling);
                    hubbard_model += _hopping_term(down_index(site), down_index(right_neighbor),
                                                    -tunneling);
                }
                if(bottom_neighbor != -1){
                    hubbard_model += _hopping_term(up_index(site), up_index(bottom_neighbor),
                                                    -tunneling);
                    hubbard_model += _hopping_term(down_index(site), down_index(bottom_neighbor),
                                                    -tunneling);
                }
                hubbard_model += _coulomb_interaction_term(n_spin_orbitals, up_index(site),
                                                            down_index(site), coulomb, 
                                                            particle_hole_symmetry);
                hubbard_model += (-chemical_potential - magnetic_field)*adag(up_index(site))*a(up_index(site));
                hubbard_model += (-chemical_potential + magnetic_field)*adag(down_index(site))*a(down_index(site));
            } //for (site in sites)
            std::shared_ptr<Observable> obs = std::make_shared<FermionOperator>(hubbard_model);
            return obs;
        }//spinfull Hubbard

        int up_index(int index){
            return 2*index;
        }
        int down_index(int index){
            return 2*index + 1;
        }
        FermionOperator _hopping_term(int i, int j, double coeff){
            auto hopping_term = coeff*adag(j)*a(i);
            hopping_term += coeff*adag(i)*a(j);
            return hopping_term;
        }
        FermionOperator _coulomb_interaction_term(int n_site, int i, int j, double coeff, 
                                        bool particle_hole_symmetry){
            auto number_op_i = adag(i)*a(i);
            auto number_op_j = adag(j)*a(j);
            if(particle_hole_symmetry){
                number_op_i -= FermionOperator(0.5);
                number_op_j -= FermionOperator(0.5);
            }
            return coeff*number_op_i*number_op_j;

        }
        int _right_neighbor(int site, int x_dimension, int y_dimension, bool periodic){
            if(x_dimension == 1){
                return -1;
            }
            if((site+1)%x_dimension == 0){
                if(periodic){
                    return site + 1 - x_dimension;
                }
                else{
                    return -1;
                }
            }
            return site + 1;
        }

        int _bottom_neighbor(int site, int x_dimension, int y_dimension, bool periodic){
            if(y_dimension == 1){
                return -1;
            }
            if((site + x_dimension + 1) > x_dimension*y_dimension){
                if(periodic){
                    return site + x_dimension - x_dimension*y_dimension;
                }
                else{
                    return -1;
                }
            }
            return site + x_dimension;
        }
    public:
        const std::string name() const override { return "hubbard"; }
        const std::string description() const override { return ""; }
    }; //HubbardProvider

}//namespace QCOR
REGISTER_PLUGIN(qcor::HubbardProvider, qcor::ObservableProvider)