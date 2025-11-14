#include <hps/src/hps.h>
#include <shci/src/det/det.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <algorithm>  // for std::find
#include <cstdio>

constexpr double SQRT2_INV = 0.7071067811865475;

class Wavefunction {
public:
  unsigned n_up = 0;

  unsigned n_dn = 0;

  double energy_hf = 0.0;

  std::vector<double> energy_var;

  bool time_sym = false;

  std::vector<Det> dets;

  std::vector<std::vector<double>> coefs;

  size_t get_n_dets() const { return dets.size(); }                                                                                                          
                                                                                                                                                             
  void unpack_time_sym() {                                                                                                                                   
    const size_t n_dets_old = get_n_dets();                                                                                                                  
    for (size_t i = 0; i < n_dets_old; i++) {                                                                                                                
      const auto& det = dets[i];                                                                                                                             
      if (det.up < det.dn) {                                                                                                                                 
        Det det_rev = det;                                                                                                                                   
        det_rev.reverse_spin();                                                                                                                              
        for (auto& state_coefs: coefs) {                                                                                                                     
          const double coef_new = state_coefs[i] * SQRT2_INV;                                                                                                
          state_coefs[i] = coef_new;                                                                                                                         
          state_coefs.push_back(coef_new);                                                                                                                   
        }                                                                                                                                                    
        dets.push_back(det_rev);                                                                                                                             
      }                                                                                                                                                      
    }                                                                                                                                                        
  }

  template <class B>
  void parse(B& buf) {
    buf >> n_up >> n_dn >> dets >> coefs >> energy_hf >> energy_var >> time_sym;
    if (time_sym) unpack_time_sym();
  }
};

int main(int argc, char *argv[]) {
  //usage:: exe wf_fileanme --threshold 1e-3 to only print coefficient>1e-3
    if (argc < 2) {
        printf("Usage: exe wf_filename\n");
        return 1;
    }
  //determinant print threshold
    bool threshold_provided = false;
    double threshold = 0.0;  // print out threshold
    bool norb_provided = false; //if not norb privded, use the occupied orbital printout
    unsigned norbs = 0; // number of orbitals(just give a randum initializaed orbital)
    std::vector<unsigned> orbs(norbs);//list of orbitals
    std::vector<int> new_orbs(norbs); //change the name of the orbitals

    for (int i = 1; i < argc - 1; ++i) {
      if (std::string(argv[i]) == "--threshold") {
          threshold = std::stod(argv[i + 1]);  // convert the next argument to double
          threshold_provided = true;
          std::cout << "Using threshold: " << threshold << std::endl;
      }
      if (std::string(argv[i])== "--norbs"){ // number of orbitals
          norbs = std::stod(argv[i + 1]);
          norb_provided = true;
          printf("norb given: %u\n", norbs);
          // Fill with values 0 to norbs - 1: {0, 1, ..., norbs-1}
          //generate a list of orbitals:
          orbs.resize(norbs);
          for (unsigned i = 0; i < norbs; ++i) {
              orbs[i] = i;
          }
          
          // Print the list
          printf("Orbital list: ");
          for (unsigned i = 0; i < orbs.size(); ++i) {
              printf("%u ", orbs[i]);
          }
          printf("\n");
          //the new orbital order is according to the out
          new_orbs = {8, 0, 1, 3, 2, 4, 5, 6, 7};
          //for (unsigned orb : orbs) {
          //    unsigned mapped = (orb == 0) ? norbs - 1 : orb - 1;
          //    new_orbs.push_back(mapped);
          //}
          // Print the list
          printf("The orbital is reordered by HF energy from shci ");
          for (unsigned i = 0; i < new_orbs.size(); ++i) {
              printf("%u ", new_orbs[i]);
          }
          printf("\n");
          //printf("the order of orbitals are reversed in bitstring\n");
          //std::reverse(new_orbs.begin(), new_orbs.end()); //reverse the orbital order
           
      }
    }

    std::ifstream serialized_wf(argv[1], std::ios::binary); 
    Wavefunction wf = hps::from_stream<Wavefunction>(serialized_wf);

    //printf("wavefunction coefficent size(number of det): %u\n",wf.coefs[0].size());
    //number of states
    int n_states = wf.coefs.size(); // number of states
    int n_det = wf.dets.size(); // number of slater determinants
    printf("n_states: %u\n", n_states);
    printf("n_det: %u\n", n_det);

    //indices of sorted coefs/dets
    std::vector<size_t> inds(n_det);


    for (int i_state = 0; i_state < n_states; ++i_state) { // print slater determinants of all states
      printf("===============State %u==================\n",i_state);
      //fill with consecutive ints
      std::iota(inds.begin(), inds.end(), 0); // index of determinant from 0 to n_det
      //sort by coef magnitude (first wavefunction)
      std::sort(inds.begin(), inds.end(), [&](const size_t &a, const size_t &b) {
          return std::abs(wf.coefs[i_state][a]) > std::abs(wf.coefs[i_state][b]);
      });
      //print with threshold: 
      size_t N_cutoff = inds.size();  // Number of det choose
      if (threshold_provided == true){ //setup threshold    
        for (size_t i = 0; i < inds.size(); ++i) {
            double coef = std::abs(wf.coefs[i_state][inds[i]]);
            if (coef < threshold) {
                N_cutoff = i;// 2. Find cutoff index: first where coef drops below threshold
                break;
            }
        }
      }

      //print the coefficient of each determinants:
      for(size_t j = 0; j < N_cutoff; ++j) {
          size_t i = inds[j];  // actual index into wf.coefs and wf.dets
          const Det &det = wf.dets[i];
          std::vector<unsigned> up_orbs(det.up.get_occupied_orbs());
          std::vector<unsigned> dn_orbs(det.dn.get_occupied_orbs());

          //print coefs of current determinants
          printf("%16.12e", wf.coefs[i_state][i]);
          //print in bit format if number of orbitals provided
          if (norb_provided){
            printf("|");
            //from left to right ordered from largest to smallest orbitals
            for (unsigned i_orb : new_orbs){
              //printf("iorb=%u \n",i_orb);
              if (std::find(up_orbs.begin(), up_orbs.end(), i_orb) != up_orbs.end()) {
                printf("1");
              }
              else {
                printf("0");
              }
            }
            printf(">|");
            for (unsigned i_orb : new_orbs){
              if (std::find(dn_orbs.begin(), dn_orbs.end(), i_orb) != dn_orbs.end()) {
                printf("1");
              }
              else {
                printf("0");
              }
            }
            printf(">\n");
          }
          //print in occupied band format if number of orbitals provided
          else{
            //print up orbs(in orbital format)
            printf("|");
            for(unsigned orb : up_orbs) {
                printf("%u ", orb);
            }
            //print bit format det(up):
            //det.up.print();

            //print in bit format
            printf(">|");
            //print dn orbs
            //printf("\t");
            for(unsigned orb : dn_orbs) {
                printf("%u ", orb);
            }
            //print bit format det(dn):
            //det.dn.print();
            printf(">\n");
          }

      }
          printf("=============================================\n");
    }
    return 0;
  };
    
