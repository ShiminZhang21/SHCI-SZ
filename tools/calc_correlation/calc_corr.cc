#include <hps/src/hps.h>
#include <src/det/det.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <algorithm>  // for std::find
#include <cstdio>
#include <cmath>

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

bool contains_orb(const std::vector<unsigned>& orbs, unsigned j)
  // condition to tell is an orbital is in the occupied list
  // usage contains_orb(orbital_list, i_orbital)
{
    return std::find(orbs.begin(), orbs.end(), j) != orbs.end();
}

int main(int argc, char *argv[]) {
  //usage:: exe <wf_fileanme name> <number of orbitals>
    if (argc < 2) {
        printf("Usage: exe <wf_fileanme name> <number of orbitals>\n");
        return 1;
    }
  //read wavefunction 
    std::ifstream serialized_wf(argv[1], std::ios::binary); 
    Wavefunction wf = hps::from_stream<Wavefunction>(serialized_wf);
    unsigned n_orbs = 0; // number of orbitals(just give a randum initializaed orbital)
    n_orbs = std::stoul(argv[2]);
    //printf("wavefunction coefficent size(number of det): %u\n",wf.coefs[0].size());
    //number of states
    int n_states = wf.coefs.size(); // number of states
    int n_det = wf.dets.size(); // number of slater determinants
    printf("Active space size: %u\n", n_orbs);
    printf("n_states: %u\n", n_states);
    printf("n_det: %u\n", n_det);

    // Compute the spin-resolved orbital occupation:
    std::vector<double> Static_correlations(n_states, 0.0);//initialize
    std::vector<double> Dynamic_correlations(n_states, 0.0);//initialize
    for (unsigned i_state = 0; i_state < n_states; i_state++) {
      printf("===============State %u==================\n",i_state);
      // Compute orb occupations.
      std::vector<double> orb_occupations_up(n_orbs, 0.0);//initialize
      std::vector<double> orb_occupations_dn(n_orbs, 0.0);//initialize
      std::vector<double> orb_occupations(n_orbs, 0.0);//initialize
      for (unsigned j = 0; j < n_orbs; j++) { //j for orbital
        for (size_t i = 0; i < n_det; i++) { //i for determinant
          const Det &det = wf.dets[i];
          std::vector<unsigned> up_orbs(det.up.get_occupied_orbs());
          std::vector<unsigned> dn_orbs(det.dn.get_occupied_orbs());
          const double coef = wf.coefs[i_state][i]; //coefficient of ith determinant
          
          if (contains_orb(up_orbs, j)) { //if determinant have orbital j, then orbital occupation + coefficient**2
            orb_occupations_up[j] += coef * coef;
          }
          if (contains_orb(dn_orbs, j)) {
            orb_occupations_dn[j] += coef * coef;
          }
        }
        orb_occupations[j] = orb_occupations_up[j] + orb_occupations_dn[j];
      }
      //Compute correlations:
      for (unsigned j = 0; j < n_orbs; j++) {
        Static_correlations[i_state] += 
              0.5 * orb_occupations_up[j]*(1-orb_occupations_up[j])
            + 0.5 * orb_occupations_dn[j]*(1-orb_occupations_dn[j]);
        Dynamic_correlations[i_state] += 
              0.25 * std::sqrt(orb_occupations_up[j]*(1-orb_occupations_up[j]))
              -0.5 * orb_occupations_up[j]*(1-orb_occupations_up[j])
              +0.25 * std::sqrt(orb_occupations_dn[j]*(1-orb_occupations_dn[j]))
              -0.5 * orb_occupations_dn[j]*(1-orb_occupations_dn[j]);
      }
      //print the orb occupations and correlation metric
      printf("----------------------------------------\n");
      printf("%-10s%12s%16s%6s%16s%6s%16s\n", "Orbital", "", "Sum c^2(tot)","", "Sum c^2(up)","", "Sum c^2(dn)");
      //for (unsigned j = 0; j < system.n_orbs && j < 50; j++) {
      for (unsigned j = 0; j < n_orbs; j++) { //change: print all orbitals
        printf("%-10u%12s%16.8f%6s%16.8f%6s%16.8f\n", j, "", orb_occupations[j],"", orb_occupations_up[j],"", orb_occupations_dn[j]);
      }
      double sum_orb_occupation = std::accumulate(orb_occupations.begin(), orb_occupations.end(), 0.0);
      double sum_orb_occupation_up = std::accumulate(orb_occupations_up.begin(), orb_occupations_up.end(), 0.0);
      double sum_orb_occupation_dn = std::accumulate(orb_occupations_dn.begin(), orb_occupations_dn.end(), 0.0);
      printf("Sum orbitals c^2 for total, up, dn: %.8f%6s%.8f%6s%.8f\n", sum_orb_occupation,"",sum_orb_occupation_up,"",sum_orb_occupation_dn);
      printf("Static correlation, Dynamic correlation: %.8f%6s%.8f\n", Static_correlations[i_state],"",Dynamic_correlations[i_state]);
    }
    //Summary of correlation: 
    printf("-----------------Summary of correlations-----------------\n");
    printf("%-10s%12s%16s%6s%16s\n", "State", "", "Static correlation","", "Dynamic correlation");
    for (unsigned i_state = 0; i_state < n_states; i_state++) {
    printf("%-10u%12s%16.8f%6s%16.8f\n", i_state, "", Static_correlations[i_state],"", Dynamic_correlations[i_state]);
    }
  };
    
