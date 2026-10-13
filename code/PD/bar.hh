#ifndef BAR_HH
#define BAR_HH

#include <cstddef>
#include <vector>

struct BondData {
  std::vector<std::size_t> i;
  std::vector<std::size_t> j;
  std::vector<double> length;
  std::vector<double> overlap;
  std::vector<double> k;
  std::vector<double> energy_coef;
  std::vector<double> s_elastic;
  std::vector<double> s_break;
  std::vector<bool> can_break;
  std::size_t count;
};

struct Bar {
  double L;
  double dx;
  double horizon;
  double m;
  std::size_t num_points;
  std::size_t num_boundary;
  std::size_t num_total;
  std::vector<double> x;
  std::vector<std::size_t> inside;
  std::vector<std::size_t> boundary;
  double area;
  double volume;
  double wave_speed;
  double dt;
  BondData bond;
  double mean_s_break;
  std::vector<double> strength_pt;
  std::vector<double> toughness_pt;
};

#endif