#include <Rcpp.h>
using namespace Rcpp;

//' Discretization of a density
//'
//' Discretizates a density function using Simposon rule
//' @param density the density function
//' @param parameters parameters of the density function
//' @param truncation_point max value for the discretization
//' @param max_probability maximum probability of a discrete point
//' @param max_deltat maximum size of interval
//' @return List with values and weights
// [[Rcpp::export]]
List discretizate_density(Function density, NumericVector parameters, double truncation_point, double max_probability, double max_deltat) {
  double deltat{};
  double t{0.00001};

  std::vector<double> prob;
  std::vector<double> values;
  double prob_aux{};


  NumericVector den_aux{};
  NumericVector den_aux2{};
  NumericVector den_aux3{};

  while (t < truncation_point) {
    den_aux = density(t, parameters);
    if (den_aux[0] < max_probability / max_deltat) {
      deltat = max_deltat;
    }
    else {
      deltat = max_probability / den_aux[0];
    }
    den_aux2 = density(t + deltat / 2, parameters);
    den_aux3 = density(t + deltat, parameters);
    prob_aux = deltat / 6 * (den_aux[0] + 4 * den_aux2[0] + den_aux3[0]);
    while (prob_aux > max_probability) {
      deltat = deltat * 0.9;
      den_aux2 = density(t + deltat / 2, parameters);
      den_aux3 = density(t + deltat, parameters);
      prob_aux = deltat / 6 * (den_aux[0]  + 4 * den_aux2[0] + den_aux3[0]);
    }
    if (prob_aux > 0) {
      values.push_back( (t * den_aux[0]  + 4 * (t + deltat / 2) * den_aux2[0] + (t + deltat) * den_aux3[0]) / (den_aux[0]  + 4 * den_aux2[0] + den_aux3[0]));
      prob.push_back(prob_aux);
    }
    t = t + deltat;
  }
  double factor{};
  factor = accumulate(prob.begin(),prob.end(),0.0);
  transform(prob.begin(), prob.end(), prob.begin(), [factor](double &c){ return c/factor; });
  return List::create(
    _["val"] = values,
    _["prob"] = prob
  );
}
