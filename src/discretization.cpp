#include <Rcpp.h>
using namespace Rcpp;

//' Discretization of an univariate density
//'
//' Discretizates an univariate density function using Simposon rule
//' @param density the density function
//' @param parameters parameters of the density function
//' @param ini_point initial value for the discretization
//' @param truncation_point max value for the discretization
//' @param max_probability maximum probability of a discrete point
//' @param max_deltat maximum size of interval
//' @return List with values and weights
// [[Rcpp::export]]
List discretizate_density(Function density, NumericVector parameters, double ini_point, double truncation_point, double max_deltat, double max_probability) {
  double deltat{};
  double t{ini_point};

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
  //double factor{};
  //factor = accumulate(prob.begin(),prob.end(),0.0);
  //transform(prob.begin(), prob.end(), prob.begin(), [factor](double &c){ return c/factor; });
  return List::create(
    _["val"] = values,
    _["prob"] = prob
  );
}




//' Discretization of an bivariate density
//'
//' Discretizates an bivariate density function using Simposon rule
//' @param density the density function
//' @param parameters parameters of the density function
//' @param ini_point1 max value for the discretization - first component
//' @param truncation_point1 max value for the discretization - first component
//' @param max_deltat1 maximum size of interval - first component
//' @param ini_point2 max value for the discretization - second component
//' @param truncation_point2 max value for the discretization - second component
//' @param max_deltat2 maximum size of interval - second component
//' @param max_probability maximum probability of a discrete point
//' @return List with values and weights
// [[Rcpp::export]]
List discretizate_bivdensity(Function density, NumericVector parameters, double ini_point1, double truncation_point1, double max_deltat1, double ini_point2, double truncation_point2, double max_deltat2, double max_probability) {
  double deltat1{};
  double t1{ini_point1};

  double deltat2{};
  double t2{ini_point2};

  std::vector<double> prob;
  std::vector<double> values1;
  std::vector<double> values2;
  double prob_aux{};

  NumericVector den_aux{};
  NumericVector den_aux2{};
  NumericVector den_aux3{};
  NumericVector den_aux4{};
  NumericVector den_aux5{};
  NumericVector den_aux6{};
  NumericVector den_aux7{};
  NumericVector den_aux8{};
  NumericVector den_aux9{};

  while (t2 < truncation_point2) {
    t1 = ini_point1;
    while (t1 < truncation_point1) {
      deltat1 = max_deltat1;
      deltat2 = max_deltat2;

      den_aux = density(t1, t2, parameters);
      den_aux2 = density(t1, t2 + deltat2 / 2, parameters);
      den_aux3 = density(t1, t2 + deltat2, parameters);
      den_aux4 = density(t1 + deltat1 / 2, t2, parameters);
      den_aux5 = density(t1 + deltat1 / 2, t2 + deltat2 / 2, parameters);
      den_aux6 = density(t1 + deltat1 / 2, t2 + deltat2, parameters);
      den_aux7 = density(t1 + deltat1, t2, parameters);
      den_aux8 = density(t1 + deltat1, t2 + deltat2 / 2, parameters);
      den_aux9 = density(t1 + deltat1, t2 + deltat2, parameters);

      prob_aux = deltat1 * deltat2 / 36 * ( den_aux[0] + 4 * den_aux2[0] + den_aux3[0] + 4 * den_aux4[0] + 16 * den_aux5[0] + 4 * den_aux6[0] + den_aux7[0] + 4 * den_aux8[0] + den_aux9[0]);
      while (prob_aux > max_probability) {
        deltat1 = deltat1 * 0.9;
        deltat2 = deltat2 * 0.9;
        den_aux2 = density(t1, t2 + deltat2 / 2, parameters);
        den_aux3 = density(t1, t2 + deltat2, parameters);
        den_aux4 = density(t1 + deltat1 / 2, t2, parameters);
        den_aux5 = density(t1 + deltat1 / 2, t2 + deltat2 / 2, parameters);
        den_aux6 = density(t1 + deltat1 / 2, t2 + deltat2, parameters);
        den_aux7 = density(t1 + deltat1, t2, parameters);
        den_aux8 = density(t1 + deltat1, t2 + deltat2 / 2, parameters);
        den_aux9 = density(t1 + deltat1, t2 + deltat2, parameters);
        prob_aux = deltat1 * deltat2 / 36 * ( den_aux[0] + 4 * den_aux2[0] + den_aux3[0] + 4 * den_aux4[0] + 16 * den_aux5[0] + 4 * den_aux6[0] + den_aux7[0] + 4 * den_aux8[0] + den_aux9[0]);
      }
      if (prob_aux > 0) {
        values1.push_back((t1 * den_aux[0] + 4 * (t1 + deltat1 / 2) * den_aux4[0] + (t1 + deltat1) * den_aux7[0]) / (den_aux[0] + 4 * den_aux4[0] + den_aux7[0]));
        values2.push_back((t2 * den_aux[0] + 4 * (t2 + deltat2 / 2) * den_aux2[0] + (t2 + deltat2) * den_aux3[0]) / (den_aux[0] + 4 * den_aux2[0] + den_aux3[0]));
        prob.push_back(prob_aux);
      }
      t1 = t1 + deltat1;
    }
    t2 = t2 + deltat2;
  }

  return List::create(
    _["val1"] = values1,
    _["val2"] = values2,
    _["prob"] = prob
  );
}
