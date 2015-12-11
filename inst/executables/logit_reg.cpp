// linear regression
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_MATRIX(x);
DATA_MATRIX(y);
DATA_VECTOR(w);
DATA_VECTOR(offset);

// parameters:
PARAMETER(a)
PARAMETER_VECTOR(b); // intercept and slopes
//PARAMETER_VECTOR(log_b); //log parameters
// we fit sigma on a log scale to keep it > 0

// procedures: (transformed parameters)
int n = w.size(); // get number of data points to loop over
Type ll = 0.0; // initialize negative log likelihood

//int nb = b.size();
vector<Type> log_b = exp(b);
vector<Type> lp = a + x*log_b + offset;
vector<Type> p = invlogit(lp); //Type(1.0)/(Type(1.0)+exp(-lp));

for(int i = 0; i < n; i++){ // C++ starts loops at 0!
  // get negative log likelihood (last argument is log = TRUE)
  ll -= dbinom(y(i,0), y(i,1), p[i], true)*w[i];
}

return ll;
}