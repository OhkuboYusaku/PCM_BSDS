data {
  // import data sizes //
    int<lower=0> N;
    int<lower=0> N_sp;
    int<lower=0> D;
   
  // import regressin model data //
    vector[N] y;        // objective variable
    matrix[N, D] X;   // explanatory variable
    matrix[N, N_sp] Z;  // spices id for random effect

}

 parameters {
  vector[D] beta;
  real<lower=0> sigma_y;
  vector[N_sp] b;
}

model {

  target+= normal_lpdf(y| X * beta + Z * b, sigma_y);

  //target+= poisson_log(y| mu);
  //target+= bernoulli_logit(y| mu);
}

generated quantities {
vector[N] log_likelihood;

  for(n in 1:N){log_likelihood[n] = normal_lpdf(y[n]| X[n,] * beta + Z[n,] * b, sigma_y);}
}


