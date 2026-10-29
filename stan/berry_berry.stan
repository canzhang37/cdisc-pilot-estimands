// AE signal mixture model (Berry and Berry 2004)
// log OR of PT j is 0 with prob pi[soc[j]]
// point mass summed out of the likelihood
data {
  int<lower=1> J;
  int<lower=1> B;
  array[J] int<lower=1, upper=B> soc;
  array[J] int<lower=0> xc;
  array[J] int<lower=0> xt;
  int<lower=1> Nc;
  int<lower=1> Nt;
}
parameters {
  real mu_gamma0;
  real<lower=0> tau_gamma;
  real<lower=0> sigma_gamma;
  vector[B] z_mg;
  vector[J] z_gamma;
  real mu_theta0;
  real<lower=0> tau_theta;
  real<lower=0> sigma_theta;
  vector[B] z_mt;
  vector[J] z_theta;
  real mu_pi;
  real<lower=0> sigma_pi;
  vector[B] z_pi;
}
transformed parameters {
  vector[B] mu_gamma = mu_gamma0 + tau_gamma * z_mg;
  vector[J] gamma = mu_gamma[soc] + sigma_gamma * z_gamma;
  vector[B] mu_theta = mu_theta0 + tau_theta * z_mt;
  vector[J] theta = mu_theta[soc] + sigma_theta * z_theta;
  vector[B] lpi = mu_pi + sigma_pi * z_pi;  // logit of the point-mass probability
}
model {
  mu_gamma0 ~ normal(-3, 2);
  tau_gamma ~ normal(0, 1);
  sigma_gamma ~ normal(0, 1);
  mu_theta0 ~ normal(0, 1);
  tau_theta ~ normal(0, 0.5);
  sigma_theta ~ normal(0, 1);
  mu_pi ~ normal(1.4, 1);
  sigma_pi ~ normal(0, 1);
  z_mg ~ std_normal();
  z_gamma ~ std_normal();
  z_mt ~ std_normal();
  z_theta ~ std_normal();
  z_pi ~ std_normal();
  xc ~ binomial_logit(Nc, gamma);
  for (j in 1:J) {
    target += log_mix(inv_logit(lpi[soc[j]]),
                      binomial_logit_lpmf(xt[j] | Nt, gamma[j]),
                      binomial_logit_lpmf(xt[j] | Nt, gamma[j] + theta[j]));
  }
}
generated quantities {
  vector[J] p_raised;  // P(in the slab and theta > 0), averaged over draws
  vector[J] effect;    // log odds ratio with the indicator drawn
  for (j in 1:J) {
    real l0 = log_inv_logit(lpi[soc[j]]) + binomial_logit_lpmf(xt[j] | Nt, gamma[j]);
    real l1 = log1m_inv_logit(lpi[soc[j]]) + binomial_logit_lpmf(xt[j] | Nt, gamma[j] + theta[j]);
    real w1 = exp(l1 - log_sum_exp(l0, l1));
    p_raised[j] = w1 * (theta[j] > 0);
    effect[j] = bernoulli_rng(w1) * theta[j];
  }
}
