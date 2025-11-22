// Discrete-time survival model with cloglog link,
// complex survey weights, and BYM2 (structured + unstructured spatial) region effects.
//
// One row per individual * per interval (person-period data).
// Time is included as covariates in X (e.g., time dummies or splines).

data {
  // ----- Basic survival data -----
  int<lower=1> N;                 // number of person-period rows
  int<lower=1> K;                 // number of covariates (includes time effects)
  int<lower=1> T;                 // number of discrete time intervals
  int<lower=1> R;                 // number of regions

  int<lower=0, upper=1> y[N];     // event indicator in this interval
  matrix[N, K] X;                 // covariate matrix (time in here)
  int<lower=1, upper=R> region_id[N];  // region index for each row

  vector<lower=0>[N] w;           // survey weights (will be normalized)

  // ----- Spatial adjacency for structured effects -----
  // Undirected adjacency represented as an edge list.
  int<lower=1> N_edges;           // number of edges in adjacency graph
  int<lower=1, upper=R> node1[N_edges];
  int<lower=1, upper=R> node2[N_edges];

  // Scaling factor for ICAR component (BYM2-style)
  // Compute this in R/Python based on adjacency so that the ICAR field
  // has approximately unit marginal variance when tau = 1.
  real<lower=0> scaling_factor;

  // ----- Prediction setup for region-level survival curves -----
  // A single "typical" covariate trajectory over time, used for all regions.
  matrix[T, K] X_pred;            // X_pred[tt, :] = covariates for interval tt

  int<lower=1, upper=T> horizon;  // e.g., 5 for survival from 0 to 5
}

transformed data {
  // Normalize weights to keep them O(1) for numerical stability
  vector[N] w_norm;
  real w_mean = mean(w);
  for (n in 1:N) {
    w_norm[n] = w[n] / w_mean;
  }
}

parameters {
  // ----- Fixed effects -----
  vector[K] beta;                 // covariate effects (including time)

  // ----- BYM2 spatial components -----
  vector[R] u_struct_raw;         // structured (ICAR) component (before centering)
  vector[R] u_unstruct;           // unstructured iid component

  real<lower=0> sigma_theta;      // overall SD of combined spatial effect
  real<lower=0, upper=1> phi;     // proportion of variance that is structured (ICAR)
}

transformed parameters {
  vector[R] u_struct;             // structured ICAR component, sum-to-zero
  vector[R] theta;                // combined BYM2 region effect

  // Enforce sum-to-zero constraint for ICAR identifiability
  u_struct = u_struct_raw - mean(u_struct_raw);

  {
    real sqrt_phi      = sqrt(phi);
    real sqrt_1mphi    = sqrt(1.0 - phi);

    for (r in 1:R) {
      theta[r] = sigma_theta * (sqrt_phi   * u_struct[r]
                              + sqrt_1mphi * u_unstruct[r]);
    }
  }
}

model {
  // ----- Priors -----
  beta ~ normal(0, 5);

  // BYM2 hyperpriors (you can tune these)
  sigma_theta ~ normal(0, 1);             // overall spatial SD
  phi         ~ beta(0.5, 0.5);           // preference for extremes (0 or 1) is modest

  // Unstructured effects: iid N(0,1)
  u_unstruct ~ normal(0, 1);

  // Structured ICAR prior on u_struct
  // Scaled by 'scaling_factor' so the ICAR field has roughly unit marginal variance.
  {
    for (e in 1:N_edges) {
      int r1 = node1[e];
      int r2 = node2[e];
      real diff = u_struct[r1] - u_struct[r2];
      target += -0.5 * square(diff) / scaling_factor;
    }
    // No explicit variance parameter here: u_struct is unit-scaled,
    // and sigma_theta + phi control the overall variance and mixing.
  }

  // ----- Likelihood with cloglog link and survey weights -----
  for (n in 1:N) {
    real eta    = dot_product(X[n], beta) + theta[region_id[n]];
    real hazard = 1 - exp(-exp(eta));   // cloglog^-1(eta), hazard in (0,1)

    target += w_norm[n] * bernoulli_lpmf(y[n] | hazard);
  }
}

generated quantities {
  // ----- Region-level hazards and survival for a typical covariate path -----
  vector[T] hazard_region[R];    // h_r(t)
  vector[T] S_region[R];         // S_r(t)
  real S0_h_region[R];           // survival from 0 to 'horizon' in each region

  for (r in 1:R) {
    // 1) Interval-specific hazards for region r
    for (tt in 1:T) {
      real eta_pred = dot_product(X_pred[tt], beta) + theta[r];
      hazard_region[r, tt] = 1 - exp(-exp(eta_pred));
    }

    // 2) Survival as product of (1 - hazard) across intervals
    // S_r(0) = 1 by definition.
    S_region[r, 1] = 1.0 - hazard_region[r, 1];
    for (tt in 2:T) {
      S_region[r, tt] = S_region[r, tt - 1] * (1.0 - hazard_region[r, tt]);
    }

    // 3) Survival from 0 to 'horizon' (e.g., 0–5)
    S0_h_region[r] = S_region[r, horizon];
    // This equals product_{tt=1..horizon} (1 - hazard_region[r, tt]).
  }
}

