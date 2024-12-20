// Q-learning model with affect and reward both affective outcome valuation
data {
  int<lower=10> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=2> n_f; //number of fractals
  array[n_s,n_t] int fA; // 2d array - one array for each subject, containing a sub-arry with the identity of
                          // fractal A on each of that subject's trials (fractals are identified with numbers 1:n_f)
  array[n_s,n_t] int fB; // diddo for fractal B identities
  array[n_s] vector[n_t] out_a; // 1d array of vectors - one for each subject - that each contains the outcomes 
                          // of fractal A and B for each of the subject's trial
  array[n_s] vector[n_t] out_b; // diddo for outcome of fractal B
  array[n_s,n_t] int choice; // 2d array - one integer array for each subject - containing subject's choices - coded
                             // 1 for fractal A and 2 for fractal B. 
  array[n_s] vector[n_t] feed_rate_0; // post-feedback valence ratings (with 0 - the subject-level mean - substituted for missing ratings)
}

parameters {
  real alpha_mu; // mean of group-level learning rate distribution
  real<lower=0> alpha_sigma; // sd of group-level learning rate distribution
  
  real beta_mu; // mean of group-level inverse temperature distribution
  real<lower=0> beta_sigma; // sd of group-level inverse temperature distribution

  real aff_sens_mu; // mean of group-level affective sensitivity distribution
  real<lower=0> aff_sens_sigma; // sd
  
  real phi_mu; //Mean of group-level autocorrelation weight distribution
  real<lower=0> phi_sigma; //sd of group-level distribution
  
  real tau_mu; //mean of group-level learning rate for C - the autocorrelation value
  real<lower=0> tau_sigma; //sd
  
  real forget_mu; //mean decay rate for the forgetting process
  real<lower=0> forget_sigma;
  
  vector[n_s] alpha_z; // z-score of subject-level alpha value (relative to group distribution)
  vector[n_s] beta_z; // z-score of subject-level beta value (relative to group distribution)
  vector[n_s] aff_sens_z; // z-score of subject-level affect sensitivity
  vector[n_s] phi_z; //z-score of subject-level phi
  vector[n_s] tau_z; //z-score of subject-level tau
  vector[n_s] forget_z;
}

transformed parameters {
  vector[n_s*n_t] log_lik;
  real Q_sd;
  real A_sd;
  real C_sd;
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q; // The Q values. These are a 2d array containing vectors: the first dimension is
                                  // the subject, the second dimension is the trial, and for each trial there's a
                                  // vector with a Q value for each fractal
    array[n_s,n_t] vector[n_f] A; // Expected affective valence, with the same structure as Q
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value, of the same structure as Q
    
    //differences between Q, A, and C values of fractal A and fractal B
    matrix[n_s,n_t] Q_diff;
    matrix[n_s,n_t] A_diff;
    matrix[n_s,n_t] C_diff;
    
    array[n_s] vector[n_t] PE_a;  // The RPEs to outcome a. This is a 1d array containing one vector per subject,
                                  // where that vector contains one value per trial
    array[n_s] vector[n_t] PE_b;  // diddo outcome B
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    
    vector[n_s] beta = beta_mu + beta_sigma*beta_z;  // get subject-level beta values using n.c.p.
    vector[n_s] aff_sens = aff_sens_mu + aff_sens_sigma*aff_sens_z; // ditto aff_sens
    vector[n_s] phi = phi_mu + phi_sigma*phi_z; //ditto phi
    vector[n_s] alpha = inv_logit(alpha_mu + alpha_sigma*alpha_z); //ditto alpha. Also, use inv_logit to get alpha between 0 and 1
    vector[n_s] tau = inv_logit(tau_mu + tau_sigma*tau_z); //likewise for tau
    vector[n_s] forget = inv_logit(forget_mu + forget_sigma*forget_z); //ditto forget
  
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          Q[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          A[s,t] = rep_vector(0,n_f); //ditto A values
          C[s,t] = rep_vector(0,n_f); // ditto C values
        }
        
        // get the likelihood of the choice made on this trial
        log_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | [beta[s]*Q[s,t,fA[s,t]] + aff_sens[s]*A[s,t,fA[s,t]] + phi[s]*C[s,t,fA[s,t]],
                                                                     beta[s]*Q[s,t,fB[s,t]] + aff_sens[s]*A[s,t,fB[s,t]] + phi[s]*C[s,t,fB[s,t]]]'); 
        
        //get differences on this trial
        Q_diff[s,t] = Q[s,t,fA[s,t]] - Q[s,t,fB[s,t]];
        A_diff[s,t] = A[s,t,fA[s,t]] - A[s,t,fB[s,t]];
        C_diff[s,t] = C[s,t,fA[s,t]] - C[s,t,fB[s,t]];
                                                                     
        PE_a[s,t] = out_a[s,t] - Q[s,t,fA[s,t]]; // get the PE to fractal A on this trial
        PE_b[s,t] = out_b[s,t] - Q[s,t,fB[s,t]]; // get the PE to fractal B on this trial

        //unlesss this is the subject's last trial, set the Q and C values for the next trial
        if(t < n_t){
          //Decay Q and A values toward 0 by default
          Q[s,t+1] = (1-forget[s])*Q[s,t]; 
          A[s,t+1] = (1-forget[s])*A[s,t];
          
          Q[s,t+1,fA[s,t]] = Q[s,t,fA[s,t]] + alpha[s]*PE_a[s,t]; // update Q value of fractal A
          Q[s,t+1,fB[s,t]] = Q[s,t,fB[s,t]] + alpha[s]*PE_b[s,t]; // update Q value of fractal B
          
          if(choice[s,t] == 1){
            // if fA was chosen, then update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
            //additionally, update the A value for fA only
            A[s,t+1,fA[s,t]] = A[s,t,fA[s,t]] + alpha[s]*(feed_rate_0[s,t] - A[s,t,fA[s,t]]);
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen
            choice_a = 0;
            choice_b = 1;
            A[s,t+1,fB[s,t]] = A[s,t,fB[s,t]] + alpha[s]*(feed_rate_0[s,t] - A[s,t,fB[s,t]]);
          }
          C[s,t+1] = C[s,t]; // Iniitalize the next trial's C values to be the same as the current trial's, as with Q values
          C[s,t+1,fA[s,t]] = C[s,t,fA[s,t]] + tau[s]*(choice_a - C[s,t,fA[s,t]]); // update C value of fA with 1 if fA was chosen and 0 if it wasn't.
          C[s,t+1,fB[s,t]] = C[s,t,fB[s,t]] + tau[s]*(choice_b - C[s,t,fB[s,t]]); // diddo fB
        }
      }
    }
    Q_sd = sd(Q_diff);
    A_sd = sd(A_diff);
    C_sd = sd(C_diff);
  }//anonymous_scope_end
}

model {
  // add the hyperpriors on the group-level parameters to the target density
  alpha_mu ~ normal(-.05,1.7); //This yields a nearly uniform 0-1 prior when logit transformed
  alpha_sigma ~ normal(0,4); //allow for large inter-subject heterogeneity without enforcing marginal subject-level priors that are excessively horeshoe-shaped
  tau_mu ~ normal(-.05,1.7); //setting tau priors same as alpha
  tau_sigma ~ normal(0,4); // ""
  beta_mu ~ normal(1,5); // Weakly informative prior
  beta_sigma ~ normal(0,5); //Setting the prior on the SD with similar logic
  aff_sens_mu ~ normal(0,4); // Weakly informative prior
  aff_sens_sigma ~ normal(0,6);
  phi_mu ~ normal(2,10); // Setting a reguarlizing prior. Matching the phi prior to the beta prior, assuming that uatocorrelation and reward are apt to have 
                         // equivalent effects on choice. But differences in phi between fA and fB will be about .5 at asymptote (based on data in which
                         // subjects tend to repeat choices 75% of the time), whereas it's more like 1.2 for beta, so the prior for phi will be set a bit wider.
  phi_sigma ~ normal(0,10);
  forget_mu ~ normal(-.05,1.7); //setting forget priors the same as alpha
  forget_sigma ~ normal(0,4);
  
  // add the "priors" on the subject-level parameters derived form the group-level parameters to the target density
  alpha_z ~ normal(0,1);
  beta_z ~ normal(0,1);
  phi_z ~ normal(0,1);
  tau_z ~ normal(0,1);
  forget_z ~ normal(0,1);
  aff_sens_z ~ normal(0,1);
  
  // Add the joint likelihood to the target density, mapping beta-adjusted Q values to choice probabilities
  // using a softmax function. Unfortunately, categorical_logit does not support vectorization. 
  target += log_lik;
}

generated quantities{
  real scd_beta_mu = Q_sd*beta_mu;
  real scd_aff_sens_mu = A_sd*aff_sens_mu;
  real scd_phi_mu = C_sd*phi_mu;
}
