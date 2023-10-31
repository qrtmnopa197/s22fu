// Q-learning model with choice autocorrelation and side bias.

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
}

parameters {
  real alpha_mu; // mean of group-level learning rate distribution
  real<lower=0> alpha_sigma; // sd of group-level learning rate distribution
  
  real beta_mu; // mean of group-level inverse temperature distribution
  real<lower=0> beta_sigma; // sd of group-level inverse temperature distribution
  
  real phi_mu; //Mean of group-level autocorrelation weight distribution
  real<lower=0> phi_sigma; //sd of group-level distribution
  
  real tau_mu; //mean of group-level learning rate for C - the autocorrelation value
  real<lower=0> tau_sigma; //sd
  
  //left-side bias
  real lsbias_mu;
  real<lower=0> lsbias_sigma;
  
  vector[n_s] alpha_z; // z-score of subject-level alpha value (relative to group distribution)
  vector[n_s] beta_z; // z-score of subject-level beta value (relative to group distribution)
  vector[n_s] phi_z; //z-score of subject-level phi
  vector[n_s] tau_z; //z-score of subject-level tau
  vector[n_s] lsbias_z; //z-score of subject-level side bias
}

transformed parameters {
  vector[n_s*n_t] log_lik;
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q; // The Q values. These are a 2d array containing vectors: the first dimension is
                                  // the subject, the second dimension is the trial, and for each trial there's a
                                  // vector with a Q value for each fractal
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value, of the same structure as Q
    
    array[n_s] vector[n_t] PE_a;  // The RPEs to outcome a. This is a 1d array containing one vector per subject,
                                  // where that vector contains one value per trial
    array[n_s] vector[n_t] PE_b;  // diddo outcome B
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    
    vector[n_s] beta = beta_mu + beta_sigma*beta_z;  // get subject-level beta values using n.c.p.
    vector[n_s] phi = phi_mu + phi_sigma*phi_z; //diddo phi
    vector[n_s] alpha = inv_logit(alpha_mu + alpha_sigma*alpha_z); //diddo alpha. Also, use inv_logit to get alpha between 0 and 1
    vector[n_s] tau = inv_logit(tau_mu + tau_sigma*tau_z); //likewise for tau
    vector[n_s] lsbias = lsbias_mu + lsbias_sigma*lsbias_z; //ditto lsbias

  
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          Q[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          C[s,t] = rep_vector(0,n_f); // diddo C values
        }
        
        log_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | beta[s]*Q[s,t,{fA[s,t],fB[s,t]}] + phi[s]*C[s,t,{fA[s,t],fB[s,t]}] + [lsbias[s],0]'); //get the likelihood of the choice made on this trial.                                                                                                                                                                        //here prevents vectorization, but you have to calculate a vector
                                                                                                                                                                        
        PE_a[s,t] = out_a[s,t] - Q[s,t,fA[s,t]]; // get the PE to fractal A on this trial
        PE_b[s,t] = out_b[s,t] - Q[s,t,fB[s,t]]; // get the PE to fractal B on this trial

        //unlesss this is the subject's last trial, set the Q and C values for the next trial
        if(t < n_t){
          Q[s,t+1] = Q[s,t]; // Initiliaze the next trial's Q values to be the same as the current trial's; all
                             // Q values will indeed be the same, except for those for which there was a PE on
                             // this trial...
          Q[s,t+1,fA[s,t]] = Q[s,t,fA[s,t]] + alpha[s]*PE_a[s,t]; // update Q value of fractal A
          Q[s,t+1,fB[s,t]] = Q[s,t,fB[s,t]] + alpha[s]*PE_b[s,t]; // update Q value of fractal B
          
          if(choice[s,t] == 1){
            // if fA was chosen, then update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen
            choice_a = 0;
            choice_b = 1;
          }
          C[s,t+1] = C[s,t]; // Iniitalize the next trial's C values to be the same as the current trial's, as with Q values
          C[s,t+1,fA[s,t]] = C[s,t,fA[s,t]] + tau[s]*(choice_a - C[s,t,fA[s,t]]); // update C value of fA with 1 if fA was chosen and 0 if it wasn't.
          C[s,t+1,fB[s,t]] = C[s,t,fB[s,t]] + tau[s]*(choice_b - C[s,t,fB[s,t]]); // diddo fB
        }
      }
    }
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
  phi_mu ~ normal(2,10); // Setting a reguarlizing prior. Matching the phi prior to the beta prior, assuming that uatocorrelation and reward are apt to have 
                         // equivalent effects on choice. But differences in phi between fA and fB will be about .5 at asymptote (based on data in which
                         // subjects tend to repeat choices 75% of the time), whereas it's more like 1.2 for beta, so the prior for phi will be set a bit wider.
  phi_sigma ~ normal(0,10);
  lsbias_mu ~ normal(0,4);
  lsbias_sigma ~ normal(0,8);

  
  // add the "priors" on the subject-level parameters derived form the group-level parameters to the target density
  alpha_z ~ normal(0,1);
  beta_z ~ normal(0,1);
  phi_z ~ normal(0,1);
  tau_z ~ normal(0,1);
  lsbias_z ~ normal(0,1);
  
  // Add the joint likelihood to the target density, mapping beta-adjusted Q values to choice probabilities
  // using a softmax function. Unfortunately, categorical_logit does not support vectorization. 
  target += log_lik;
}
