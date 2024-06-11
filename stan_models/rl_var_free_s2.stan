// Q-learning model with choice autocorrelation and forgetting (independent decay rate).
// Additionally, this includes a linear model of valence ratings predicted by variables from all theories.
// Non-centered parameterization is used for all hierarchical parameters.

data {
  //Choice model data
  int<lower=10> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=2> n_f; //number of fractals
  array[n_s,n_t] int fA; // 2d array - one array for each subject, containing a sub-arry with the identity of
                          // fractal A (i.e., the lefthand fractal) on each of that subject's trials (fractals are identified with indices 1:n_f).
  array[n_s,n_t] int fB; // ditto for fractal B (the righthand fractal) identities
  array[n_s] vector[n_t] out_a; // 1d array of vectors - one for each subject - that each contains the outcomes 
                          // of fractal A and B for each of the subject's trial
  array[n_s] vector[n_t] out_b; // ditto for outcome of fractal B
  array[n_s,n_t] int choice; // 2d array - one integer array for each subject - containing subject's choices - coded
                             // 1 for fractal A and 2 for fractal B. 
 
  array[n_s,n_t] int chosen_frac; //indices of the chosen fractal for each trial
  array[n_s,n_t] int unchosen_frac; //ditto unchosen fractals     
  
  array[n_s] vector[n_t] chosen_out; //the outcomes of the chosen fractals
  array[n_s] vector[n_t] unchosen_out;  //ditto unchosen fractals
}

transformed data{
  int n_hp = 8;
}

parameters {
  //Choice parameters
  real alpha_mu; // mean of group-level learning rate distribution
  real<lower=0> alpha_sigma; // sd of group-level learning rate distribution
  vector[n_s] alpha_z; // z-score of subject-level alpha value (relative to group distribution)

  
  real beta_mu; // mean of group-level inverse temperature distribution
  real<lower=0> beta_sigma; // sd of group-level inverse temperature distribution
  vector[n_s] beta_z; // z-score of subject-level beta value (relative to group distribution)
  
  real phi_mu; //Mean of group-level autocorrelation weight distribution
  real<lower=0> phi_sigma; //sd of group-level distribution
  vector[n_s] phi_z; //z-score of subject-level phi
  
  real tau_mu; //mean of group-level learning rate for C - the autocorrelation value
  real<lower=0> tau_sigma; //sd
  vector[n_s] tau_z; //z-score of subject-level tau
  
  //the rate for the decay toward 0 of fractals not presented
  real forget_mu; //mean 
  real<lower=0> forget_sigma; //sd 
  vector[n_s] forget_z; //z-score
  
  real vt_w_mu; 
  real<lower=0> vt_w_sigma;
  vector[n_s] vt_w_z; 
  
  real ev_w_mu; 
  real<lower=0> ev_w_sigma;
  vector[n_s] ev_w_z; 
  
  real ru_w_mu; 
  real<lower=0> ru_w_sigma;
  vector[n_s] ru_w_z; 
}


transformed parameters {
  vector[n_s*n_t] choice_lik; //log likelihood of each choice made, according to the model

  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q; // The Q values. These are a 2d array containing vectors: the first dimension is
                                  // the subject, the second dimension is the trial, and for each trial there's a
                                  // vector with a Q value for each fractal
    array[n_s,n_t] vector[n_f] X; // Expected affective valence, with the same structure as Q
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value, of the same structure as Q
    
    array[n_s] vector[n_t] PE_a;  // The RPEs to outcome a. This is a 1d array containing one vector per subject,
                                  // where that vector contains one value per trial
    array[n_s] vector[n_t] PE_b;  // ditto outcome B
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    

    
    vector[n_s] beta = beta_mu + beta_sigma*beta_z;  // get subject-level beta values
    vector[n_s] phi = phi_mu + phi_sigma*phi_z; //ditto phi
    vector[n_s] alpha = inv_logit(alpha_mu + alpha_sigma*alpha_z); //ditto alpha. Also, use inv_logit to constrain alpha to be between 0 and 1
    vector[n_s] tau = inv_logit(tau_mu + tau_sigma*tau_z); //likewise for tau
    vector[n_s] forget = inv_logit(forget_mu + forget_sigma*forget_z); //likewise for forgetting rate
  
    vector[n_s] vt_w = vt_w_mu + vt_w_sigma*vt_w_z;
    vector[n_s] ev_w = ev_w_mu + ev_w_sigma*ev_w_z;
    vector[n_s] ru_w = ru_w_mu + ru_w_sigma*ru_w_z;
    
    real vt;
    real ev;
    real ru;
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          Q[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          X[s,t] = rep_vector(0,n_f); //ditto A values
          C[s,t] = rep_vector(0,n_f); // ditto C values
        }

        choice_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | beta[s]*Q[s,t,{fA[s,t],fB[s,t]}] + X[s,t,{fA[s,t],fB[s,t]}] + phi[s]*C[s,t,{fA[s,t],fB[s,t]}]); //get the likelihood of the choice made on this trial.
        
        PE_a[s,t] = out_a[s,t] - Q[s,t,fA[s,t]]; // get the PE for outcome a on this trial
        PE_b[s,t] = out_b[s,t] - Q[s,t,fB[s,t]]; // get the PE for outcome b on this trial
        
        vt = chosen_out[s,t]*exp(choice_lik[(s-1)*n_t+t]) + unchosen_out[s,t]*(1-exp(choice_lik[(s-1)*n_t+t]));
        ev = Q[s,t,chosen_frac[s,t]]; //assign Q_chosen,feed.
        ru = unchosen_out[s,t];
        
        //unless this is the subject's last trial, set the Q, A, and C values for the next trial
        if(t < n_t){
          Q[s,t+1] = (1-forget[s])*Q[s,t]; // Decay Q values of fractals not presented toward 0
          X[s,t+1] = (1-forget[s])*X[s,t]; //ditto A
          
          Q[s,t+1,fA[s,t]] = Q[s,t,fA[s,t]] + alpha[s]*PE_a[s,t]; // update Q value of fractal A
          Q[s,t+1,fB[s,t]] = Q[s,t,fB[s,t]] + alpha[s]*PE_b[s,t]; // update Q value of fractal B
          
          if(choice[s,t] == 1){
            // if fA was chosen, then update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
            //additionally, update the A value for fA only
            X[s,t+1,fA[s,t]] = X[s,t,fA[s,t]] + alpha[s]*(vt_w[s]*vt + ev_w[s]*ev + ru_w[s]*ru - X[s,t,fA[s,t]]);
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen
            choice_a = 0;
            choice_b = 1;
            X[s,t+1,fB[s,t]] = X[s,t,fB[s,t]] + alpha[s]*(vt_w[s]*vt + ev_w[s]*ev + ru_w[s]*ru - X[s,t,fB[s,t]]);
          }
          C[s,t+1] = C[s,t]; // Iniitalize the next trial's C values to be the same as the current trial's
          C[s,t+1,fA[s,t]] = C[s,t,fA[s,t]] + tau[s]*(choice_a - C[s,t,fA[s,t]]); // update C value of fA with 1 if fA was chosen and 0 if it wasn't.
          C[s,t+1,fB[s,t]] = C[s,t,fB[s,t]] + tau[s]*(choice_b - C[s,t,fB[s,t]]); // ditto fB
        }
      }
    }
  }//anonymous_scope_end
}

model {
  matrix[n_s,n_hp] m_zs; //a matrix for the subject-level z values of every hirearchical parameter
  vector[n_s*n_hp] v_zs; //a vector for the subject-level z values of every hierarchical parameter;
  
  //Choice model group-level priors
  // add weakly informative hyperpriors on the group-level parameters to the target density
  alpha_mu ~ normal(-.05,1.7); //This yields a nearly uniform 0-1 prior when logit transformed
  alpha_sigma ~ normal(0,4); //allow for large inter-subject heterogeneity without enforcing marginal subject-level priors that are excessively horeshoe-shaped
  tau_mu ~ normal(-.05,1.7);//setting tau priors same as alpha
  tau_sigma ~ normal(0,4);
  forget_mu ~ normal(-.05,1.7);//setting forgetting rate priors same as alpha
  forget_sigma ~ normal(0,4);
  
  beta_mu ~ normal(1,5); 
  beta_sigma ~ normal(0,5); 
  phi_mu ~ normal(2,10);
  phi_sigma ~ normal(0,10);
  
  ru_w_mu ~ normal(0,5);
  ru_w_sigma ~ normal(0,5);
  ev_w_mu ~ normal(0,5);
  ev_w_sigma ~ normal(0,5);
  vt_w_mu ~ normal(0,5);
  vt_w_sigma ~ normal(0,5);
  
  //Assign normal(0,1) priors to all z-scores. To speed sampling, stack z-scores into one big vector and use a sampling statement.
  m_zs = append_col(append_col(append_col(append_col(append_col(append_col(append_col(vt_w_z,ru_w_z),ev_w_z),alpha_z),beta_z),tau_z),phi_z),forget_z);
  v_zs = to_vector(m_zs); //convert matrix to vector
  v_zs ~ normal(0,1); //sample
  
  //Add the joint choice likelihood to the target density
  target += choice_lik;
}

