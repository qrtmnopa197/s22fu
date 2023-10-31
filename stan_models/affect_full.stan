// Q-learning model with choice autocorrelation, plus an affect model that
// predicts valence using the full set of predictors.

data {
  //Choice model data
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
 
 
  //Affect model data
  
  //For V calculation
  int<lower=1> n_tb; // number of trials in each block
  array[n_s] vector[n_t] t_in_b; //trial-in-block
   
  array[n_s,n_t] int chosen_frac; //contains the index of the chosen fractal for each trial
  array[n_s,n_t] int unchosen_frac; //ditto unchosen fractals     
  
  array[n_s] vector[n_t] chosen_out; //the outcomes of the chosen fractals
  array[n_s] vector[n_t] unchosen_out;  //ditto unchosen fractals

  array[n_s] vector[n_t] block_c; //block number, centered, for each trial
  array[n_s] vector[n_t] tinb_c; //trial in block, centered

  array[n_s,n_t] int dp_num; //overall number of decision probe
  array[n_s,n_t] int fp_num; //overall number of feedback probe
  
  int n_dp; // total number of decision probes completed across all subjects
  int n_fp; // ditto feedback probes
  
  vector[n_dp] dec_rate;
  vector[n_fp] feed_rate;
  
  array[n_s] int decrate_first; // 1/0 indicator for each subject indicating whether they made decision ratings in the first block
}

transformed data{
  int n_dpi = 3; //Number of decision preditors of interest
  int n_fpi = 4; //Number of feedback predictors of interest
  int n_hp = n_dpi + n_fpi + 7; //number of hirearchical parameters
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
  
  
  //Affect parameters
  
  //baseline valence 
  real base_mu; 
  real<lower=0> base_sigma; 
  vector[n_s] base_z;
  
  //effect of block number
  real blk_num_mu;
  //effect of making a decision vs. feedback rating
  real dec_mu;
  real<lower=0> bnum_dec_sigma; //total SD for blk_num and dec effects
  vector[n_s] blk_z; //overall effect of block for the subject

  //effect of trial number
  real tr_mu;
  real<lower=0> tr_sigma;
  vector[n_s] tr_z;
  
  //decision predictors of interest weights
  vector[n_dpi] dpiw_mu;
  vector<lower=0>[n_dpi] dpiw_sigma;
  matrix[n_s,n_dpi] dpiw_z;
  
  //feedback predictors of interest weights
  vector[n_fpi] fpiw_mu;
  vector<lower=0>[n_fpi] fpiw_sigma;
  matrix[n_s,n_fpi] fpiw_z;
  
  //Predictors of interest; Bayesian ridge shrinkage prior is applied to these weights
  real<lower=0> shrink; //sigma on the shrinkage prior
  
  real<lower=0> d_resid;  //residual variance in decision ratings
  real<lower=0> f_resid;  //ditto feedback ratings
}


transformed parameters {
  vector[n_s*n_t] choice_lik;
  vector[n_dp] dec_pred; //predicted post-decision valence ratings
  vector[n_fp] feed_pred; //predicted post-feedback ratings
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q; // The Q values. These are a 2d array containing vectors: the first dimension is
                                  // the subject, the second dimension is the trial, and for each trial there's a
                                  // vector with a Q value for each fractal
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value, of the same structure as Q
    array[n_s] vector[n_t] V; //V_block (block-level expected value)

    
    array[n_s] vector[n_t] PE_a;  // The RPEs to outcome a. This is a 1d array containing one vector per subject,
                                  // where that vector contains one value per trial
    array[n_s] vector[n_t] PE_b;  // diddo outcome B
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    
    vector[n_s] beta = beta_mu + beta_sigma*beta_z;  // get subject-level beta values using n.c.p.
    vector[n_s] phi = phi_mu + phi_sigma*phi_z; //diddo phi
    vector[n_s] alpha = inv_logit(alpha_mu + alpha_sigma*alpha_z); //diddo alpha. Also, use inv_logit to get alpha between 0 and 1
    vector[n_s] tau = inv_logit(tau_mu + tau_sigma*tau_z); //likewise for tau
    
    //Affect model parameters
    row_vector[n_dpi] dpi; //decision predictors of interest on the current trial
    row_vector[n_fpi] fpi; //ditto feedback
    
    vector[n_s] base = base_mu + base_sigma*base_z;
   
    vector[n_s] tr = tr_mu + tr_sigma*tr_z; //effect of trial
    vector[n_s] blk; //effect of block
    matrix[n_s,n_dpi] dpiw; //weights for decision predictors of interest
    matrix[n_s,n_fpi] fpiw; //ditto feedback
    
    
    //to set the weight on the block predictor, loop through each subject...
    for(s in 1:n_s){
      if(decrate_first[s] == 0){
        blk[s] = blk_num_mu + dec_mu + bnum_dec_sigma*blk_z[s]; //if this subject completed decision ratings in the second block, then the decision effect should be added to the block effect
      } else if(decrate_first[s] == 1){
        blk[s] = blk_num_mu - dec_mu + bnum_dec_sigma*blk_z[s]; //if the subject completed decision ratings in the first block, the decision effect should be subtracted
      }
    }
    
    for(w in 1:n_dpi){
      dpiw[,w] = dpiw_mu[w] + dpiw_sigma[w]*dpiw_z[,w];
    }
    for(w in 1:n_fpi){
      fpiw[,w] = fpiw_mu[w] + fpiw_sigma[w]*fpiw_z[,w];
    }
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          Q[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          C[s,t] = rep_vector(0,n_f); // diddo C values
        }
        if(t_in_b[s,t] == 1){
          V[s,t] = 0; //if this is the first trial of the block, then V_block should be 0
        }
        
        choice_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | beta[s]*Q[s,t,{fA[s,t],fB[s,t]}] + phi[s]*C[s,t,{fA[s,t],fB[s,t]}]); //get the likelihood of the choice on this trial. Doing the log lik
                                                                                                                                                                        //here prevents vectorization, but you have to calculate a vector
                                                                                                                                                                        //of log lik for loo_CV anyway. Plus, it can't be done for choice data.
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
          
          V[s,t+1] = V[s,t] + alpha[s]*(chosen_out[s,t] - V[s,t]); //update V_block with whatever the outcome was
        }
        if(dp_num[s,t] != 0){
          dpi[1] = V[s,t]; //assign V_block
          dpi[2] = Q[s,t,chosen_frac[s,t]];
          dpi[3] = Q[s,t,unchosen_frac[s,t]];
          dec_pred[dp_num[s,t]] = base[s] + blk[s]*block_c[s,t] + tr[s]*tinb_c[s,t] + dot_product(dpiw[s],dpi);
        } else if(fp_num[s,t] != 0){
          fpi[1] = chosen_out[s,t]*exp(choice_lik[(s-1)*n_t+t]) + unchosen_out[s,t]*(1-exp(choice_lik[(s-1)*n_t+t]));
          fpi[2] = Q[s,t,chosen_frac[s,t]];
          fpi[3] = chosen_out[s,t];
          fpi[4] = unchosen_out[s,t];
          feed_pred[fp_num[s,t]] = base[s] + blk[s]*block_c[s,t] + tr[s]*tinb_c[s,t] + dot_product(fpiw[s],fpi);
        }
      }
    }
  }//anonymous_scope_end
}

model {
  matrix[n_s,n_hp] m_zs; //a matrix for the subject-level z values of every hirearchical parameter
  vector[n_s*n_hp] v_zs; //a vector for the subject-level z values of every hierarchical parameter;
  
  //Choice model group-level priors
  // add the hyperpriors on the group-level parameters to the target density
  alpha_mu ~ normal(-.05,1.7); //This yields a nearly uniform 0-1 prior when logit transformed
  alpha_sigma ~ normal(0,4); //allow for large inter-subject heterogeneity without enforcing marginal subject-level priors that are excessively horeshoe-shaped
  beta_mu ~ normal(1,5); // Weakly informative prior
  beta_sigma ~ normal(0,5); //Setting the prior on the SD with similar logic
  tau_mu ~ normal(-.05,1.7); //setting tau priors same as alpha
  tau_sigma ~ normal(0,4); // ""
  phi_mu ~ normal(2,10); // Setting a reguarlizing prior. Matching the phi prior to the beta prior, assuming that uatocorrelation and reward are apt to have 
                         // equivalent effects on choice. But differences in phi between fA and fB will be about .5 at asymptote (based on data in which
                         // subjects tend to repeat choices 75% of the time), whereas it's more like 1.2 for beta, so the prior for phi will be set a bit wider.
  phi_sigma ~ normal(0,10);
  
  //Affect model group-level priors
  base_mu ~ normal(0,3); //set prior for ad intercept; weakly informative - 1SD mean that baseline valal is about as highest/lowest affect the subject experiences in the task
  base_sigma ~ normal(0,5); //1 SD out suggests that if the average subject is near an affective low at baseline, some subjects may still be near an affective high, more or less
  dec_mu ~ normal(0,3);
  blk_num_mu ~ normal(0,3);
  bnum_dec_sigma ~ normal(0,8);
  tr_mu ~ normal(0,0.1);
  tr_sigma ~ normal(0,0.2);
  
  d_resid ~ normal(0,2); //1SD out would mean the model explains near-0 variability, which it well may. Extending to 2 since 1SD is a fairly typical case.
  f_resid ~ normal(0,2);

  dpiw_mu ~ normal(0,shrink);//weakly informative
  dpiw_sigma ~ normal(0,1);
  
  fpiw_mu ~ normal(0,shrink);
  fpiw_sigma ~ normal(0,1);
  
  shrink ~ normal(0,1);
  
  //Assign normal(0,1) priors to all z-scores. To speed sampling, stack z-scores into one big vector and use a sampling statement.
  m_zs = append_col(append_col(append_col(append_col(append_col(append_col(append_col(append_col(dpiw_z,fpiw_z),base_z),blk_z),alpha_z),beta_z),tau_z),phi_z),tr_z);
  v_zs = to_vector(m_zs); //convert matrix to vector
  v_zs ~ normal(0,1); //sample
  
  //Add the joint choice likelihood to the target density
  target += choice_lik;
  
  //Add the affect prediction likelihoods to the target density
  dec_rate ~ normal(dec_pred, d_resid);
  feed_rate ~ normal(feed_pred, f_resid);
}

generated quantities{
  vector[n_dp+n_fp] affect_lik; //log likelihoods of all affect rating
  vector[n_s*n_t+n_dp+n_fp] log_lik; //log likelihoods of all choices followed by log likelihoods of all affect ratings
  {//anonymous scope start
    vector[n_dp] dec_lik; //log liks of post-dec ratings
    vector[n_fp] feed_lik; //ditto post-feed
    
    for(i in 1:n_dp){
      dec_lik[i] = normal_lpdf(dec_rate[i] | dec_pred[i], d_resid);
    }
    for(i in 1:n_fp){
      feed_lik[i] = normal_lpdf(feed_rate[i] | feed_pred[i], f_resid);
    }
    
    affect_lik = append_row(dec_lik,feed_lik);
    
    log_lik = append_row(choice_lik,affect_lik);
  }//anonymous scope end
}

