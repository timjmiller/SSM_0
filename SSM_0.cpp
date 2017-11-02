//This is a striped down version of the model used by Miller and Hyun (in press at CJFAS)
//and includes code to generated simulated data sets.

#include <TMB.hpp>
#include <iostream>

template <class Type>
matrix<Type> get_selblocks(int n_ages, int n_selblocks, int n_estimated_selpars, int n_other_selpars, vector<int> selblock_models, vector<int> estimated_pointers, vector<int> other_pointers, vector<Type> estimated_selpars, vector<Type> other_selpars, vector<Type> selpars_lower, vector<Type> selpars_upper)
{
  Type zero = Type(0);
  Type one = Type(1);
  Type half = Type(0.5);
  Type two = Type(2);
  int n_selpars = n_estimated_selpars + n_other_selpars;
  Type a50_1 = zero, k_1 = zero, a50_2 = zero, k_2 = zero;
  matrix<Type> selectivity_blocks(n_selblocks,n_ages);
  vector<Type> selpars(n_selpars);
  if(n_estimated_selpars>0) for(int i = 0; i < n_estimated_selpars; i++)
  {
    selpars(estimated_pointers(i)-1) = selpars_lower(i) + (selpars_upper(i) - selpars_lower(i))/(1 + exp(-estimated_selpars(i)));
  }
  if(n_other_selpars>0) for(int i = 0; i < n_other_selpars; i++) selpars(other_pointers(i)-1) = other_selpars(i);
  int count = 0;
  for(int i = 0; i < n_selblocks; i++) 
  {
    if (selblock_models(i)==1)
    { //proportions at age
      for(int a = 0; a < n_ages; a++) 
      {
        selectivity_blocks(i,a) = selpars(count);
        count++;
      }
    }
    else
    { //logistic or double-logistic
      a50_1 = selpars(count); // a50 parameter
      count++;
      k_1 = selpars(count); //  1/slope
      count++;
      if (selblock_models(i)==2) 
      { //increasing logistic
        Type age = zero;
        for (int a = 0; a < n_ages; a++) 
        {
          age += one;
          selectivity_blocks(i,a) = one/(one + exp(-(age - a50_1)/k_1));
        }
        for (int a = 0; a < n_ages; a++) selectivity_blocks(i,a) = selectivity_blocks(i,a)/selectivity_blocks(i,n_ages-1);
      }
      else
      { //double logistic
        if(selblock_models(i) == 3)
        {
          a50_2 = selpars(count);
          count++;
          k_2 = selpars(count);
          count++;
          Type age = zero;
          for (int a = 0; a < n_ages; a++)
          {
            age += one;
   	        selectivity_blocks(i,a) = one/(one + exp(-(age - a50_1)/k_1));
            selectivity_blocks(i,a) *= one/(one + exp((age - a50_2)/k_2)); //1-p
          }
        }
        else //model 4: declining logistic
        {
          Type age = zero;
          for (int a = 0; a < n_ages; a++) 
          {
            age += one;
            selectivity_blocks(i,a) = one/(one + exp((age - a50_1)/k_1));
          }
          for (int a = 0; a < n_ages; a++) selectivity_blocks(i,a) = selectivity_blocks(i,a)/selectivity_blocks(i,0);
        }
      }
    }
  }
  return selectivity_blocks;
}

template<class Type>
Type mydmultinom(int n_ages, Type Neff, vector<Type> paa_obs, vector<Type> paa_pred, int do_log)
{
  //multinomial
  vector<Type> n = Neff * paa_obs;
  Type ll = lgamma(Neff + Type(1.0));
  for(int a = 0; a < n_ages; a++) ll += -lgamma(n(a) + Type(1.0)) + n(a) * log(paa_pred(a));
  if(do_log == 1) return ll;
  else return exp(ll);
}
template<class Type>
vector<Type> rmultinom(int n_ages, Type Neff, vector<Type> paa_pred)
{
  //multinomial
  vector<Type> x(n_ages);
  int N = asInteger(asSEXP(Neff));
  vector<Type> cumsum(n_ages);
  x.setZero();
  for(int i = 0; i < N; i++)
  {
    Type y = runif(Type(0.0),Type(1.0));
    for(int a = 0; a < n_ages; a++) if(y < paa_pred.head(a+1).sum()) 
    {
      x(a) += Type(1.0);
      break;
    }
  }
  return x;
}
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n_years);
  DATA_INTEGER(n_ages);
  DATA_INTEGER(n_fleets);
  DATA_INTEGER(n_indices);
  DATA_INTEGER(n_selblocks);
  DATA_IVECTOR(selblock_models);
  DATA_IMATRIX(selblock_pointer_fleets);
  DATA_IMATRIX(selblock_pointer_indices);
  DATA_VECTOR(fracyr_SSB);
  DATA_MATRIX(mature);
  DATA_IVECTOR(waa_pointer_fleets);
  DATA_INTEGER(waa_pointer_totcatch);
  DATA_IVECTOR(waa_pointer_indices);
  DATA_INTEGER(waa_pointer_ssb);
  DATA_INTEGER(waa_pointer_jan1);
  DATA_ARRAY(waa);
  DATA_MATRIX(agg_catch);
  DATA_MATRIX(agg_catch_sigma);
  DATA_MATRIX(catch_paa);
  DATA_IMATRIX(use_catch_paa);
  DATA_MATRIX(catch_Neff);
  DATA_IVECTOR(units_indices);
  DATA_MATRIX(fracyr_indices);
  DATA_MATRIX(agg_indices);
  DATA_IMATRIX(use_indices);
  DATA_MATRIX(agg_index_sigma);
  DATA_IVECTOR(units_index_paa);
  DATA_MATRIX(index_paa);
  DATA_IMATRIX(use_index_paa);
  DATA_MATRIX(index_Neff);
  DATA_VECTOR(q_lower);
  DATA_VECTOR(q_upper);
  DATA_INTEGER(n_estimated_selpars);
  DATA_INTEGER(n_other_selpars);
  DATA_VECTOR(other_selpars);
  DATA_IVECTOR(estimated_selpar_pointers);
  DATA_IVECTOR(other_selpar_pointers);
  DATA_VECTOR(selpars_lower);
  DATA_VECTOR(selpars_upper);
  DATA_INTEGER(n_NAA_sigma);
  DATA_IVECTOR(NAA_sigma_pointers);
  DATA_IVECTOR(MAA_pointer); //n_ages
  DATA_IVECTOR(M_sigma_par_pointer); //n_M_re
  DATA_INTEGER(M_model); //0: just age-specific M, 1: Lorenzen M decline with age/weight
  DATA_INTEGER(N1_model); //0: just age-specific numbers at age, 1: 2 pars: log_N_{1,1}, log_F0, age-structure defined by equilibrium NAA calculations
  DATA_INTEGER(use_NAA_re);
  DATA_INTEGER(random_recruitment);
  
  PARAMETER_VECTOR(mean_rec_pars);
  PARAMETER_VECTOR(logit_q);
  PARAMETER_VECTOR(log_F1);
  PARAMETER_MATRIX(F_devs);
  PARAMETER_VECTOR(log_N1_pars); //length = n_ages or 2
  PARAMETER_VECTOR(log_NAA_sigma);
  PARAMETER_VECTOR(estimated_selpars);
  PARAMETER_MATRIX(log_NAA);
  PARAMETER_VECTOR(M_pars1);
  PARAMETER(log_b);
  PARAMETER_VECTOR(log_R); //n_years-1, if used.
  PARAMETER(log_R_sigma);

  Type zero = Type(0);
  Type one = Type(1);
  Type half = Type(0.5);
  Type two = Type(2);
  vector<int> any_index_age_comp(n_indices);
  vector<int> any_fleet_age_comp(n_fleets);
  for(int i = 0; i < n_indices; i++)
  {
    any_index_age_comp(i) = 0;
    for(int y = 0; y < n_years; y++) if(use_index_paa(y,i) == 1) any_index_age_comp(i) = 1;
  }
  for(int i = 0; i < n_fleets; i++)
  {
    any_fleet_age_comp(i) = 0;
    for(int y = 0; y < n_years; y++) if(use_catch_paa(y,i) == 1) any_fleet_age_comp(i) = 1;
  }
  vector<Type> sigma2_log_NAA = exp(log_NAA_sigma*two);
  matrix<Type> MAA(n_years,n_ages);
  for(int i = 0; i < n_ages; i++) 
  {
    if(M_model == 0) //M by age
    {
      MAA(0,i) = exp(M_pars1(MAA_pointer(i)-1));
      for(int y = 1; y < n_years; y++) MAA(y,i) = MAA(0,i);
    }
    else //M_model == 1, Lorenzen
    {
      MAA(0,i) = exp(M_pars1(0) - exp(log_b) * log(waa(waa_pointer_jan1-1,0,i)));
      for(int y = 1; y < n_years; y++) MAA(y,i) = exp(M_pars1(0) - exp(log_b) * log(waa(waa_pointer_jan1-1,y,i)));
    }
  }
  vector<Type> SSB(n_years);
  vector<Type> log_SSB(n_years);
  matrix<Type> F(n_years,n_fleets);
  matrix<Type> log_F(n_years,n_fleets);
  array<Type> pred_CAA(n_years,n_fleets,n_ages);
  array<Type> pred_catch_paa(n_years,n_fleets,n_ages);
  matrix<Type> pred_catch(n_years,n_fleets);
  array<Type> pred_IAA(n_years,n_indices,n_ages);
  array<Type> pred_index_paa(n_years,n_indices,n_ages);
  matrix<Type> pred_indices(n_years,n_indices);
  matrix<Type> log_pred_catch(n_years,n_fleets);
  matrix<Type> NAA(n_years,n_ages);
  matrix<Type> pred_NAA(n_years,n_ages);
  array<Type> FAA(n_years,n_fleets,n_ages);
  matrix<Type> FAA_tot(n_years,n_ages);
  matrix<Type> ZAA(n_years,n_ages);
  array<Type> QAA(n_years,n_indices,n_ages);
  matrix<Type> selblocks(n_selblocks,n_ages);
  vector<Type> q(n_indices);
  Type nll = zero; //negative log-likelihood
  vector<Type> t_paa(n_ages), t_pred_paa(n_ages);

  selblocks = get_selblocks(n_ages, n_selblocks, n_estimated_selpars, n_other_selpars, selblock_models, estimated_selpar_pointers, other_selpar_pointers, estimated_selpars, other_selpars, selpars_lower, selpars_upper);
  for(int i = 0; i < n_indices; i++)
  {
    q(i) = q_lower(i) + (q_upper(i) - q_lower(i))/(1 + exp(-logit_q(i)));
    for(int y = 0; y < n_years; y++) 
    {
      for(int a = 0; a < n_ages; a++) QAA(y,i,a) = q(i) * selblocks(selblock_pointer_indices(y,i)-1,a);
    }
  }
  FAA_tot.setZero();
  for(int f = 0; f < n_fleets; f++)
  {
    log_F(0,f) = log_F1(f);
    F(0,f) = exp(log_F(0,f));
    for(int a = 0; a < n_ages; a++) 
    {
      FAA(0,f,a) = F(0,f) * selblocks(selblock_pointer_fleets(0,f)-1,a);
      FAA_tot(0,a) = FAA_tot(0,a) + FAA(0,f,a);
    }
    for(int y = 1; y < n_years; y++) 
    {
      log_F(y,f) = log_F(y-1,f) + F_devs(y-1,f);
      F(y,f) = exp(log_F(y,f));
      for(int a = 0; a < n_ages; a++) 
      {
        FAA(y,f,a) = F(y,f) * selblocks(selblock_pointer_fleets(y,f)-1,a);
        FAA_tot(y,a) = FAA_tot(y,a) + FAA(y,f,a);
      }
    }
  }
  ZAA = FAA_tot + MAA;
  
  SSB.setZero();
  //year 1, no random effects.
  for(int a = 0; a < n_ages; a++) 
  {
    pred_NAA(0,a) = NAA(0,a) = exp(log_N1_pars(a));
    SSB(0) += NAA(0,a) * waa(waa_pointer_ssb-1,0,a) * mature(0,a) * exp(-ZAA(0,a)*fracyr_SSB(0));
  }
  
  //after year 1
  Type nll_NAA = zero;
  Type nll_recruit = zero;
  for(int y = 1; y < n_years; y++) 
  {
    //expected recruitment
    pred_NAA(y,0) = exp(mean_rec_pars(0)); //random about mean, whether or not using NAA_re.
    
    //random effects NAA, state-space model for all numbers at age
    if(use_NAA_re == 1)
    {
      SIMULATE log_NAA(y-1,0) = rnorm(log(pred_NAA(y,0)), exp(log_NAA_sigma(NAA_sigma_pointers(0)-1)));
      NAA(y,0) = exp(log_NAA(y-1,0));
      //expected numbers at age after recruitment
      for(int a = 1; a < n_ages-1; a++) 
      {
        pred_NAA(y,a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
        SIMULATE log_NAA(y-1,a) = rnorm(log(pred_NAA(y,a)), exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)));
        NAA(y,a) = exp(log_NAA(y-1,a));
      }
      pred_NAA(y,n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));
      SIMULATE log_NAA(y-1,n_ages-1) = rnorm(log(pred_NAA(y,n_ages-1)), exp(log_NAA_sigma(NAA_sigma_pointers(n_ages-1)-1)));
      NAA(y,n_ages-1) = exp(log_NAA(y-1,n_ages-1));
      for(int a = 0; a < n_ages; a++) nll_NAA -= dnorm(log(NAA(y,a)), log(pred_NAA(y,a)), exp(log_NAA_sigma(NAA_sigma_pointers(a)-1)), 1);
    }

    //no random effects for NAA, SCAA-like model, recruitments are either fixed or random effects
    if(use_NAA_re != 1) 
    {
      if(random_recruitment == 1) //estimate recruitment as random effects, otherwise as fixed effects. pred_NAA(y,0) should be properly specified above in any case.
      {
        nll_recruit -= dnorm(log_R(y-1), log(pred_NAA(y,0)), exp(log_R_sigma), 1);
        SIMULATE log_R(y-1) = rnorm(log(pred_NAA(y,0)), exp(log_R_sigma));
        //when random effects not used for all numbers at age, survival is deterministic.
        //NAA for recruits is already filled in above in this case
      }
      NAA(y,0) = exp(log_R(y-1));
      //expected numbers at age after recruitment
      for(int a = 1; a < n_ages-1; a++) pred_NAA(y,a) = NAA(y,a) = NAA(y-1,a-1) * exp(-ZAA(y-1,a-1));
      pred_NAA(y,n_ages-1) = NAA(y,n_ages-1) = NAA(y-1,n_ages-2) * exp(-ZAA(y-1,n_ages-2)) + NAA(y-1,n_ages-1) * exp(-ZAA(y-1,n_ages-1));
    }
    for(int a = 0; a < n_ages; a++) SSB(y) += NAA(y,a) * waa(waa_pointer_ssb-1,y,a) * mature(y,a) * exp(-ZAA(y,a)*fracyr_SSB(y));
  }
  nll += nll_NAA;
  nll += nll_recruit;
  SIMULATE
  {
    if(use_NAA_re == 1) REPORT(log_NAA);
    if(random_recruitment == 1) REPORT(log_R);
  }

  vector<Type> nll_agg_catch(n_fleets), nll_catch_acomp(n_fleets);
  nll_agg_catch.setZero(); nll_catch_acomp.setZero();
  for(int y = 0; y < n_years; y++)
  {
    int acomp_par_count = 0;
    for(int f = 0; f < n_fleets; f++)
    {
      pred_catch(y,f) = zero;
      Type tsum = zero;
      for(int a = 0; a < n_ages; a++) 
      {
        pred_CAA(y,f,a) =  NAA(y,a) * FAA(y,f,a) * (1 - exp(-ZAA(y,a)))/ZAA(y,a);
        pred_catch(y,f) += waa(waa_pointer_fleets(f)-1,y,a) * pred_CAA(y,f,a);
        tsum += pred_CAA(y,f,a);
      }
      nll_agg_catch(f) -= dnorm(log(agg_catch(y,f)), log(pred_catch(y,f)), agg_catch_sigma(y,f),1);
      SIMULATE agg_catch(y,f) = exp(rnorm(log(pred_catch(y,f)), agg_catch_sigma(y,f)));
      if(any_fleet_age_comp(f) == 1)
      {
        if(use_catch_paa(y,f) == 1) 
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_catch_paa(y,f,a) = pred_CAA(y,f,a)/tsum;
            t_pred_paa(a) = pred_catch_paa(y,f,a);
            t_paa(a) = catch_paa(f * n_years + y,a);
          }
          nll_catch_acomp(f) -= mydmultinom(n_ages, catch_Neff(y,f), t_paa, t_pred_paa,1);
          SIMULATE
          {
            t_paa = rmultinom(n_ages, catch_Neff(y,f), t_pred_paa)/catch_Neff(y,f);
            for(int a = 0; a < n_ages; a++) catch_paa(f * n_years + y, a) = t_paa(a);
          }
        }
      }
    }
  }
  nll += nll_agg_catch.sum();
  nll += nll_catch_acomp.sum();
  SIMULATE
  {
    REPORT(catch_paa);
    REPORT(agg_catch);  
  }
  vector<Type> nll_agg_indices(n_indices), nll_index_acomp(n_indices);
  nll_agg_indices.setZero(); nll_index_acomp.setZero(); pred_indices.setZero();
  for(int y = 0; y < n_years; y++)
  {
    int acomp_par_count = 0;
    for(int i = 0; i < n_indices; i++) 
    {
      Type tsum = zero;
      for(int a = 0; a < n_ages; a++) 
      {
        pred_IAA(y,i,a) =  NAA(y,a) * QAA(y,i,a) * exp(-ZAA(y,a) * fracyr_indices(y,i));
        if(units_indices(i) == 1) pred_indices(y,i) += waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        else pred_indices(y,i) += pred_IAA(y,i,a);
      }
      for(int a = 0; a < n_ages; a++) 
      {
        if(units_index_paa(i) == 1) pred_IAA(y,i,a) = waa(waa_pointer_indices(i)-1,y,a) * pred_IAA(y,i,a);
        tsum += pred_IAA(y,i,a);
      }
      
      if(use_indices(y,i) == 1)
      {
        nll_agg_indices(i) -= dnorm(log(agg_indices(y,i)), log(pred_indices(y,i)), agg_index_sigma(y,i), 1);
        SIMULATE agg_indices(y,i) = exp(rnorm(log(pred_indices(y,i)), agg_index_sigma(y,i)));
      }
      if(any_index_age_comp(i) == 1)
      {
        if(use_index_paa(y,i) > 0)
        {
          for(int a = 0; a < n_ages; a++)
          {
            pred_index_paa(y,i,a) = pred_IAA(y,i,a)/tsum;
            t_pred_paa(a) = pred_index_paa(y,i,a);
            t_paa(a) = index_paa(i * n_years + y,a);
          }
          nll_index_acomp(i) -= mydmultinom(n_ages, index_Neff(y,i), t_paa, t_pred_paa,1);
          SIMULATE
          {
            t_paa = rmultinom(n_ages, index_Neff(y,i), t_pred_paa)/index_Neff(y,i);
            for(int a = 0; a < n_ages; a++) index_paa(i * n_years + y, a) = t_paa(a);
          }
        }
      }
    }
  }
  nll += nll_agg_indices.sum();
  nll += nll_index_acomp.sum();
  SIMULATE
  {
    REPORT(index_paa);
    REPORT(agg_indices);  
  }
  
  log_SSB = log(SSB);
  
  REPORT(NAA);
  REPORT(MAA);
  REPORT(F);
  REPORT(FAA_tot);
  REPORT(pred_NAA);
  REPORT(SSB);
  REPORT(selblocks);
  REPORT(q);
  ADREPORT(F);
  ADREPORT(log_SSB);
  
  return nll;
}

