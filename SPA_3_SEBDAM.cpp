#include <TMB.hpp>

template <class Type>
Type sqr(Type x) {
  return x * x;
}

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  //Data
  DATA_VECTOR(I); //commercial biomass survey index by tow (dim n_i)
  DATA_VECTOR(IR); //recruit survey index by tow (dim n_i)
  DATA_VECTOR(area); //area covered by each knot (dim n_s)
  DATA_MATRIX(C); //commercial catch (dim n_s,n_t)
  DATA_VECTOR(L); //number of clappers (dim n_i)
  DATA_VECTOR(n_bin); //number of shell caught in tow (dim n_i)
  DATA_VECTOR(n_tows); //number of tows in year (dim n_t)
  DATA_VECTOR(pos_tows_I); //number of tows that captured commercial biomass (dim n_t)
  DATA_VECTOR(pos_tows_IR); //number of tows that captured recruits (dim n_t)
  
  //Sizes
  DATA_INTEGER(n_i); //Number of observations per year
  DATA_INTEGER(n_t); //Number of years
  DATA_INTEGER(n_s); //Number of knots
  DATA_INTEGER(n_m); //number of vertex in mesh
  
  //Indices
  DATA_FACTOR(s_i); //Indexing for knot (dim n_i)
  DATA_FACTOR(t_i); //Indexing for year (dim n_i)
  DATA_FACTOR(v_i); //Indexing for specific mesh location to match knots with right vertex (dim n_m)
  
  //SPDE objects
  DATA_STRUCT(spde_aniso, spde_aniso_t); //INLA anisotropic mesh structure
  
  //Fixed covariates
  DATA_VECTOR(gI); //commercial biomass growth
  DATA_VECTOR(gR); //recruitment growth
  
  //Include prior
  DATA_INTEGER(prior_q);
  
  //Parameters
  PARAMETER(log_sigma_epsilon); //obs sd survey index commercial biomass
  PARAMETER(log_sigma_upsilon); //obs sd survey index recruits
  PARAMETER(log_S); //clapper catchability
  PARAMETER(log_kappa_B);//commercial size range parameter
  PARAMETER(log_tau_B); //commercial size spatio-temporal variability parameter
  PARAMETER(log_kappa_R); //recruit range parameter
  PARAMETER(log_tau_R); //recruit spatio-temporal variability parameter
  PARAMETER(log_kappa_m); //mortality range parameter
  PARAMETER(log_tau_m); //mortality spatio-temporal variability parameter
  PARAMETER(log_R0);//initial recruit mean value
  PARAMETER(log_B0);//initial commercial biomass mean value
  PARAMETER(log_m0);//initial mortality mean value
  PARAMETER(logit_p_I); //probability of capturing commercial biomass
  PARAMETER(logit_p_IR); //probability of capturing recruits
  PARAMETER_VECTOR(log_qI); //commerical biomass catchability
  // PARAMETER(log_qI);
  PARAMETER(log_qR); //recruit catchability
  PARAMETER_VECTOR(log_H_input_B); //commercial size anisotropy parameters
  PARAMETER_VECTOR(log_H_input_m); //mortality anisotropy parameters
  
  Type pi = 3.141593;
  Type Range_B = sqrt(8)/exp(log_kappa_B);
  Type Range_R = sqrt(8)/exp(log_kappa_R); 
  Type Range_m = sqrt(8)/exp(log_kappa_m);
  Type sigma_epsilon = exp(log_sigma_epsilon);
  Type sigma_upsilon = exp(log_sigma_upsilon);
  Type S = exp(log_S); 
  Type qR = exp(log_qR);
  // Type qI = exp(log_qI);
  vector <Type> qI(n_s); qI = exp(log_qI);
  Type kappa_B = exp(log_kappa_B);
  Type tau_B = exp(log_tau_B);
  Type kappa_R = exp(log_kappa_R);
  Type tau_R = exp(log_tau_R);
  Type kappa_m = exp(log_kappa_m);
  Type tau_m = exp(log_tau_m);
  Type R0 = exp(log_R0);
  Type B0 = exp(log_B0);
  Type m0 = exp(log_m0);
  Type p_I = invlogit(logit_p_I);
  Type p_IR = invlogit(logit_p_IR);
  
  // Calculate marginal field variances based on relation described in Lindgren et al., 2012
  Type SigmaO_B = 1 / (sqrt(4*pi)*exp(log_tau_B)*exp(log_kappa_B));
  Type SigmaO_R = 1 / (sqrt(4*pi)*exp(log_tau_R)*exp(log_kappa_R));
  Type SigmaO_m = 1 / (sqrt(4*pi)*exp(log_tau_m)*exp(log_kappa_m));
  
  //Random Effects
  PARAMETER_ARRAY(omega_B);
  PARAMETER_ARRAY(omega_R);
  PARAMETER_ARRAY(omega_m);
  
  // Set up matrices for processes of interest
  matrix <Type> log_B(n_s,(n_t+1));
  matrix <Type> B(n_s,(n_t+1));
  matrix <Type> resid_B(n_s,(n_t+1));
  matrix <Type> log_areaB(n_s,(n_t+1));
  matrix <Type> log_areaR(n_s,(n_t+1));
  matrix <Type> areaB(n_s,(n_t+1));
  matrix <Type> areaR(n_s,(n_t));
  matrix <Type> log_R(n_s,(n_t));
  matrix <Type> R(n_s,(n_t));
  matrix <Type> log_m(n_s,(n_t+1));
  matrix <Type> m(n_s,(n_t+1));
  vector <Type> resid_I(n_i);
  Type sum_area; sum_area = 0;
  //
  
  
  //Set up initial states and other elements
  SIMULATE{
    n_bin = rpois(Type(100));
    n_tows = Type(120);
    for (int t = 0; t < n_t; t++){
      gI(t) = 1.1;
      gR(t) = 1.5;
    }
    REPORT(gI);
    REPORT(gR);
    REPORT(n_bin);
    REPORT(n_tows);
    
    for (int s = 0; s < n_s; s++){
      sum_area += area(s);
    }
    
  }
  
  //Setup for simulations and derived values
  vector <Type> bern_I(n_i);
  vector <Type> bern_IR(n_i);
  vector <Type> totB(n_t+1); totB.setZero();
  vector <Type> totR(n_t); totR.setZero();
  vector <Type> log_totB(n_t+1);
  vector <Type> log_totR(n_t);
  vector <Type> mean_m(n_t+1); mean_m.setZero();
  matrix <Type> mean_pro_m(n_s,n_t+1);
  matrix <Type> mean_pro_B(n_s,n_t+1);
  // vector <Type> no_spat_m(n_t+1); no_spat_m.setZero();
  // vector <Type> no_spat_R(n_t); no_spat_R.setZero();
  // vector <Type> no_spat_B(n_t+1); no_spat_B.setZero();
  // vector <Type> no_spat_C(n_t+1); no_spat_C.setZero();
  
  // ----------------------------------------------
  // nll
  vector <Type> nll_comp(14); nll_comp.setZero();
  
  
  //Anisotropic fields, so must have a matrix H
  //Parameterize so that det(H)=1 to preserve volume as seen in Lindgren 2011
  //commercial size anisotropy matrix
  matrix<Type> H_B(2,2);
  H_B(0,0) = exp(log_H_input_B(0));
  H_B(1,0) = log_H_input_B(1);
  H_B(0,1) = log_H_input_B(1);
  H_B(1,1) = (1+log_H_input_B(1)*log_H_input_B(1)) / exp(log_H_input_B(0));
  REPORT(H_B);
  
  //mortality anisotropy matrix
  matrix<Type> H_m(2,2);
  H_m(0,0) = exp(log_H_input_m(0));
  H_m(1,0) = log_H_input_m(1);
  H_m(0,1) = log_H_input_m(1);
  H_m(1,1) = (1+log_H_input_m(1)*log_H_input_m(1)) / exp(log_H_input_m(0));
  REPORT(H_m);
  
  //Set up GMRF for commercial size biomass
  SparseMatrix <Type> Q_B = Q_spde(spde_aniso, kappa_B, H_B);
  for (int t = 0; t < (n_t+1); t++){
    if (t == 0) {nll_comp(0) += SCALE(GMRF(Q_B),1/tau_B)(omega_B.col(t));}
    if (t >= 1) {nll_comp(0) += SCALE(GMRF(Q_B),1/tau_B)(omega_B.col(t));}
  }
  
  SIMULATE {
    SparseMatrix <Type> Q_B = Q_spde(spde_aniso, kappa_B, H_B);
    density::GMRF_t<Type> gmrf1(Q_B);
    for (int t = 0; t < (n_t+1); t++){
      vector <Type> temp_omega(n_m);
      SCALE(gmrf1,1/tau_B).simulate(temp_omega);
      omega_B.col(t)=temp_omega;
    }
    REPORT(Q_B);
    REPORT(H_B);
    REPORT(omega_B);
  }
  
  //Set up GMRF for recruits
  SparseMatrix <Type> Q_R = Q_spde(spde_aniso, kappa_R, H_B);
  for (int t = 0; t < (n_t); t++){
    nll_comp(1) += SCALE(GMRF(Q_R),1/tau_R)(omega_R.col(t));
  }
  
  SIMULATE {
    SparseMatrix <Type> Q_R = Q_spde(spde_aniso, kappa_R, H_B);
    density::GMRF_t<Type> gmrf2(Q_R);
    for (int t = 0; t < n_t; t++){
      vector <Type> temp_omega(n_m);
      SCALE(gmrf2,1/tau_R).simulate(temp_omega);
      omega_R.col(t)=temp_omega;
    }
    REPORT(Q_R);
    REPORT(omega_R);
  }
  
  //Set up GMRF for mortality
  SparseMatrix <Type> Q_m = Q_spde(spde_aniso, kappa_m, H_m);
  // nll_comp(11) += SCALE(GMRF(Q_m),1/tau_m)(balurp_m);
  for (int t = 0; t < (n_t+1); t++){
    if (t == 0) {nll_comp(2) += SCALE(GMRF(Q_m),1/tau_m)(omega_m.col(t));}
    if (t >= 1) {nll_comp(2) += SCALE(GMRF(Q_m),1/tau_m)(omega_m.col(t));}
  }
  
  SIMULATE {
    SparseMatrix <Type> Q_m = Q_spde(spde_aniso, kappa_m, H_m);
    density::GMRF_t<Type> gmrf3(Q_m);
    for (int t = 0; t < (n_t+1); t++){
      vector <Type> temp_omega(n_m);
      SCALE(gmrf3,1/tau_m).simulate(temp_omega);
      omega_m.col(t)=temp_omega;
    }
    REPORT(Q_m);
    REPORT(H_m);
    REPORT(omega_m);
  }
  
  
  //Derived values
  //For mortality
  for (int s = 0; s < n_s; s++){
    log_m(s,0) = log(m0 * exp(omega_m(v_i(s),0)));
    m(s,0) = exp(log_m(s,0));
    for (int t = 1; t < (n_t); t++){
      log_m(s,t) = log( m(s,t-1) * exp(omega_m(v_i(s),t)) );
      m(s,t) = exp(log_m(s,t));
    }
    //Predict 1 year ahead, assume same spatial pattern for error but propagate error
    log_m(s,n_t) = log( m(s,n_t-1)*exp(omega_m(v_i(s),n_t)));
    m(s,n_t) = exp(log_m(s,n_t));
  }
  
  //Simulate Natural mortality
  SIMULATE {
    // no_spat_m(0) = m0*exp(rnorm(Type(0.0),SigmaO_m));
    // for (int t = 1; t < (n_t+1); t++){
    //   no_spat_m(t) = no_spat_m(t-1)*exp(rnorm(Type(0.0),SigmaO_m));
    // }
    
    
    for (int s = 0; s < n_s; s++){
      m(s,0) = m0 * exp(omega_m(v_i(s),0));
      // m(s,0) = m0 * exp(rnorm(Type(0.0),SigmaO_m));
      // m(s,0) = no_spat_m(0);
      log_m(s,0) = log(m(s,0));
      for (int t = 1; t < (n_t); t++){
        // m(s,t) = no_spat_m(t);
        m(s,t) = m(s,t-1) * exp(omega_m(v_i(s),t)) ;
        // m(s,t) = m(s,t-1) * exp(rnorm(Type(0.0),SigmaO_m));
        log_m(s,t) = log(m(s,t));
      }
      // m(s,n_t) = no_spat_m(n_t);
      m(s,n_t) = m(s,n_t-1) * exp(omega_m(v_i(s),n_t));
      // m(s,n_t) = m(s,n_t-1) * exp(rnorm(Type(0.0),SigmaO_m));
      // pred_m(s) = log(m(s,n_t-1)) - sqr(SigmaO_m)/2.0 +rnorm(Type(0.0),SigmaO_m);
      // m(s,n_t) = exp(pred_m(s));
      // m(s,n_t) = m(s,n_t-1) * exp(omega_m(v_i(s),n_t-1)) ;
      log_m(s,n_t) = log(m(s,n_t));
    }
    // REPORT(pred_m);
    // REPORT(no_spat_m);
    REPORT(m);
    REPORT(log_m);
  }
  
  //Recruit derivation
  for (int s = 0; s < n_s; s++){
    log_R(s,0) = log(R0 * exp(omega_R(v_i(s),0)));
    R(s,0) = exp(log_R(s,0));
    for (int t = 1; t < (n_t); t++){
      log_R(s,t) = log(R(s,t-1)*exp(omega_R(v_i(s),t)));
      R(s,t) = exp(log_R(s,t));
    }
  }
  
  //Simulate recruitment
  SIMULATE{
    // no_spat_R(0) = (R0*sum_area/1000);
    // for (int t = 1; t < n_t; t++){
    //   no_spat_R(t) = no_spat_R(t-1)*exp(rnorm(Type(0.0),SigmaO_R));
    // }
    for (int s = 0; s < n_s; s++){
      // R(s,0) = no_spat_R(0)/sum_area*1000;
      R(s,0) = R0 * exp(omega_R(v_i(s),0));
      // R(s,0) = R0 * exp(rnorm(Type(0.0),SigmaO_R));
      log_R(s,0) = log(R(s,0));
      for (int t = 1; t < (n_t); t++){
        // R(s,t) = no_spat_R(t)/sum_area*1000;
        R(s,t) = R(s,t-1)*exp(omega_R(v_i(s),t));
        // R(s,t) = R(s,t-1) * exp(rnorm(Type(0.0),SigmaO_R));
        log_R(s,t) = log(R(s,t));
      }
    }
    REPORT(R);
    REPORT(log_R);
  }
  
  
  //Biomass derivation
  for (int s = 0; s < n_s; s++) {
    log_B(s,0) = log(B0*exp(omega_B(v_i(s),0)));
    B(s,0) = exp(log_B(s,0));
    resid_B(s,0) = B(s,0) - exp(log_B(s,0));
    for (int t = 1; t < (n_t); t++) {
      Type mean_pro = (exp(-m(s,t))*gI(t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
      log_B(s,t) = log(mean_pro);
      B(s,t) = exp(log_B(s,t));
      resid_B(s,t) = log_B(s,t) - log(mean_pro);
    }
    //Predict 1 year ahead, assume same spatial pattern as last year but propagate error
    Type mean_pro =(exp(-m(s,n_t))*gI(n_t-1)*(B(s, (n_t - 1)) - C(s,(n_t-1))) + exp(-m(s,n_t))*gR(n_t-1)*R(s,(n_t-1)))*exp(omega_B(v_i(s),n_t));
    log_B(s,n_t) = log(mean_pro);
    B(s,n_t) = exp(log_B(s,n_t));
    resid_B(s,n_t) = log_B(s,n_t) - log(mean_pro);
  }
  
  //Simulate biomass and commercial catch
  SIMULATE {
    // no_spat_B(0) = (B0*sum_area/1000)*exp(rnorm(Type(0.0),SigmaO_B));
    // no_spat_C(0) = no_spat_B(0)*0.1 * exp(rnorm(Type(0.0),Type(0.2)));
    // for (int t = 1; t < (n_t+1); t++){
    //   Type mean_no_spat = (exp(-no_spat_m(t))*gI(t-1)*(no_spat_B(t - 1) - no_spat_C(t-1)) + exp(-no_spat_m(t))*gR(t-1)*no_spat_R(t-1));
    //   no_spat_B(t) = mean_no_spat*exp(rnorm(Type(0.0),SigmaO_B));
    //   no_spat_C(t) = no_spat_B(t)*0.1*exp(rnorm(Type(0.0),Type(0.2)));
    // }
    vector <Type> counter_B(n_t); counter_B.setZero();
    for (int t = 0; t < n_t; t++){
      for (int s = 0; s < n_s; s++){
        if (t == 0){
          mean_pro_B(s,t) = B0 * exp(omega_B(v_i(s),t));
          // mean_pro_B(s,t) = B0 * exp(rnorm(Type(0.0),SigmaO_B));
          // B(s,t) = exp(log(mean_pro_B(s,t))-(sqr(sigma_tau)/2) + rnorm(Type(0.0),sigma_tau));
          // B(s,t) = no_spat_B(t)/sum_area*1000;
          B(s,t) = exp(log(mean_pro_B(s,t)));
          log_B(s,t) = log(B(s,t));
          counter_B(t) = counter_B(t) + (B(s,t)*area(s));
        }
        else {
          vector <Type> which_B(n_s); which_B.setZero();
          vector <Type> prev_B(n_s); prev_B = B.col(t-1);
          Type sum_prop; sum_prop = 0;
          Type num_B; num_B = 0;
          for (int i = 0; i < n_s; i++){
            if ((prev_B(i)*area(i)) >  (counter_B(t-1)/n_s)) {
              num_B += 1;
              sum_prop = sum_prop + (prev_B(i)*area(i));
            }
          }
          
          Type prop_B; prop_B = (B(s,t-1)*area(s))/sum_prop;
          // Type exploitation = Type(0.0);
          Type exploitation =  exp(log(((counter_B(t-1)/1000)*0.1)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
          if ( (B(s,t-1)*area(s)) > (counter_B(t-1)/n_s) ) C(s,t-1) = (exp(log((exploitation/area(s))*prop_B))+ rnorm(Type(0.0),Type(0.2)));
          // else C(s,t-1)=B(s,t-1)*0.02*exp(rnorm(Type(0.0),Type(0.5)));
          else C(s,t-1) = Type(0.0);
          // C(s,t-1) = exp(log((exploitation*prop_B)/area(s)) + rnorm(Type(0.0),Type(0.2)));
          // C(s,t-1) = B(s,t-1)*0.1 *exp(rnorm(Type(0.0),Type(0.2)));
          // C(s,t-1) = no_spat_C(t-1)/sum_area*1000;
          mean_pro_B(s,t) = (exp(-m(s,t))*gI(t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(t-1)*R(s,t-1))*exp(omega_B(v_i(s),t));
          // mean_pro_B(s,t) = (exp(-m(s,t))*gI(t-1)*(B(s, t - 1) - C(s,t-1)) + exp(-m(s,t))*gR(t-1)*R(s,t-1))*exp(rnorm(Type(0.0),SigmaO_B));
          // B(s,t) = no_spat_B(t)/sum_area*1000;
          B(s,t) = exp(log(mean_pro_B(s,t)));
          log_B(s,t) = log(B(s,t));
          counter_B(t) = counter_B(t) + (B(s,t)*area(s));
        }
      }
    }
    for (int s = 0; s < n_s; s++){
      vector <Type> which_B(n_s); which_B.setZero();
      vector <Type> prev_B(n_s); prev_B = B.col(n_t-1);
      Type sum_prop; sum_prop = 0;
      Type num_B; num_B = 0;
      for (int i = 0; i < n_s; i++){
        if ((prev_B(i)*area(i)) >  (counter_B(n_t-1)/n_s)) {
          num_B += 1;
          sum_prop = sum_prop + (prev_B(i)*area(s));
        }
      }
      Type prop_B; prop_B = (B(s,n_t-1)*area(s))/sum_prop;
      // Type exploitation = Type(0.0);
      Type exploitation = exp(log(((counter_B(n_t-1)/1000)*0.1)))*exp(rnorm(Type(0.0),Type(0.2)))*1000;
      if ( (B(s,n_t-1)*area(s)) > (counter_B(n_t-1)/n_s) ) C(s,n_t-1) = exp(log((exploitation/area(s))*prop_B) + rnorm(Type(0.0),Type(0.2)));
      // else C(s,n_t-1)=B(s,n_t-1)*0.02*exp(rnorm(Type(0.0),Type(0.2)));
      else C(s,n_t-1) = Type(0.0);
      // C(s,n_t-1) = exp(log((exploitation*prop_B)/area(s)) + rnorm(Type(0.0),Type(0.2)));
      // C(s,n_t-1) = B(s,n_t-1)*0.1 *exp(rnorm(Type(0.0),Type(0.2)));
      // C(s,n_t-1) = no_spat_C(n_t-1)/sum_area*1000;
      mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(n_t-1)*R(s,n_t-1))*exp(omega_B(v_i(s),n_t));
      // mean_pro_B(s,n_t) = (exp(-m(s,n_t))*gI(n_t-1)*(B(s, n_t - 1) - C(s,n_t-1)) + exp(-m(s,n_t))*gR(n_t-1)*R(s,n_t-1))*exp(rnorm(Type(0.0),SigmaO_B));
      // pred_B(s) = log(mean_pro_B(s,n_t)) -sqr(SigmaO_B)/2.0 + rnorm(Type(0.0),SigmaO_B);
      B(s,n_t) = mean_pro_B(s,n_t);
      // B(s,n_t) = no_spat_B(n_t)/sum_area*1000;
      // B(s,n_t) = exp(pred_B(s));
      log_B(s,n_t) = log(B(s,n_t));
      C(s,n_t) = 0;
    }
    // REPORT(pred_B);
    REPORT(B);
    REPORT(log_B);
    REPORT(C);
  }
  
  
  //Calculating predicted biomass and recruitment over area covered by each knots
  //Commercial biomass
  for (int s = 0; s < n_s; s++){
    for (int t = 0; t < (n_t+1); t++){
      areaB(s,t) = B(s,t) * area(s);
      log_areaB(s,t) = log(areaB(s,t));
    }
  }
  
  SIMULATE{
    for (int s = 0; s < n_s; s++){
      for (int t = 0; t < (n_t+1); t++){
        areaB(s,t) = B(s,t) * area(s);
        log_areaB(s,t) = log(areaB(s,t));
      }
    }
    REPORT(areaB);
    REPORT(log_areaB);
  }
  
  //Recruits
  for (int s = 0; s < n_s; s++){
    for (int t = 0; t < (n_t); t++){
      areaR(s,t) = R(s,t) * area(s);
      log_areaR(s,t) = log(areaR(s,t));
    }
  }
  
  SIMULATE{
    for (int s = 0; s < n_s; s++){
      for (int t = 0; t < (n_t); t++){
        areaR(s,t) = R(s,t) * area(s);
        log_areaR(s,t) = log(areaR(s,t));
      }
    }
    REPORT(areaR);
    REPORT(log_areaR);
  }
  
  //Calculate mean natural mortality, and total biomass and recruitment
  for (int t = 0; t < (n_t+1); t++){
    for (int s = 0; s < (n_s); s++){
      mean_m(t) = mean_m(t) + m(s,t);
    }
    mean_m(t) = mean_m(t) / n_s;
  }
  
  //Divide by 1000 to represent metric tonnes
  for (int t = 0; t < (n_t+1); t++){
    for (int s = 0; s < n_s; s++){
      totB(t) = totB(t) + (areaB(s,t)/1000);
    }
    log_totB(t) = log(totB(t));
  }
  
  for (int t = 0; t < (n_t); t++){
    for (int s = 0; s < n_s; s++){
      totR(t) = totR(t)+(areaR(s,t)/1000);
    }
    log_totR(t) = log(totR(t));
  }
  
  // Observation equations
  
  //Probability of capturing commercial biomass
  for (int t = 0; t <n_t; t++){
    nll_comp[5] -= dbinom_robust(pos_tows_I(t), n_tows(t), logit(p_I),true);
  }
  
  SIMULATE{
    for (int t = 0; t < n_t; t++){
      pos_tows_I(t) = rbinom(n_tows(t),p_I);
    }
    REPORT(pos_tows_I);
  }
  
  // Set prior on q
  if (prior_q == 1){
    for (int s = 0; s < n_s; s++){
      nll_comp[12] -= dbeta(qI(s),Type(10.0),Type(12.0),true);
    }
  }
  // if (prior_q == 1){
  //   nll_comp[12] -= dbeta(qI,Type(10.0),Type(12.0),true);
  // }
  
  // Commercial Index
  for (int i = 0; i < n_i; i++){
    if( !isNA(I(i) )) {
      Type mean_B = qI(s_i(i))*B(s_i(i),t_i(i))/p_I;
      nll_comp(6) -= dnorm(log(I(i)), log (mean_B)- sqr(sigma_epsilon)/Type(2.0), sigma_epsilon, true);
      resid_I(i) = log(I(i)) - log(mean_B);
    }
  }
  
  SIMULATE{
    for (int i = 0; i < n_i; i++){
      Type mean_I = qI(s_i(i))*B(s_i(i),t_i(i))/p_I;
      bern_I(i) = rbinom(Type(1.0),p_I);
      if (bern_I(i)>0) I(i) = exp(log(mean_I) - (sqr(sigma_epsilon)/2.0) + rnorm(Type(0.0),sigma_epsilon));
      else I(i) = NA_REAL;
    }
    REPORT(bern_I);
    REPORT(I);
  }
  
  //Probability of capturing recruits
  for (int t = 0; t <n_t; t++){
    nll_comp[7] -= dbinom_robust(pos_tows_IR(t), n_tows(t), logit(p_IR),true);
  }
  
  SIMULATE {
    for (int t = 0; t < n_t; t++){
      pos_tows_IR(t) = rbinom(n_tows(t),p_IR);
    }
    REPORT(pos_tows_IR);
  }
  
  // Recruit Index
  for (int i = 0; i < n_i; i++){
    if ( !isNA(IR(i) )) {
      Type mean_R = qR*R(s_i(i),t_i(i))/p_IR;
      nll_comp(8) -= dnorm(log(IR(i)), log(mean_R)-sqr(sigma_upsilon)/Type(2.0), sigma_upsilon, true);
    }
  }
  
  SIMULATE{
    for (int i = 0; i < n_i; i++){
      Type mean_IR = qR*R(s_i(i),t_i(i))/p_IR;
      bern_IR(i) = rbinom(Type(1.0),p_IR);
      if(bern_IR(i)>0) IR(i) = exp(log(mean_IR) - (sqr(sigma_upsilon)/2.0) + rnorm(Type(0.0),sigma_upsilon));
      else IR(i) = NA_REAL;
    }
    REPORT(bern_IR);
    REPORT(IR);
  }
  
  
  // Observations to natural mortality
  for (int i = 0; i < n_i; i++){
    if (!isNA(L(i))) {
      nll_comp(9) -= dbinom_robust(L(i), n_bin(i),logit(m(s_i(i),t_i(i))*S),true);
    }
  }
  
  SIMULATE {
    for (int i = 0; i < n_i; i++){
      L(i) = rbinom(n_bin(i), m(s_i(i),t_i(i))*S);
    }
    REPORT(L);
  }
  
  //Reporting
  ADREPORT(kappa_B);
  ADREPORT(tau_B);
  ADREPORT(kappa_R);
  ADREPORT(tau_R);
  ADREPORT(kappa_m);
  ADREPORT(tau_m);
  ADREPORT(sigma_epsilon);
  ADREPORT(sigma_upsilon);
  ADREPORT(S);
  ADREPORT(R0);
  ADREPORT(B0);
  ADREPORT(m0);
  ADREPORT(qI);
  ADREPORT(qR);
  ADREPORT(p_I);
  ADREPORT(p_IR);
  ADREPORT(log_H_input_B);
  ADREPORT(log_H_input_m);
  ADREPORT(H_B);
  ADREPORT(H_m);
  REPORT(SigmaO_B);
  REPORT(SigmaO_R);
  REPORT(SigmaO_m);
  REPORT(Range_B);
  REPORT(Range_R);
  REPORT(Range_m);
  ADREPORT(Range_B);
  ADREPORT(Range_R);
  ADREPORT(Range_m);
  REPORT(omega_R);
  // ADREPORT(omega_R);
  REPORT(omega_B);
  // ADREPORT(omega_B);
  REPORT(omega_m);
  // ADREPORT(omega_m);
  REPORT(log_B);
  REPORT(B);
  REPORT(areaB);
  // ADREPORT(areaB);
  REPORT(log_R);
  REPORT(R);
  REPORT(areaR);
  // ADREPORT(areaR);
  REPORT(log_m); 
  REPORT(m);
  ADREPORT(R);
  ADREPORT(m);
  ADREPORT(B);
  REPORT(totB); 
  ADREPORT(totB);
  REPORT(totR);
  ADREPORT(totR);
  REPORT(log_totB);
  ADREPORT(log_totB);
  REPORT(log_totR);
  ADREPORT(log_totR);
  REPORT(mean_m);
  ADREPORT(mean_m);
  ADREPORT(SigmaO_B);
  ADREPORT(SigmaO_R);
  ADREPORT(SigmaO_m);
  REPORT(resid_B);
  REPORT(resid_I);
  REPORT(nll_comp);
  
  Type nll = nll_comp.sum();
  return nll;
}
