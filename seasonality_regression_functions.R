#This script contains the functions used for the seasonality regression analyses

###################################################
#Seasonality regression models
###################################################

#Model A: sinusoidal(t) + linear(t)
model_A_func = function(data, model_matrix, omega, on) {
  
  #get number of clinical per year combinations
  n_cys = length(unique(data[[on]]))
  
  #get first guess for clinic per year intercepts: use mean values
  start_intercepts = data %>%
    group_by(!!sym(on)) %>%
    summarize(y = mean(y)) %>%
    arrange(!!sym(on)) %T>%
    #check that order matches model matrix columns
    {stopifnot(all(.[[on]] == colnames(model_matrix)))} %>%
    pull(y)
  
  #get first guess for clinic per year slopes: set as 0
  start_slopes = rep(0, n_cys)
  
  #fit the model
  model = nls(y ~ amplitude * cos(omega * (month - phase)) + drop(model_matrix %*% slope) * month + drop(model_matrix %*% intercept),
              start = list(amplitude = 0.1, phase = 0, slope = start_slopes, intercept = start_intercepts),
              data = data) 
  
  return(model)
}

#Model B: MIC ~ sinusoidal(t) + linear(t) + age + sex
model_B_func = function(data, model_matrix, omega, on) {
  
  #get number of clinical per year combinations
  n_cys = length(unique(data[[on]]))
  
  #get first guess for clinic per year intercepts: use mean values
  start_intercepts = data %>%
    group_by(!!sym(on)) %>%
    summarize(y = mean(y)) %>%
    arrange(!!sym(on)) %T>%
    #check that order matches model matrix columns
    {stopifnot(all(.[[on]] == colnames(model_matrix)))} %>%
    pull(y)
  
  #get first guess for clinic per year slopes: set as 0
  start_slopes = rep(0, n_cys)
  
  #fit the model
  model = nls(y ~ amplitude * cos(omega * (month - phase)) + beta_age * age + beta_sex * is_male + drop(model_matrix %*% slope) * month + drop(model_matrix %*% intercept),
              start = list(amplitude = 0.1, phase = 0, beta_age = 0, beta_sex = 0, slope = start_slopes, intercept = start_intercepts),
              data = data) 
  
  return(model)
}

#Model C: MIC ~ sinusoidal(t) + linear(t) + age + sex + site_of_infection
model_C_func = function(data, model_matrix, omega, on) {
  
  #get number of clinical per year combinations
  n_cys = length(unique(data[[on]]))
  
  #get first guess for clinic per year intercepts: use mean values
  start_intercepts = data %>%
    group_by(!!sym(on)) %>%
    summarize(y = mean(y)) %>%
    arrange(!!sym(on)) %T>%
    #check that order matches model matrix columns
    {stopifnot(all(.[[on]] == colnames(model_matrix)))} %>%
    pull(y)
  
  #get first guess for clinic per year slopes: set as 0
  start_slopes = rep(0, n_cys)
  
  #fit the model
  model = nls(y ~ amplitude * cos(omega * (month - phase)) +
                beta_age * age + beta_sex * is_male + beta_SST * is_SST + beta_AB * is_AB + beta_BL * is_BL + beta_RT * is_RT +
                drop(model_matrix %*% slope) * month + drop(model_matrix %*% intercept),
              start = list(amplitude = 0.1, phase = 0, beta_age = 0, beta_sex = 0, beta_SST = 0, beta_AB = 0, beta_BL = 0, beta_RT = 0, slope = start_slopes, intercept = start_intercepts),
              data = data) 
  
  return(model)
}

###################################################
#Other functions 
###################################################

#Function to create a model matrix: where rows represent data rows and columns represent
#each year or clinic/year. The model matrix consists of 0s and 1s, where a 1 means that this 
#data row is from the corresponding year or clinic/year. The model matrix is used in the 
#nls regression to allow for a separate slope and intercept term to be estimated for each
#year or clinic/year.
make_model_matrix_func = function(data, on) {
  #get number of clinical per year combinations
  n_cys = length(unique(data[[on]]))
  
  #create model matrix - same number of rows as data rows, 0/1s, 1 if this row in clinic year
  model_matrix = model.matrix(formula(paste("~ 0 +", on)), data) %>%
    set_colnames(str_remove(colnames(.), on))
  
  #check that num rows in model matrix same as data
  stopifnot(dim(model_matrix) == c(nrow(data), n_cys))
  
  return(model_matrix)
  
}

#Function to make table of parameter values from regression model output
get_model_summary_func = function(model, model_matrix, on, num_covariates) {
  
  names = colnames(model_matrix)
  
  model_values1 = summary(model)$coefficients %>%
    as_tibble(rownames = "term") %>%
    set_colnames(c("term","estimate","std.error","statistic","p.value")) 
  
  model_values2 = confint.default(model) %>%
    set_colnames(c("ci.lower", "ci.upper")) %>%
    as_tibble(rownames = "term")
  
  model_values = left_join(model_values1, model_values2, by="term") %>%
    #add hospital/year to intercept/slope terms
    mutate(!!on := c(rep("", 2 + num_covariates), names, names)) %>%
    mutate(term = case_when(str_detect(term, "slope") ~ "slope",
                            str_detect(term, "intercept") ~ "intercept",
                            TRUE ~ term))
  
  return(model_values)
}

#Function to convert negative amplitudes to positive and phases to range from 0 months to model period. 
#A sinusoid with amplitude -A is the same as a sinusoid with amplitude A and phase + 0.5*period.  
#Similarly, a sinusoid with phase p is the same a sinuoid with phase p-period or p+period. 
convert_a_phases_func = function(a_estimate, a_ci.lower, a_ci.upper, phase_estimate, phase_ci.lower, phase_ci.upper, period) {
  
  #if the amplitude is negative, add (-2*a) to a, a_ci.lower, a_ci.upper; add 1/2*period to phase, phase_ci.lower, phase_ci.upper
  if (a_estimate < 0) {
    a_ci.lower = a_ci.lower + -2*a_estimate
    a_ci.upper = a_ci.upper + -2*a_estimate
    a_estimate = a_estimate + -2*a_estimate
    
    phase_estimate = phase_estimate + period/2
    phase_ci.lower = phase_ci.lower + period/2
    phase_ci.upper = phase_ci.upper + period/2
  }
  
  #while phase is not between 1-period, add or subtract period 
  while (!(phase_estimate >= 0 & phase_estimate <=period)) {
    if(phase_estimate < 0) {
      phase_estimate = phase_estimate + period
      phase_ci.lower = phase_ci.lower + period
      phase_ci.upper = phase_ci.upper + period
    }
    
    if(phase_estimate > period) {
      phase_estimate = phase_estimate - period
      phase_ci.lower = phase_ci.lower - period
      phase_ci.upper = phase_ci.upper - period
    }
  }
  
  dat = data.frame(amplitude_estimate = a_estimate, amplitude_ci.lower = a_ci.lower,
                   amplitude_ci.upper = a_ci.upper, phase_estimate = phase_estimate,
                   phase_ci.lower = phase_ci.lower, phase_ci.upper = phase_ci.upper)
  return(dat)
}


