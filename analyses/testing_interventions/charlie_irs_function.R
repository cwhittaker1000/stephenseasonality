
charlie_irs_function <- function(EIR_prev_reference_df, #Reference df
                                 years_run, #Number of years to run the simulation for
                                 year_irs_implemented, #The year irs coverage is input
                                 day_irs_implemented, #The day irs coverage is input, accepts values 1:365
                                 irs_half_life, #The irs half life
                                 IRS_interval, #How often irs is implemented in years, default is 1
                                 irs_cov, #irs coverage
                                 match_prev, #Prevalence we want to attempt to match to
                                 custom_seasonality # your custom seasonality vector, must be 365 values in length
                                 ){
  
 
  #This was built for using a LHC of parameter draws, so we're just going to use the closest values to the defaults
  EIR_density_reference_subset_1 <- EIR_prev_reference_df[which(EIR_prev_reference_df$itn_cov == 0 & EIR_prev_reference_df$irs_cov == 0 & EIR_prev_reference_df$delayMos == 10), ]
  EIR_density_reference_subset_2 <- EIR_density_reference_subset_1[which.min(abs(EIR_density_reference_subset_1$prev - match_prev)), ]

  #Run the model
  density_vec <- rep(EIR_density_reference_subset_2$mv0, 365 * years_run)
  
  #IRS parameters
  t_vector_irs <- c(-100, 2, (365 * year_irs_implemented) + day_irs_implemented) #Ignore the first 2 arguments, this has to start with a negative, then close to zero, then whatever you want after that
  irs_cov <- c(0, 0, irs_cov) #Additionally we start with 0 here

  #ITN parameters, not going to use but need to specify 9
  t_vector_itn <- t_vector_irs
  itn_cov <- c(0, 0, 0)
  
  #Resistance parameters for ITN's, just need to specify these
  itn_resistance <- subset(fread("//fi--didef3.dide.ic.ac.uk/Malaria/Arran/stephensi_east_africa/data/net_irs/pyrethroid_net_input.csv"), type == "Pyrethroid" & resistance_level == 1)
  
  #These
  model_run_fit <- create_r_model_epidemic(
    
    #Odin model to use 
    odin_model_path = "//fi--didef3.dide.ic.ac.uk/Malaria/Arran/stephensi_east_africa/odin/odin_model_seasonality_future_intervention.R",
    
    #Number of interventions, heterogeneity and age groups
    num_int = 4,
    het_brackets = 5,
    age = c(0, 1, 2, 5, 10, 20, 80),#c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80),
    
    #Initial prevalence and where it ends converted to density
    init_EIR = EIR_density_reference_subset_2$EIR,
    density_vec = density_vec,
    time_length = length(density_vec),
    
    #Initial nets, irs and treatment 
    init_ft = 0.5,
    itn_cov = 1e-25,      #These have to be a very small non-zero number if they are "zero"
    irs_cov = last(irs_cov),
    
    #Coverage of interventions
    ITN_IRS_on = min(c(t_vector_irs[which(t_vector_irs > 0)], t_vector_itn[which(t_vector_itn > 0)])),
    IRS_on = first(t_vector_irs[which(t_vector_irs > 0)]),
    
    irs_vector = irs_cov,
    itn_vector = itn_cov,            #Making these the same at 0 coverage
    larvicide_vector = itn_cov,      #Making these the same at 0 coverage
    
    #Timing of interventions
    t_vector_irs = t_vector_irs,
    t_vector_itn = t_vector_itn,            #Making these the same as they have 0 coverage
    t_vector_larvicide = t_vector_itn,      #Making these the same as they have 0 coverage
    
    #Seasonality
    admin_unit = NULL,
    country = NULL,
    custom_seasonality = custom_seasonality,     
    
    #Mosquito bionomic guesstimates - these are the defaults for stephensi, feel free to change
    mu0 = 0.12345679, 
    Q0 = 0.25, 
    chi = 0.5, 
    bites_Bed = 0.4776, 
    bites_Indoors = 0.52186,
    delayMos = 10,
    
    #IRS half-life and implementation parameters
    IRS_interval = IRS_interval, #How often IRS is applied in years, the default is every year
    irs_half_life = irs_half_life, #IRS halflife, default is 0.5
    
    #Resistance parameters - Ignore anything to do with ITN's only need to change the IRS parameters
    d_ITN0 = itn_resistance$ERG_d_ITN0,
    r_ITN0 = itn_resistance$ERG_r_ITN0,
    baseline_itn_half_life = itn_resistance$itn_half_life,
    baseline_d_ITN0 = itn_resistance$ERG_d_ITN0,
    baseline_r_ITN0 = itn_resistance$ERG_r_ITN0,
    itn_half_life = itn_resistance$itn_half_life,
    time_resistance = 365 * 3                       #Burn in time as it takes a while to reach eqt_vector_irs[which(t_vector_irs != 0)]uilibrium
    
  )
  
  #There is some post processing because the deterministic model doesnt like putting in coverage at different times, you can ignore this
  # Edits equilibrium condition with new coverage split. Using the default of 0.25 to all because its the middle ground
  mod <- model_run_fit
  mod <- edit_equilibrium_varying_nets_irs(mod,  as.numeric(c(0.25, 0.25, 0.25, 0.25)))
  mod$state$t_vector_irs <- as.numeric(mod$state$t_vector_irs)
  
  mod <- mod$generator(user = mod$state, use_dde = TRUE)
  mod_run <- mod$run(t = 1:length(density_vec))
  
  out <- mod$transform_variables(mod_run)
  model_ran <- as.data.frame(out)
  
  #Unpack results
  these_columns <- c("Incidence" = "Incidence",
                     "Prevalence" = "prev",
                     "Vector density" = "mv")
  
  unpack_take_these <- strsplit("Incidence;Prevalence;Vector density", ";")[[1]]
  extract_values <- unlist(model_ran[, these_columns[which(names(these_columns) %in% unpack_take_these)]])
  
  result_df <- data.frame( years_run, #Number of years to run the simulation for
                           year_irs_implemented, #The year irs coverage is input
                           day_irs_implemented, #The day irs coverage is input, accepts values 1:365
                           irs_half_life, #The irs half life
                           irs_cov, #irs coverage
                           prevalence_calibrated = match_prev,
                           t = model_ran$t,
                           type = rep(unpack_take_these, 
                                      each = length(model_ran$t)),
                           value = as.numeric(extract_values),
                           stringsAsFactors = FALSE)
  result_df
  
  
 
  
  
}