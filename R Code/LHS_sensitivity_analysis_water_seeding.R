#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Team 3: Animal - Environment
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# libraries
library("deSolve")
library("ggplot2")
library("dplyr")
library("lhs") #library for Latin Hypercube Sampling

set.seed(42) 


# Model
AE_model <- function(time, state, parameters) {  
  with(as.list(c(state, parameters)), { 
    
    ########### Forces of infection:
    lambdaFishS <- beta_FishS*(FishCs)/(FishN+FishCs+FishCr)
    lambdaFishR <- beta_FishS*fc*(FishCr)/(FishN+FishCs+FishCr)
    
    lambdaFarmS <- beta_FarmS*(FarmCs)/(FarmN+FarmCs+FarmCr)
    lambdaFarmR <- beta_FarmS*fc*(FarmCr)/(FarmN+FarmCs+FarmCr)
    
    lambdaWildS <- beta_WildS*(WildCs)/(WildN+WildCs+WildCr)
    lambdaWildR <- beta_WildS*fc*(WildCr)/(WildN+WildCs+WildCr)
    
    lambdaPetS <- beta_PetS*(PetCs)/(PetN+PetCs+PetCr)
    lambdaPetR <- beta_PetS*fc*(PetCr)/(PetN+PetCs+PetCr)
    
    
    ########################################## Model equations: #############################################
    
    dFishN <- -(lambdaFishS + lambdaFishR)*FishN + (gammaFish+gammaFishABX*FishExp)*FishCs + gammaFish*FishCr
    
    dFishCs <- lambdaFishS*FishN - (gammaFish+gammaFishABX*FishExp)*FishCs + FishLoss*FishCr -
      FishCs*(FarmtoFish*FarmCr/(FarmCr+FarmCs+FarmN) + 
                PettoFish*PetCr/(PetCr+PetCs+PetN) +
                WildtoFish*WildCr/(WildCr+WildCs+WildN) +
                WatertoFish*(WaterCr/WaterCarrying) +
                SoiltoFish*(SoilCr/SoilCarrying))
    
    dFishCr <- lambdaFishR*FishN - gammaFish*FishCr - FishLoss*FishCr +
      FishCs*(FarmtoFish*FarmCr/(FarmCr+FarmCs+FarmN) + 
                PettoFish*PetCr/(PetCr+PetCs+PetN) +
                WildtoFish*WildCr/(WildCr+WildCs+WildN) +
                WatertoFish*(WaterCr/WaterCarrying) +
                SoiltoFish*(SoilCr/SoilCarrying))
    
    
    dFarmN <- -(lambdaFarmS + lambdaFarmR)*FarmN + (gammaFarm+gammaFarmABX*FarmExp)*FarmCs + gammaFarm*FarmCr
    
    dFarmCs <- lambdaFarmS*FarmN - (gammaFarm+gammaFarmABX*FarmExp)*FarmCs + FarmLoss*FarmCr -
      FarmCs*(PettoFarm*PetCr/(PetCr+PetCs+PetN) +
                WildtoFarm*WildCr/(WildCr+WildCs+WildN) +
                FishtoFarm*FishCr/(FishCr+FishCs+FishN) +
                WatertoFarm*(WaterCr/WaterCarrying) +
                SoiltoFarm*(SoilCr/SoilCarrying))
    
    dFarmCr <- lambdaFarmR*FarmN - gammaFarm*FarmCr - FarmLoss*FarmCr +
      FarmCs*(PettoFarm*PetCr/(PetCr+PetCs+PetN) +
                WildtoFarm*WildCr/(WildCr+WildCs+WildN) +
                FishtoFarm*FishCr/(FishCr+FishCs+FishN) +
                WatertoFarm*(WaterCr/WaterCarrying) +
                SoiltoFarm*(SoilCr/SoilCarrying))
    
    
    dWildN <- -(lambdaWildS + lambdaWildR)*WildN + (gammaWild+gammaWildABX*WildExp)*WildCs + gammaWild*WildCr
    
    dWildCs <- lambdaWildS*WildN - (gammaWild+gammaWildABX*WildExp)*WildCs + WildLoss*WildCr - 
      WildCs*(FarmtoWild*FarmCr/(FarmCr+FarmCs+FarmN) +
                PettoWild*PetCr/(PetCr+PetCs+PetN) +
                FishtoWild*FishCr/(FishCr+FishCs+FishN) +
                WatertoWild*(WaterCr/WaterCarrying) +
                SoiltoWild*(SoilCr/SoilCarrying))
    
    dWildCr <- lambdaWildR*WildN - gammaWild*WildCr - WildLoss*WildCr +
      WildCs*(FarmtoWild*FarmCr/(FarmCr+FarmCs+FarmN) +
                PettoWild*PetCr/(PetCr+PetCs+PetN) +
                FishtoWild*FishCr/(FishCr+FishCs+FishN) +
                WatertoWild*(WaterCr/WaterCarrying) +
                SoiltoWild*(SoilCr/SoilCarrying))
    
    
    dPetN <- -(lambdaPetS + lambdaPetR)*PetN + (gammaPet+gammaPetABX*PetExp)*PetCs + gammaPet*PetCr
    
    dPetCs <- lambdaPetS*PetN - (gammaPet+gammaPetABX*PetExp)*PetCs + PetLoss*PetCr -
      PetCs*(FarmtoPet*FarmCr/(FarmCr+FarmCs+FarmN) +
               WildtoPet*WildCr/(WildCr+WildCs+WildN) +
               FishtoPet*FishCr/(FishCr+FishCs+FishN) +
               WatertoPet*(WaterCr/WaterCarrying) +
               SoiltoPet*(SoilCr/SoilCarrying))
    
    dPetCr <- lambdaPetR*PetN - gammaPet*PetCr - PetLoss*PetCr +
      PetCs*(FarmtoPet*FarmCr/(FarmCr+FarmCs+FarmN) +
               WildtoPet*WildCr/(WildCr+WildCs+WildN) +
               FishtoPet*FishCr/(FishCr+FishCs+FishN) +
               WatertoPet*(WaterCr/WaterCarrying) +
               SoiltoPet*(SoilCr/SoilCarrying))
    
    
    dWaterCs <- WaterCs*WaterGrowthCs*(1-(WaterCs+WaterCr)/WaterCarrying) - WaterCs*WaterDecay - WaterCs*WaterHGT*(WaterCr/(WaterCs+WaterCr)) + WaterLoss*WaterCr - WaterCs*WaterAbx -
      WaterCs*(FarmtoWater*FarmCr/(FarmCr+FarmCs+FarmN) +
                 PettoWater*PetCr/(PetCr+PetCs+PetN) +
                 WildtoWater*WildCr/(WildCr+WildCs+WildN) +
                 FishtoWater*FishCr/(FishCr+FishCs+FishN) +
                 SoiltoWater*(SoilCr/SoilCarrying))
    
    dWaterCr <- WaterCr*WaterGrowthCr*fc*(1-(WaterCs+WaterCr)/WaterCarrying) - WaterCr*WaterDecay + WaterCs*WaterHGT*(WaterCr/(WaterCs+WaterCr)) - WaterLoss*WaterCr +
      WaterCs*(FarmtoWater*FarmCr/(FarmCr+FarmCs+FarmN) +
                 PettoWater*PetCr/(PetCr+PetCs+PetN) +
                 WildtoWater*WildCr/(WildCr+WildCs+WildN) +
                 FishtoWater*FishCr/(FishCr+FishCs+FishN) +
                 SoiltoWater*(SoilCr/SoilCarrying))
    
    
    dSoilCs <- SoilCs*SoilGrowthCs*(1-(SoilCs+SoilCr)/SoilCarrying) - SoilCs*SoilDecay - SoilCs*SoilHGT*(SoilCr/(SoilCs+SoilCr)) + SoilLoss*SoilCr - SoilCs*SoilAbx -
      SoilCs*(FarmtoSoil*FarmCr/(FarmCr+FarmCs+FarmN) +
                PettoSoil*PetCr/(PetCr+PetCs+PetN) +
                FishtoSoil*FishCr/(FishCr+FishCs+FishN) +
                WildtoSoil*WildCr/(WildCr+WildCs+WildN) +
                WatertoSoil*(WaterCr/WaterCarrying))
    
    dSoilCr <- SoilCr*SoilGrowthCr*fc*(1-(SoilCs+SoilCr)/SoilCarrying) - SoilCr*SoilDecay + SoilCs*SoilHGT*(SoilCr/(SoilCs+SoilCr)) - SoilLoss*SoilCr +
      SoilCs*(FarmtoSoil*FarmCr/(FarmCr+FarmCs+FarmN) +
                PettoSoil*PetCr/(PetCr+PetCs+PetN) +
                FishtoSoil*FishCr/(FishCr+FishCs+FishN) +
                WildtoSoil*WildCr/(WildCr+WildCs+WildN) +
                WatertoSoil*(WaterCr/WaterCarrying))
    
    
    return(list(c(dFishN, dFishCs, dFishCr, 
                  dFarmN, dFarmCs, dFarmCr, 
                  dWildN, dWildCs, dWildCr, 
                  dPetN, dPetCs, dPetCr, 
                  dWaterCs, dWaterCr, 
                  dSoilCs, dSoilCr))) 
  })
  
}



extract_outcomes <- function(output) {
  #How long until the ARG is above a threshold in 50% of compartments
  #(1 animal for farm animals, 1 fish for fish, 1 animal for wildlife, 1 animal for peridomestic, 100 cfu/mL or g for water and soil)
  min_time_Fish  <- ifelse(any(output$FishCr > 1), min(which(output$FishCr > 1)), NA)
  min_time_Farm  <- ifelse(any(output$FarmCr > 1), min(which(output$FarmCr > 1)), NA)
  min_time_Wild  <- ifelse(any(output$WildCr > 1), min(which(output$WildCr > 1)), NA)
  min_time_Pet   <- ifelse(any(output$PetCr > 1), min(which(output$PetCr > 1)), NA)
  min_time_Water <- ifelse(any(output$WaterCr > 100), min(which(output$WaterCr > 100)), NA)
  min_time_Soil  <- ifelse(any(output$SoilCr > 100), min(which(output$SoilCr > 100)), NA)
  
  # Time when 50% (3 out of 6) compartments cross threshold
  times_vec <- sort(c(min_time_Fish, min_time_Farm, min_time_Wild, min_time_Pet, min_time_Water, min_time_Soil))
  min_time_all <- if(length(times_vec) >= 3) times_vec[3] else NA
  
  # Prevalence at 1 week (Index 8 because seq starts at 0)
  idx1w <- 8
  idx1m <- 31
  
  prev_animals_1week <- sum(output$FishCr[idx1w], output$FarmCr[idx1w], output$WildCr[idx1w], output$PetCr[idx1w]) /
    sum(output[idx1w, c("FishN", "FishCs", "FishCr", "FarmN", "FarmCs", "FarmCr", "WildN", "WildCs", "WildCr", "PetN", "PetCs", "PetCr")])
  
  prev_envir_1week   <- sum(output$WaterCr[idx1w], output$SoilCr[idx1w]) /
    sum(output[idx1w, c("WaterCs", "WaterCr", "SoilCs", "SoilCr")])
  
  # Prevalence at 1 month
  prev_animals_1month <- sum(output$FishCr[idx1m], output$FarmCr[idx1m], output$WildCr[idx1m], output$PetCr[idx1m]) /
    sum(output[idx1m, c("FishN", "FishCs", "FishCr", "FarmN", "FarmCs", "FarmCr", "WildN", "WildCs", "WildCr", "PetN", "PetCs", "PetCr")])
  
  prev_envir_1month   <- sum(output$WaterCr[idx1m], output$SoilCr[idx1m]) /
    sum(output[idx1m, c("WaterCs", "WaterCr", "SoilCs", "SoilCr")])
  
  return(data.frame(min_time_all, prev_animals_1week, prev_envir_1week, prev_animals_1month, prev_envir_1month))
}



# kept 120 total (number of each animal group)
initial_state <- c(FishN = 96, FishCs = 24, FishCr = 0, # fish seeding
                   FarmN = 12, FarmCs = 108, FarmCr = 0, # farm seeding or 
                   WildN = 12, WildCs = 108, WildCr = 0, 
                   PetN = 12, PetCs = 108, PetCr = 0,
                   WaterCs = 10e5, WaterCr = 1, # water seeding 
                   SoilCs = 10e10, SoilCr = 0)   



sim_time <- seq(from = 0, to = 365, by = 1) 



# fixed parameters
fixed_params <- c(gammaFarm = 1/30,
               gammaFarmABX = 1/7,
               gammaPet = 1/30,
               gammaPetABX = 1/7,
               gammaWild = 1/30,
               gammaWildABX = 1/7,
               gammaFish = 1/10,
               gammaFishABX = 1/7,
               FarmLoss = 10e-6,
               PetLoss = 10e-6,
               WildLoss = 10e-6,
               FishLoss = 10e-6,
               
               WaterCarrying = 10e5, # carrying capacity of bacteria in water
               WaterHGT = 10e-6,     # rate of HGT in water
               WaterLoss = 0.05,     # rate of ARG loss in water
               
               SoilCarrying = 10e10, # carrying capacity of bacteria in soil
               SoilHGT = 10e-7,      # rate of HGT in soil
               SoilLoss = 0.05,       # rate of ARG loss in soil
               
               # contact rates that are 0 are kept constant
               FarmtoFish = 0, 
               PettoFish = 0,
               FishtoFarm = 0,
               FishtoPet = 0,
               FishtoSoil = 0,
               SoiltoFish = 0

)





# Latin Hypercube Sampling              
                
n_samples <- 1000
lhs_cube <- randomLHS(n_samples, 41) # nb of varying parameters 41 (listed below)


# varying parameters
params_set <- data.frame(
  beta_FarmS  = lhs_cube[,1] * (0.4 - 0.1) + 0.1, # within-compartment transmission rate
  beta_PetS = lhs_cube[,2] * (0.2 - 0.08) + 0.08,
  beta_WildS = lhs_cube[,3] * (0.1 - 0) + 0,
  beta_FishS = lhs_cube[,4] * (1.2 - 0.1) + 0.1,
  fc = lhs_cube[,5] * (0.99 - 0.85) + 0.85,  # relative fitness vary btwn 0.85 an 0.99
  FarmExp = lhs_cube[,6] * (0.99 - 0) + 0, # proportion of animals exposed to ABX
  PetExp = lhs_cube[,7] * (0.3 - 0) + 0,
  WildExp = lhs_cube[,8] * (0.05 - 0) + 0,
  FishExp = lhs_cube[,9] * (0.99 - 0) + 0,
  WaterGrowthCs = lhs_cube[,10] * (0.22 - 0.18) + 0.18,   # 0.2 #growth rate of S bacteria in water, VARY +/-10%
  WaterGrowthCr = lhs_cube[,11] * (0.22 - 0.18) + 0.18,   # 0.2 #growth rate of R bacteria in water, VARY +/-10%
  WaterDecay = lhs_cube[,12] * (0.022 - 0.018) + 0.018,   # 0.02 #decay rate of bacteria in water, VARY +/-10%
  WaterAbx = lhs_cube[,13] * (0.0011 - 0.0009) + 0.0009,  # 0.001 #rate of antibiotic exposure in water, VARY +/-10%
  SoilGrowthCs = lhs_cube[,14] * (0.11 - 0.09) + 0.09, # 0.1,   # growth rate of S bacteria in soil, VARY +/-10%
  SoilGrowthCr = lhs_cube[,15] * (0.11 - 0.09) + 0.09, # 0.1,   # growth rate of R bacteria in soil, VARY +/-10%
  SoilDecay = lhs_cube[,16] * (0.011 - 0.009) + 0.009, #0.01,      # decay rate of bacteria in soil, VARY +/-10%
  SoilAbx = lhs_cube[,17] * (0.0011 - 0.0009) + 0.0009, #0.001,      # rate of antibiotic exposure in soil, VARY +/-10%
  
  # Contact rates btwn ecological compartments, ### VARY all +/-10%
  FarmtoPet = lhs_cube[,18] * (0.022 - 0.018) + 0.018, #0.02, 
  FarmtoWild = lhs_cube[,19] * (0.011 - 0.009) + 0.009, #0.01,
  FarmtoWater = lhs_cube[,20] * (0.011 - 0.009) + 0.009,  #0.01,
  FarmtoSoil = lhs_cube[,21] * (0.011 - 0.009) + 0.009,  #0.01,
  #FarmtoFish = FishtoFarm, 0 zero
  
  PettoFarm = lhs_cube[,22] * (0.022 - 0.018) + 0.018,
  PettoWild = lhs_cube[,23] * (0.011 - 0.009) + 0.009,#0.01,
  PettoWater = lhs_cube[,24] * (0.011 - 0.009) + 0.009,#0.01,
  PettoSoil = lhs_cube[,25] * (0.011 - 0.009) + 0.009,#0.01,
  #PettoFish = FishtoPet, 0 zero
  
  WildtoFarm = lhs_cube[,26] * (0.011 - 0.009) + 0.009,
  WildtoPet = lhs_cube[,27] * (0.011 - 0.009) + 0.009,
  WildtoFish = lhs_cube[,28] * (0.0055 - 0.0045) + 0.0045, #0.005,
  WildtoWater = lhs_cube[,29] * (0.011 - 0.009) + 0.009, #0.01,
  WildtoSoil = lhs_cube[,30] * (0.011 - 0.009) + 0.009, #0.01,
  
  # FishtoFarm = FarmtoFish, 0 zero
  # FishtoPet = PettoFish, 0 zero
  FishtoWild = lhs_cube[,31] * (0.0055 - 0.0045) + 0.0045,
  FishtoWater = lhs_cube[,32] * (0.022 - 0.018) + 0.018, #0.02,
  #FishtoSoil = SoiltoFish, 0 zero
  
  WatertoFarm = lhs_cube[,33] * (0.011 - 0.009) + 0.009,
  WatertoPet = lhs_cube[,34] * (0.011 - 0.009) + 0.009,
  WatertoWild = lhs_cube[,35] * (0.011 - 0.009) + 0.009,
  WatertoFish = lhs_cube[,36] * (0.022 - 0.018) + 0.018,
  WatertoSoil = lhs_cube[,37] * (0.011 - 0.009) + 0.009,  #0.01,
  
  SoiltoFarm = lhs_cube[,38] * (0.011 - 0.009) + 0.009,
  SoiltoPet = lhs_cube[,39] * (0.011 - 0.009) + 0.009,
  SoiltoWild = lhs_cube[,40] * (0.011 - 0.009) + 0.009,
  # SoiltoFish = FishtoSoil, 0 zero
  SoiltoWater = lhs_cube[,41] * (0.011 - 0.009) + 0.009
  
)



results_list <- list()

cat("Starting LHS simulations...\n")


for (i in 1:n_samples) {
  # merge varying and fixed params
  current_params <- c(as.list(params_set[i, ]), as.list(fixed_params))
  
  # run ODE
  out <- as.data.frame(ode(y = initial_state, 
                           times = sim_time, 
                           func = AE_model, 
                           parms = current_params))
  
  # extract and store outcomes
  results_list[[i]] <- extract_outcomes(out)
  
  if(i %% 100 == 0) cat("Simulation", i, "of", n_samples, "complete...\n")
}

# df of results
lhs_final_results <- do.call(rbind, results_list)
summary(lhs_final_results)



# plot 5 outcomes above, 2 seeding scenarios Farm and Water 10 plots total
# Farm seeding and Water Seeding
lhs_final_results


library(sensitivity)

outcome_names <- c("min_time_all", "prev_animals_1week", "prev_envir_1week", 
                   "prev_animals_1month", "prev_envir_1month")

group_labels <- c(
  "min_time_all" = "Time to 50% Threshold",
  "prev_animals_1week" = "Animal Prevalence (1 Week)",
  "prev_envir_1week" = "Envir. Prevalence (1 Week)",
  "prev_animals_1month" = "Animal Prevalence (1 Month)",
  "prev_envir_1month" = "Envir. Prevalence (1 Month)"
)

all_prcc_data <- data.frame()


for (out_name in outcome_names) {
  # select the outcome column
  Y <- lhs_final_results[[out_name]]
  X <- params_set
  
  # remove any NA rows
  valid_rows <- !is.na(Y)
  
  # calculate PRCC
  # rank=TRUE calculates Partial Rank Correlation Coefficient (PRCC)
  pcc_res <- pcc(X[valid_rows, ], Y[valid_rows], rank = TRUE)
  
  # extract PRCC values
  df <- as.data.frame(pcc_res$PRCC)
  df$Parameter <- rownames(df)
  df$Group <- group_labels[out_name]
  colnames(df)[1] <- "PRCC"
  
  all_prcc_data <- rbind(all_prcc_data, df)
}


# add "Gap" rows so the plot respects spacing listed in param_order
gaps <- data.frame(
  PRCC = 0, 
  Parameter = c("gap1", "gap2", "gap3", "gap4", "gap5"),
  Group = rep(unique(all_prcc_data$Group), each = 5)
)
all_prcc_data <- rbind(all_prcc_data, gaps)






# named vector for the formatted labels
parameter_labels <- c(
  "gap1" = "", "gap2" = "", "gap3" = "", "gap4" = "", "gap5" = "",
  
  "beta_FarmS" = expression(beta[Farm]),
  "beta_PetS"  = expression(beta[Pet]),
  "beta_WildS" = expression(beta[Wild]),
  "beta_FishS" = expression(beta[Fish]),
  "fc"         = expression(f[r]),
  "FarmExp"    = "ABXFarm", # proportion of farm animals exposed to abx
  "PetExp"     = "ABXPet",
  "WildExp"    = "ABXWild",
  "FishExp"    = "ABXFish",
  "WaterGrowthCs" = "GSWater", # 0.2 #growth rate of S bacteria in water,
  "WaterGrowthCr" = "GRWater", # 0.2 #growth rate of R bacteria in water,
  "WaterDecay" = "DecayWater",   # 0.02 #decay rate of bacteria in water, VARY +/-10%
  "WaterAbx" = "ABXWater",  # 0.001 #rate of antibiotic exposure in water, VARY +/-10%
  "SoilGrowthCs" = "GSSoil", # 0.1,   # growth rate of S bacteria in soil, VARY +/-10%
  "SoilGrowthCr" = "GRSoil", # 0.1,   # growth rate of R bacteria in soil, VARY +/-10%
  "SoilDecay" = "DecaySoil", #0.01,      # decay rate of bacteria in soil, VARY +/-10%
  "SoilAbx" = "ABXSoil", #0.001,      # rate of antibiotic exposure in soil, VARY +/-10%
  
  # Contact rates btwn ecological compartments, ### VARY all +/-10%
  "FarmtoPet" = "FarmPet", #0.02, 
  "FarmtoWild" = "FarmWild", #0.01,
  "FarmtoWater" = "FarmWater",  #0.01,
  "FarmtoSoil" = "FarmSoil",  #0.01,
  #"FarmtoFish" = "cFarmFish", #0, using 1e-6 instead
  
  
  "PettoFarm" = "PetFarm",
  "PettoWild" = "PetWild",#0.01,
  "PettoWater" = "PetWater",#0.01,
  "PettoSoil" = "PetSoil",#0.01,
  #"PettoFish" = "cPetFish",#0, using 1e-6 instead
  
  
  "WildtoFarm" = "WildFarm",
  "WildtoPet" = "WildPet",
  "WildtoFish" = "WildFish", #0.005,
  "WildtoWater" = "WildWater", #0.01,
  "WildtoSoil" = "WildSoil", #0.01,
  
  # FishtoFarm = FarmtoFish,
  # FishtoPet = PettoFish,
  "FishtoWild" = "FishWild",
  "FishtoWater" = "FishWater", #0.02,
  #"FishtoSoil" = "cFishSoil",#0, using 1e-6 instead
  
  "WatertoFarm" = "WaterFarm",
  "WatertoPet" = "WaterPet",
  "WatertoWild" = "WaterWild",
  "WatertoFish" = "WaterFish",
  "WatertoSoil" = "WaterSoil",  #0.01,
  
  "SoiltoFarm" = "SoilFarm",
  "SoiltoPet" = "SoilPet",
  "SoiltoWild" = "SoilWild",
  # SoiltoFish = FishtoSoil,
  "SoiltoWater" = "SoilWater"
  
)








# order of params
param_order <- c(
  # relative fitness of R strains
  "fc",
  "gap1",
  # transmission rates
  "beta_FarmS", "beta_FishS", "beta_WildS", "beta_PetS", 
  "gap2",
  # ABX exposure
  "FarmExp", "FishExp", "WildExp", "PetExp", "WaterAbx", "SoilAbx",
  "gap3",
  # Growth of S
  "WaterGrowthCs", "SoilGrowthCs", 
  # Growth of R
  "WaterGrowthCr", "SoilGrowthCr", 
  "gap4",
  # Decay
  "WaterDecay", "SoilDecay", 
  "gap5",
  # contact rates
  "FarmtoPet", "FarmtoWild", "FarmtoWater", "FarmtoSoil",
  "PettoFarm", "PettoWild", "PettoWater", "PettoSoil",
  "WildtoFarm", "WildtoPet", "WildtoFish", "WildtoWater", "WildtoSoil",
  "FishtoWild", "FishtoWater", 
  "WatertoFarm", "WatertoPet", "WatertoWild", "WatertoFish", "WatertoSoil",
  "SoiltoFarm", "SoiltoPet", "SoiltoWild", "SoiltoWater"
)


all_prcc_data$Parameter <- factor(all_prcc_data$Parameter, levels = param_order) # rev
all_prcc_data$Group <- factor(all_prcc_data$Group, levels = group_labels)


#----------------------------------------- PLOT --------------------------------------------
all5_wide_row <- ggplot(all_prcc_data, aes(x = Parameter, y = PRCC, fill = PRCC)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = ifelse(abs(PRCC) > 0.1, round(PRCC, 2), "")), 
            hjust = ifelse(all_prcc_data$PRCC > 0, -0.2, 1.2), 
            fontface = "bold", size = 4) +
  facet_wrap(~Group, nrow = 1) + 
  coord_flip() +
  scale_fill_gradient2(low = "#d7191c", mid = "#ffffbf", high = "#2c7bb6") +
  ylim(-1.2, 1.2) + 
  scale_x_discrete(labels = parameter_labels, drop = FALSE) +
  labs(title = "Sensitivity Analysis Across Ecological Compartments (Water Seeding)",
       subtitle = "Comparing Parameter Impact On Resistance Prevalence, PRCC = Partial Rank Correlation Coefficient",
       x = "Model Parameter", 
       y = "Correlation Strength (PRCC Value)") +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14), 
    axis.text.x = element_text(size = 14), 
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    # Space between panels
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "gray95", color = NA)
  )

print(all5_wide_row)
ggsave("Sensitivity_Water_Seeding.png", plot = all5_wide_row, height = 8, width = 15)



# plot only top 10 parameters

sel_prcc_data <- all_prcc_data %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise(max_abs_prcc = max(abs(PRCC), na.rm = TRUE)) %>%
  dplyr::slice_max(order_by = max_abs_prcc, n = 10) %>%
  dplyr::pull(Parameter)

plot_df <- all_prcc_data %>%
  dplyr::filter(Parameter %in% sel_prcc_data)

# keep selected parameters in original param_order order
sel_levels <- param_order[param_order %in% sel_prcc_data]

# reset factor levels to only selected ones
plot_df$Parameter <- factor(plot_df$Parameter, levels = sel_levels)

# keep the same labels, but only for selected params
sel_labels <- parameter_labels[sel_levels]

sel5_wide_row <- ggplot(plot_df, aes(x = Parameter, y = PRCC, fill = PRCC)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = ifelse(abs(PRCC) > 0.1, round(PRCC, 2), "")), 
            hjust = ifelse(plot_df$PRCC > 0, -0.2, 1.2), 
            fontface = "bold", size = 4) +
  facet_wrap(~Group, nrow = 1) + 
  coord_flip() +
  scale_fill_gradient2(low = "#d7191c", mid = "#ffffbf", high = "#2c7bb6") +
  ylim(-1.2, 1.2) + 
  scale_x_discrete(labels = sel_labels) +
  labs(title = "Sensitivity Analysis Across Ecological Compartments (Water Seeding)",
       subtitle = "Comparing Parameter Impact On Resistance Prevalence, PRCC = Partial Rank Correlation Coefficient",
       x = "Model Parameter", 
       y = "Correlation Strength (PRCC Value)") +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14), 
    axis.text.x = element_text(size = 14), 
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    # Space between panels
    panel.spacing = unit(1, "lines"),
    strip.background = element_rect(fill = "gray95", color = NA)
  )

print(sel5_wide_row)
ggsave("Sensitivity_Water_Seeding_TOP10.png", plot = sel5_wide_row, height = 8, width = 15)
