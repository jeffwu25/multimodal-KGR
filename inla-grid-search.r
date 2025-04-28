### Load packages
library(dplyr) 
library(tidyverse)
library(lubridate)
library(stringr)
library(zoo)
library(ggplot2)
library(urbnmapr)
library(devtools)
library(readxl)
library(spdep)
library(sp)
library(huge)
library(INLA)
library(HMMpa)
library(invgamma)
library(brinla)
library(reshape2)
library(patchwork)
library(jsonlite)
library(geosphere)
library(urbnmapr)
library(RAQSAPI)
library(con2aqi)
library(pscl)
library(randtoolbox)
library(lhs)
library(scales)

### Load INLA workspace 
load("Multimodal-KGR-2-4.24.RData")

multimodal_df = multimodal_df %>% dplyr::filter(time <= 115)

### Excess mortality model
kgr_model1_deaths = function(dataset, rho_time_rbf = 1, rho_time_periodic = 1, sigma2_time = 1, link=1){
  
  #Calculating gram matrix K_time
  K_time = time_kernel(time_span = t,rho_rbf = rho_time_rbf, 
                       rho_periodic = rho_time_periodic, sigma2 = sigma2_time)
  
  #Heatmap of resulting K 
  K_time_heatmap = matrix_heatmap(K_time,title = "K_time heatmap")
  
  #Calculate trace norm of gram matrix
  K_time_weight = norm((1/t)*K_time,type = "F")
  
  #Calculate proposed kernel
  covGP = kronecker((K_time/t),(H1^2/n))
  
  #Need to ensure precision matrix is not computationally singular i.e det > 0
  covGP_jittered = desingularize(covGP,threshold = 1e-2,increment = 0.5)
  covGP = covGP_jittered[[1]]
  
  inv_covGP = solve(covGP)
  
  #Heatmap of resulting inv_covGP2 
  inv_covGP_heatmap = matrix_heatmap(inv_covGP,title = "")
  
  ###Fit INLA model 
  # kgr_formula1 = cumdelta_deaths ~ -1 + Intercept_1 + Intercept_2 + Intercept_3 + Intercept_4 + 
  #   Intercept_5 + Intercept_6 + Intercept_7 + Intercept_8 + Intercept_9 + f(id2,model = "generic0",Cmatrix = inv_covGP)
  
  kgr_formula1 = excess_deaths ~ -1 + month + Intercept_1 + Intercept_2 + Intercept_3 + Intercept_4 +
    Intercept_5 + Intercept_6 + Intercept_7 + Intercept_8 + Intercept_9 + f(id2,model = "generic0",Cmatrix = inv_covGP)
  
  model = inla(formula = kgr_formula1,family = "poisson",data = dataset, num.threads = 10,
               control.compute = list(dic=TRUE,waic=TRUE,
                                      return.marginals.predictor=TRUE),
               control.inla = list(strategy = "simplified.laplace"),
               control.predictor = list(compute = TRUE, link = link))
  
  ###Extract relevant information and store in the list
  # model_summary <- model$summary.fixed
  # bri_hyperpar_summary <- bri.hyperpar.summary(model)
  model_DIC <- model$dic$dic
  model_WAIC <- model$waic$waic
  # preds_model <- model$summary.fitted.values
  # preds_model <- cbind(dataset$region, dataset$date, preds_model)
  # colnames(preds_model) <- c("region", "date", "mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
  # marginal_fvs <- model$marginals.fitted.values
  # 
  # #Exponentiating parameter to get better interpretation of estimates 
  # multeff <- exp(model$summary.fixed$mean)
  # names(multeff) <- model$names.fixed
  # 
  # #Plot of each parameters' posterior density 
  # mf <- melt(model$marginals.fixed)
  # cf <- spread(mf,Var2,value)
  # names(cf)[2] <- 'parameter'
  # param_plot = ggplot(cf,aes(x=x,y=y)) + geom_line()+facet_wrap(~ parameter, 
  #                                                               scales="free") + geom_vline(xintercept=0) + ylab("density")
  # 
  # #Plot of precision of random effect (main hyperparameter of interest)
  # sden <- data.frame(bri.hyper.sd(model$marginals.hyperpar[[1]]))
  # hyperparam_plot = ggplot(sden,aes(x,y)) + geom_line() + ylab("density") + 
  #   xlab("linear predictor")
  
  #Store the results in the list
  kgr_model1_results = list(
    # K_time_heatmap = K_time_heatmap,
    # K_time_weight = K_time_weight/(K_time_weight + gfilter_weight1),
    # gfilter_weight = gfilter_weight1/(K_time_weight + gfilter_weight1),
    # covmatrix = covGP,
    # prec = inv_covGP,
    # num_jitters = covGP_jittered[[2]],
    # prec_heatmap = inv_covGP_heatmap,
    # model_summary = model_summary,
    # bri_hyperpar_summary = bri_hyperpar_summary,
    # exp_effects = multeff,
    # param_plot = param_plot,
    # hyperparam_plot = hyperparam_plot,
    model_DIC = model_DIC,
    model_WAIC = model_WAIC
    # fitted_values = preds_model,
    # marg_fitted_values = marginal_fvs
  )
  
  return(kgr_model1_results)
}

### Define grid for 3 parameters in K_time
hyper_grid = randomLHS(100,3)
colnames(hyper_grid) = c("rho_time_rbf","rho_time_periodic","sigma2_time")
hyper_grid[,1:2] = hyper_grid[,1:2]*0.1
hyper_grid[,3] = hyper_grid[,3]*2

### Perform the grid search

apply_function <- function(row, index) {
  rho_time_rbf <- row[1]
  rho_time_periodic <- row[2]
  sigma2_time <- row[3]
  
  result <- tryCatch(
    {
      kgr_model1_deaths(dataset = multimodal_df, rho_time_rbf = rho_time_rbf,
                 rho_time_periodic = rho_time_periodic, sigma2_time = sigma2_time)
    },
    error = function(e) {
      # Handle the error by printing a message and returning NULL
      cat("Error in model", index, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  # If result is not NULL, process the result
  if (!is.null(result)) {
    
    kgr_model1_deaths_running_results <<- c(kgr_model1_deaths_running_results, result$model_WAIC)
    # Print and save WAIC every 10th model
    if (index %% 10 == 0) {
      cat("Model", index, "WAIC:", result$model_WAIC, "\n")
      saveRDS(kgr_model1_deaths_running_results, paste0("kgr_model1_deaths_running_results_WAIC.rds"))
    }
  }
  
  # Return the result (which could be NULL if there was an error)
  return(result)
}

kgr_model1_deaths_running_results <- vector()  # Initialize an empty vector to store WAIC results

kgr_model1_deaths_results_list <- lapply(seq_len(nrow(hyper_grid)), function(i) {
  apply_function(hyper_grid[i, ], i)
})

# Remove NULL results from the list if needed
kgr_model1_deaths_results_list <- Filter(Negate(is.null), kgr_model1_deaths_results_list)

#Extracting WAIC values 
kgr_model1_deaths_results_WAIC = c()

for (i in 1:length(kgr_model1_deaths_results_list)){
  kgr_model1_deaths_results_WAIC[i] = pred_data = kgr_model1_deaths_results_list[[i]]$model_WAIC
}


saveRDS(kgr_model1_deaths_results_WAIC, paste0("kgr_model1_deaths_results_WAIC.rds"))


hist(kgr_model1_deaths_results_WAIC)
top5 = head(sort(kgr_model1_deaths_results_WAIC))
top5
top5_idx = which(kgr_model1_deaths_results_WAIC <= top5[5])

kgr_model1_deaths_results_WAIC = cbind(hyper_grid,kgr_model1_deaths_results_WAIC)
colnames(kgr_model1_deaths_results_WAIC) = c(colnames(hyper_grid),"WAIC")
kgr_model1_deaths_results_WAIC[top5_idx,]



### Avg max temperature
kgr_model1_temp = function(dataset, rho_time_rbf = 1, rho_time_periodic = 1, sigma2_time = 1, link=1){
  
  #Calculating gram matrix K_time
  K_time = time_kernel(time_span = t,rho_rbf = rho_time_rbf, 
                       rho_periodic = rho_time_periodic, sigma2 = sigma2_time)
  
  #Heatmap of resulting K 
  K_time_heatmap = matrix_heatmap(K_time,title = "K_time heatmap")
  
  #Calculate trace norm of gram matrix
  K_time_weight = norm((1/t)*K_time,type = "F")
  
  #Calculate proposed kernel
  covGP = kronecker((K_time/t),(H2^2/n))
  
  #Need to ensure precision matrix is not computationally singular i.e det > 0
  covGP_jittered = desingularize(covGP,threshold = 1e-2,increment = 0.5)
  covGP = covGP_jittered[[1]]
  
  inv_covGP = solve(covGP)
  
  #Heatmap of resulting inv_covGP2 
  inv_covGP_heatmap = matrix_heatmap(inv_covGP,title = "")
  
  ###Fit INLA model 
  # kgr_formula1 = cumpctdelta_temp ~ -1 + Intercept_1 + Intercept_2 + Intercept_3 + Intercept_4 + 
  #   Intercept_5 + Intercept_6 + Intercept_7 + Intercept_8 + Intercept_9 + f(id2,model = "generic0",Cmatrix = inv_covGP)
  
  kgr_formula1 = avg_max_temp ~ -1 + month + Intercept_1 + Intercept_2 + Intercept_3 + Intercept_4 +
    Intercept_5 + Intercept_6 + Intercept_7 + Intercept_8 + Intercept_9 + f(id2,model = "generic0",Cmatrix = inv_covGP)
  
  model = inla(formula = kgr_formula1,family = "gaussian",data = dataset, num.threads = 10,
               control.compute = list(dic=TRUE,waic=TRUE,
                                      return.marginals.predictor=TRUE),
               control.inla = list(strategy = "simplified.laplace"),
               control.predictor = list(compute = TRUE))
  
  ###Extract relevant information and store in the list
  # model_summary <- model$summary.fixed
  # bri_hyperpar_summary <- bri.hyperpar.summary(model)
  model_DIC <- model$dic$dic
  model_WAIC <- model$waic$waic
  # preds_model <- model$summary.fitted.values
  # preds_model <- cbind(dataset$region, dataset$date, preds_model)
  # colnames(preds_model) <- c("region", "date", "mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
  # marginal_fvs <- model$marginals.fitted.values
  # 
  # #Exponentiating parameter to get better interpretation of estimates 
  # multeff <- exp(model$summary.fixed$mean)
  # names(multeff) <- model$names.fixed
  # 
  # #Plot of each parameters' posterior density 
  # mf <- melt(model$marginals.fixed)
  # cf <- spread(mf,Var2,value)
  # names(cf)[2] <- 'parameter'
  # param_plot = ggplot(cf,aes(x=x,y=y)) + geom_line()+facet_wrap(~ parameter, 
  #                                                               scales="free") + geom_vline(xintercept=0) + ylab("density")
  # 
  # #Plot of precision of random effect (main hyperparameter of interest)
  # sden <- data.frame(bri.hyper.sd(model$marginals.hyperpar[[1]]))
  # hyperparam_plot = ggplot(sden,aes(x,y)) + geom_line() + ylab("density") + 
  #   xlab("linear predictor")
  
  #Store the results in the list
  kgr_model1_results = list(
    # K_time_heatmap = K_time_heatmap,
    # K_time_weight = K_time_weight/(K_time_weight + gfilter_weight2),
    # gfilter_weight = gfilter_weight2/(K_time_weight + gfilter_weight2),
    # covmatrix = covGP,
    # prec = inv_covGP,
    # num_jitters = covGP_jittered[[2]],
    # prec_heatmap = inv_covGP_heatmap,
    # model_summary = model_summary,
    # bri_hyperpar_summary = bri_hyperpar_summary,
    # exp_effects = multeff,
    # param_plot = param_plot,
    # hyperparam_plot = hyperparam_plot,
    model_DIC = model_DIC,
    model_WAIC = model_WAIC
    # fitted_values = preds_model,
    # marg_fitted_values = marginal_fvs
  )
  
  return(kgr_model1_results)
}

### Define grid for 3 parameters in K_time
hyper_grid = randomLHS(100,3)
colnames(hyper_grid) = c("rho_time_rbf","rho_time_periodic","sigma2_time")
hyper_grid[,1:2] = hyper_grid[,1:2]*100 
hyper_grid[,3] = hyper_grid[,3]*100 


### Perform the grid search
apply_function <- function(row, index) {
  rho_time_rbf <- row[1]
  rho_time_periodic <- row[2]
  sigma2_time <- row[3]
  
  result <- tryCatch(
    {
      kgr_model1_temp(dataset = multimodal_df, rho_time_rbf = rho_time_rbf,
                        rho_time_periodic = rho_time_periodic, sigma2_time = sigma2_time)
    },
    error = function(e) {
      # Handle the error by printing a message and returning NULL
      cat("Error in model", index, ":", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  # If result is not NULL, process the result
  if (!is.null(result)) {
    
    kgr_model1_temp_running_results <<- c(kgr_model1_temp_running_results, result$model_WAIC)
    # Print and save WAIC every 10th model
    if (index %% 10 == 0) {
      cat("Model", index, "WAIC:", result$model_WAIC, "\n")
      saveRDS(kgr_model1_temp_running_results, paste0("kgr_model1_temp_running_results_WAIC.rds"))
    }
  }
  
  # Return the result (which could be NULL if there was an error)
  return(result)
}

kgr_model1_temp_running_results <- vector()  # Initialize an empty vector to store WAIC results

kgr_model1_temp_results_list <- lapply(seq_len(nrow(hyper_grid)), function(i) {
  apply_function(hyper_grid[i, ], i)
})

# Remove NULL results from the list if needed
kgr_model1_temp_results_list <- Filter(Negate(is.null), kgr_model1_temp_results_list)

#Extracting WAIC values 
kgr_model1_temp_results_WAIC = c()

for (i in 1:length(kgr_model1_deaths_results_list)){
  kgr_model1_temp_results_WAIC[i] = pred_data = kgr_model1_temp_results_list[[i]]$model_WAIC
}


saveRDS(kgr_model1_temp_results_WAIC, paste0("kgr_model1_temp_results_WAIC.rds"))


hist(kgr_model1_temp_results_WAIC)
top5 = head(sort(kgr_model1_temp_results_WAIC))
top5
top5_idx = which(kgr_model1_temp_results_WAIC <= top5[5])

kgr_model1_temp_results_WAIC = cbind(hyper_grid,kgr_model1_temp_results_WAIC)
colnames(kgr_model1_temp_results_WAIC) = c(colnames(hyper_grid),"WAIC")
kgr_model1_temp_results_WAIC[top5_idx,]