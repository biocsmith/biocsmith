### glmmTMB functions
### Make note:

### I designed these functions to be used with a melted dataframe with rows representing samples 
### and columns representing taxa, counts of those taxa, 
### and relevant metadata to explore how counts of those taxa vary in abundance across treatments or groups.

##### `fml` and `fml_list` are written in R as:

# fml1 <- as.formula(paste("count", "~", paste(c("Age","offset(log(Library.Size))"), collapse = "+")))
# fml2 <- as.formula(paste("count", "~", paste(c("Age","offset(log(Library.Size))", "(1|Clutch)"), collapse = "+")))
# fml_list <- list(fml1, fml2)

##### with formulas ordered in the list by increasing complexity. 
##### Currently, the code has only been tested for with formulas with one response and one conditional variable.

library(tidyverse)
library(glmmTMB) ### GLMM
library(DHARMa) ### Zero inflation tests
library(MuMIn) ### Model Selection

countZeros <- function(x) sum(x == 0)

### This function runs glmmTMB with given formula `fml` and `data`, with no zero-inflation parameter. 
### It then uses `simulateResiduals` and `testGeneric` with the user-written function `countZeros` 
### to test for zero-inflation in the model fit.
### If p-value < 0.05, it reruns glmmTMB with the zero-inflation parameter. 
### In both scenarios, it returns a glmmTMB fit object.

glmm.fun <- function(fml, data) {
  fit <- glmmTMB(fml, zi = ~0, family = nbinom2, data = data)
  sim <- simulateResiduals(fit)
  if(testGeneric(sim, countZeros, plot = FALSE)$p.value > 0.05){
    return(fit)
  }else{
    fit <- glmmTMB(fml, zi = ~., family = nbinom2, data = data)
    return(fit)
  }
}

### This function runs the `glmm.fun` function above. 
### It takes a formula list `fml_list`, `data` and max delta AIC `delta.max` as inputs.
### It then runs `model.sel` and `get.models` to from DHARMa to select the best (below `delta.max`) and simplest models, 
### in this case, "fit1" or the firts forumal in the formula list.

model.fun <- function(fml_list, data, delta.max) {
  models <- as.list(1:length(fml_list))
  names(models) <- paste0("fit", 1:length(fml_list))
  for (f in 1:length(fml_list)) {
    models[[f]] <- glmm.fun(fml_list[[f]], data)
  }
  ms <- get.models(model.sel(models), subset = delta <= delta.max)
  if(length(ms) == 1){
    return(ms)
  }else{
    ms <- ms[names(ms) %in% "fit1" == TRUE]
    return(ms)
  }
  return(ms)
}

### After running `model.fun` the results are stored in a list of glmmTMB objects, where the names are the taxa.
### The following function summarizes each object in the list using the `summary` function, 
### extracts relevant coefficients, model info, and p-values and stores the information in a dataframe.
### The dataframe is now usable for making plots in ggplot.

convert_list_to_df <- function(results_list) {
  res_list <- list()
  for(p in names(results_list)){
    tryCatch({
      res <- summary(model_list[[p]][[1]])
      res_model <- as.data.frame(res$coefficients$cond) %>% rownames_to_column("Response")
      res_model$model <- "conditional"
      if(!is.null(res$coefficients$zi)){
        res_model_zi <- as.data.frame(res$coefficients$zi) %>% rownames_to_column("Response")
        res_model_zi$model = "zero-inflated"
        res_model = rbind(res_model, res_model_zi)
      }
      res_model$TaxaGroup <- p
      res_model$zi <- paste(res[["call"]][["ziformula"]][2])
      res_model$formula <- paste(as.formula(res$call)[3])
      res_list[[p]] <- list(res_model)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  res.df <- bind_rows(res_list)
  return(res.df)
}
