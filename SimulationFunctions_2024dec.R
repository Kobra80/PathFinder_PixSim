
### load libraries
library(data.table)
library(fst)
library("torch")


### load regeneration table that is used in growth projections where needed
load("~/PathFinder/Reg_Data.RData")


####### Simulation functions, models and parameters
# Create a list of models and parameters

ModelsAndParameters <- list(
  value = list(
    Model.1.biomass_I = "exp(-0.5 * ((biomass1 - a1) / a2)^2) * (a3 * npp)",  # Fixed expression
    Params.1.biomass_I = c(a1 = 356.6489, a2 = 296.2323, a3 = 0.004454764),
    
    Model.2.biomass_I = "exp(-0.5 * ((biomass1 - a1) / a2)^2) * (a3 * npp)",  # Fixed expression
    Params.2.biomass_I = c(a1 = 273.7304, a2 = 163.1729, a3 = 0.003230687),
    
    Model.3.biomass_I = "exp(-0.5 * ((biomass1 - a1) / a2)^2) * (a3 * npp)",  # Fixed expression
    Params.3.biomass_I = c(a1 = 176.0427, a2 = 115.0772, a3 = 0.002893258)
  ),
  visible = FALSE
)


######################################
## Growth models
GrowthModels <- function (Data, ModelsAndParameters, nSpecies) {
  
  if (!inherits(ModelsAndParameters, "list")) {
    stop("ModelsAndParameters object must be a list")
  }
  for (i in seq_along(ModelsAndParameters)) {
    assign(names(ModelsAndParameters)[i], ModelsAndParameters[[i]])
  }
  invisible(lapply(nSpecies, function(specie) {
    
    biomass_I.P <- paste0("Params.", specie, ".biomass_I")
    biomass_I.M <- paste0("Model.", specie, ".biomass_I")
    
    Data[code == 1 & Species == specie, `:=`(biomass_I, {
      for (i in seq_along(ModelsAndParameters[[biomass_I.P]])) {
        assign(names(ModelsAndParameters[[biomass_I.P]])[i], ModelsAndParameters[[biomass_I.P]][i])
      }
      Model.biomass_I <- parse(text = ModelsAndParameters[[biomass_I.M]])
      eval(Model.biomass_I)
    })]
    
    Data[code == 1 & Species == specie, `:=`(biomass2, biomass1 + biomass_I)]
    Data[, biomass_I:=NULL]
    
  }))
  Data[npp <=0 | is.na(npp) | biomass2 < 0, `:=`(code, 2)]
}

#######################################
# load the deep learning model
model <- torch_load( "~/PathFinder/Final_model_SP123.pth")


Deep_learning_model <- function(Model, Data) {
  
  #  mean and standard deviation vectors 
  x_mean <- c(biomass1 = 88.49560,  biomass2 = 98.73302, qmd1 = 16.02152, hl1 = 12.59579,
              ba1 = 20.32857,  vol1 = 135.58661,  npp = 4014.72350)
  
  x_sd <- c(biomass1 = 62.087784, biomass2 = 68.114217, qmd1 = 5.316511, hl1 = 3.843544,
            ba1 = 11.242124, vol1 = 105.364240,  npp = 1049.717469)
  
  
  y_mean <- c(qmd2 = 16.57135, hl2 = 13.24812, ba2 = 22.05845, vol2 = 153.46528)
  y_sd <- c(qmd2 = 5.345567, hl2 = 3.982301, ba2 = 11.909186, vol2 = 117.235923)
  
  
  
  # Define the input variables and species columns
  x_var <- c("biomass1", "biomass2", "qmd1", "hl1", "ba1", "vol1", "npp")
  sp_cols <- c("sp1", "sp2", "sp3")
  
  # Standardization function (for temporary use)
  standardize <- function(x, mean, sd) {
    (x - mean) / sd
  }
  
  
  # Destandardization function
  destandardize <- function(x, mean, sd) {
    return ((x * sd) + mean)
  }
  
  # Create a copy of Data to avoid modifying the original data
  x <- Data[, c(x_var, sp_cols), with = FALSE]  # Only select x_var and sp_cols
  
  # Standardize only the x_var columns (biomass1, qmd1, hl1, ba1, vol1, npp)
  for (var in x_var) {
    x[[var]] <- standardize(x[[var]], x_mean[var], x_sd[var])
  }
  
  # Prepare the input tensor (including species columns) from the standardized data
  x_tensor <- torch_tensor(as.matrix(x))
  
  # Set the model to evaluation mode
  Model$eval()
  
  # Make predictions (without calculating gradients)
  with_no_grad({
    y_pred_tensor <- Model(x_tensor)  # Pass the input tensor through the model
  })
  
  # Convert the tensor predictions to an array
  y_pred_array <- as.array(y_pred_tensor)
  
  
  # Destandardize the predicted values and add them to the original Data
  Data[, qmd2 := destandardize(y_pred_array[, 1], y_mean["qmd2"], y_sd["qmd2"])]
  Data[, hl2 := destandardize(y_pred_array[, 2], y_mean["hl2"], y_sd["hl2"])]
  Data[, ba2 := destandardize(y_pred_array[, 3], y_mean["ba2"], y_sd["ba2"])]
  Data[, vol2 := destandardize(y_pred_array[, 4], y_mean["vol2"], y_sd["vol2"])]
  
  
  # Classify and assign the 'code' column based on prediction values
  Data[hl2 < 0 | qmd2 < 0 | ba2 < 0 | vol2 < 0 | biomass2 < 0, `:=`(code, 2)]  # Invalid values
  Data[hl2 > 0 | qmd2 > 0 | ba2 > 0 | vol2 > 0 | biomass2 > 0, `:=`(code, 1)]  # Valid predictions
  Data[hl2 == 0 | qmd2 == 0 | ba2 == 0 | vol2 == 0 | biomass2 == 0, `:=`(code, 3)]  # Zero values
  
  # Return the updated Data table (DA) with predictions and actual values
  return(Data)
}


#####################################
## Regeneration function

RegFunction <- function (Data, Reg_Data, XX) {
  Data[code==3 & qmd1 == 0 & vol1 == 0, `:=`(Code, 1)]
  Reg_Data[, `:=`(Code, 1)]
  Data[Reg_Data, on = .(npp_class, Species, Code), `:=`(Interval.reg = i.Interval)]
  Data[, `:=`(Code, NULL)]
  Data[ !is.na(Interval.reg) & Interval == Interval.reg, `:=`(Code,1)]
  Data[Reg_Data, on = .(npp_class, Species, Code), `:=`(biomass1 = i.biomass1, qmd1 = i.qmd1, hl1 = i.hl1, ba1 = i.ba1, vol1=i.vol1, code=1 )]
  Reg_Data[, `:=`(Code, NULL)]
  Data[, `:=`(Code, NULL)]
}

########################################
#####PostRegeneration function

PostRegFunction <- function (Data) {
  Data[code==3 & Interval < Interval.reg, `:=`(biomass2 = NA, qmd2 = NA, hl2 = NA, ba2 = NA, vol2 = NA)]
}


########################################

## A local folder where simulation results should be written.
Fold <- tempfile()
dir.create(Fold)

#####################################################
## Simulation function

PixSim <- function (Data, Np, nSpecies, functions, WriteOut = FALSE, LocalFldr = NULL, ...) 
{
  Ellipsis <- list(...)
  
  if (!data.table::is.data.table(Data)) {
    stop("Data must be a data.table object")
  }
  
  if (!all(unique(Data$Species) %in% nSpecies)) {
    stop("Species in Data and nSpecies must be the same")
  }
  
  if (!("code" %in% colnames(Data))) {
    Data[, code := ifelse(biomass1== 0 & qmd1 == 0 & hl1 == 0 & ba1 == 0 & vol1 == 0, 3, 
                          ifelse(hl1 < 0 | qmd1 < 0 | ba1 < 0 | vol1 < 0 | npp <=0 | is.na(npp), 2, 1))]
  }    
  
  if (WriteOut) {
    fst::write_fst(Data, file.path(LocalFldr, "DataPred_000.fst"), compress = 100)
    gc()
  }
  
  for (XX in 1:Np) {
    Data[, Interval := XX]
    
    if ("RegFunction" %in% names(functions)) {
      functions$RegFunction(Data = Data, Reg_Data = Ellipsis$Reg_Data, XX = XX)
      # print(paste("RegFunction results for XX =", XX))
      #print(head(Data))
    }
    
    functions$GrowthModels(Data = Data, ModelsAndParameters = Ellipsis$ModelsAndParameters, nSpecies = nSpecies)
    #print(paste("GrowthModels results for XX =", XX))
    #print(head(Data))
    
    functions$Deep_learning_model(Model = model, Data = Data)
    
    
    if ("PostRegFunction" %in% names(functions)) {
      functions$PostRegFunction(Data = Data)
      # print(paste("PostRegFunction results for XX =", XX))
      # print(head(Data))
    }
    
    if (WriteOut) {
      ProjTosave <- Data[, .(biomass2, qmd2, hl2, ba2, vol2)]
      data.table::setnames(ProjTosave, old = c("biomass2", "qmd2", "hl2", "ba2", "vol2"), new = c("biomass1", "qmd1", "hl1", "ba1", "vol1"))
      ProjTosave <- ProjTosave[, round(.SD, 3), .SDcols = 1:5]
      nameTosave <- file.path(LocalFldr, paste0("DataPred_", sprintf("%04d", XX), ".fst"))
      fst::write_fst(ProjTosave, nameTosave, compress = 100)
      rm(ProjTosave, nameTosave)
      gc()
    }
    
    Data[, c("biomass1", "qmd1", "hl1", "ba1", "vol1") := NULL]
    data.table::setnames(Data, old = c("biomass2", "qmd2", "hl2", "ba2", "vol2"), new = c("biomass1", "qmd1", "hl1", "ba1", "vol1"))
    gc()
  }
}
############################################################################################################################

################################################################################################
