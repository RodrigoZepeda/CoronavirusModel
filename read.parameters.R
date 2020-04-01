#Function for reading dataset resources from .csv format
read.resources <- function(file = "resources.csv"){
  try({
    data.resources           <- read.table(file, header = FALSE, sep = ",", row.names = 1)
    colnames(data.resources) <- c("Amount")
    return(data.resources)
  })
}

#Function for reading dataset parameters from .csv format
read.parameters <- function(params.file = "parameters.csv", estimate.saturation = FALSE, resource.file = "resources.csv"){
  try({
    
    #Read dataset
    data.parameters           <- read.table(params.file, header = TRUE, sep = ",", 
                                            row.names = 1, stringsAsFactors = F)
    colnames(data.parameters) <- as.character(data.parameters["Category interpretation",])
    
    #Get interpretation category
    interpretation.position   <- which(row.names(data.parameters) == "Category interpretation")
    data.parameters           <- data.parameters[-c(interpretation.position),] 
    
    #Get position of R0
    R0.position               <- which(row.names(data.parameters) == "R0")
    latent.position           <- which(row.names(data.parameters) == "Days of latency (with virus but not contagious)")
    
    #Loop through each rate and estimate from R0 if necessary
    for (rate in c("mild","severe","critical","asymptomatic")){
      rate.position   <- which(row.names(data.parameters) == paste("Contact rate for", rate, "cases"))
      for (category.position in 1:ncol(data.parameters)){
        if (data.parameters[rate.position, category.position] == "Estimate from R0"){
          data.parameters[rate.position, category.position]  <-  
            as.numeric(data.parameters[R0.position, category.position]) / as.numeric(data.parameters[latent.position, category.position])
        }
      }
    }
    
    #Change all variables to numeric
    data.parameters[] <- lapply(data.parameters, function(x) as.numeric(as.character(x)))
    
 
    
    #Create parameter list and state list
    parameter.list <- list(
      age_cats                = ncol(data.parameters),
      beta.1to2               = as.numeric(1/data.parameters["Days on mild case before becoming severe",]),  
      beta.AtoR               = as.numeric(1/data.parameters["Days asymptomatic contagious before recovery",]),
      beta.1toR               = as.numeric(1/data.parameters["Days of mild case before recovery",]),
      beta.3toR               = as.numeric(1/data.parameters["Days of critical case before recovery",]),
      beta.3toM               = as.numeric(1/data.parameters["Days of critical case before death",]),
      gamma.1                 = as.numeric(data.parameters["Contact rate for mild cases",]), 
      gamma.2                 = as.numeric(data.parameters["Contact rate for severe cases",]),
      gamma.3                 = as.numeric(data.parameters["Contact rate for critical cases",]), 
      gamma.A                 = as.numeric(data.parameters["Contact rate for asymptomatic cases",]),
      gamma.E                 = as.numeric(data.parameters["Days of latency (with virus but not contagious)",]),
      p.I2                    = as.numeric(data.parameters["Probability of severe given mild",]),
      p.A                     = as.numeric(data.parameters["Proportion of asymptomatic cases",]),
      p.M.saturated           = as.numeric(data.parameters["Probability of death given critical without hospitalization",]),
      p.M.not.saturated       = as.numeric(data.parameters["Probability of death given critical with hospitalization",]),
      p.I3.saturated          = as.numeric(data.parameters["Probability of critical given severe without hospitalization",]),
      p.I3.not.saturated      = as.numeric(data.parameters["Probability of critical given severe with hospitalization",]),
      beta.2to3.saturated     = as.numeric(data.parameters["Days of severe case before critical without hospitalization",]),
      beta.2to3.not.saturated = as.numeric(data.parameters["Days of severe case before critical with hospitalization",]),
      beta.2toR.saturated     = as.numeric(data.parameters["Days of severe case before recovery without hospitalization",]),
      beta.2toR.not.saturated = as.numeric(data.parameters["Days of severe case before recovery with hospitalization",]),
      saturationI2            = Inf, 
      saturationI3            = Inf
    )
    
    #Get resource file
    if (estimate.saturation){
      data.resources <- read.resources(resource.file)
      params$saturationI2 <- (data.resources["Default limit of severe cases","Amount"] + data.resources["Extended limit of severe cases","Amount"])/data.resources["Total population","Amount"]
      params$saturationI3 <- (data.resources["Default limit of critical cases","Amount"] + data.resources["Extended limit of critical cases","Amount"])/data.resources["Total population","Amount"]
    }
    
    state.vec <- c(
      S  = as.numeric(data.parameters["Initial susceptible cases",]),
      E  = as.numeric(data.parameters["Initial latent cases",]),
      I1 = as.numeric(data.parameters["Initial mild cases",]),
      I2 = as.numeric(data.parameters["Initial severe cases",]),
      I3 = as.numeric(data.parameters["Initial critical cases",]),
      A  = as.numeric(data.parameters["Initial asymptomatic cases",]),
      M  = rep(0, ncol(data.parameters)),
      Q  = rep(0, ncol(data.parameters)),  
      QA = rep(0, ncol(data.parameters)), 
      QI = rep(0, ncol(data.parameters))
    )
    
    return(list(params = parameter.list, state = state.vec))
    
  })
  
}
