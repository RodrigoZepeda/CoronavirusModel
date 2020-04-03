#Function to estimate R nought from AR
R0.from.AR  <- function(x, S = 1){
  return(-log((1-x)/S)/(x - (1 - S)))
}

#Function to estimate AR from Rnought
AR.from.R0 <- function(R0, S = 1, lims = c(0.001,0.9999)){
  R0.temp <- function(x){R0.from.AR(x, S = S) - R0}
  return(uniroot(R0.temp, lims)$root)
}

#Function to aggregare data
model.aggregate.data <- function(dats, age_cats = ncol(dats)/11, 
                                 k = 1:age_cats, 
                                 variable = "S"){
  
  #Get column names
  col.names       <- paste0(variable, k)
  dats[,variable] <- rowSums(dats[,col.names])
  return(dats)
  
}

#Function to run the model
run.model.continuous <- function(params, state, init.time = 0, end.time = 20, 
                                 t.step = 0.1, method = "lsoda"){
  
  #Get times vector
  times          <- seq(init.time, end.time, by = t.step)
  
  #Get model
  out            <- ode(y = as.numeric(state), times = times, func = model, 
                          parms = params, method = method)
  dats           <- as.data.frame(out)
  
  #Get colnames to dataset
  colnames(dats) <- c("time", names(state))
  
  #Get ending state
  end.state      <- dats[nrow(dats),] %>% select(-time)
  
  #Loop agregating data
  age_cats <- ncol(dats)/11
  for (var in c("S","E","A","I1","I2","I3","Q","QA","QI","M")){
    dats <- model.aggregate.data(dats, variable = var, age_cats = age_cats)
  }
  
  return(list(dats = dats, state = end.state))
}

#Function to quarantine kth category
quarantine_k <- function(state, age_cats = length(state)/11, k = 1:age_cats,
                         quarantine.type = "Susceptible", quarantine.proportion = 1){
  
  #Check if quarantine proportion is vector or numner
  if (length(quarantine.proportion) == 1 & length(k) > 1){
    quarantine.proportion <- rep(quarantine.proportion, length(k))
  }
  
  p.col.name <- switch(quarantine.type,
                       Susceptible  = "S",
                       Exposed      = "E",
                       Asymptomatic = "A",
                       Mild = "I1")
  
  q.col.name <- switch(quarantine.type,
                       Susceptible  = "Q",
                       Exposed      = "QE",
                       Asymptomatic = "QA",
                       Mild = "QI")
  
  #Get column names
  q.col.name <- paste0(q.col.name,k)
  p.col.name <- paste0(p.col.name,k)
  
  #Loop assigning p.names to q.nm
  for (i in 1:length(k)){
    state[,q.col.name[i]] <- quarantine.proportion[i]*state[,p.col.name[i]] + state[,q.col.name[i]]
    state[,p.col.name[i]] <- (1 - quarantine.proportion[i])*state[,p.col.name[i]]
  }
  
  return(state)
}


#Function to quarantine all
quarantine_all <- function(state, age_cats = length(state)/11, quarantine.proportion = 1){
  
  qtype <- c("Susceptible","Exposed","Asymptomatic","Mild")
  
  if (length(quarantine.proportion) == 1){
    quarantine.proportion <- rep(quarantine.proportion, 3)
  }
  
  for (i in 1:3){
    state <- quarantine_k(state, quarantine.type = qtype[i], quarantine.proportion = quarantine.proportion[i])
  }
  
  return(state)
  
}

#Function to change from quarantine to normal
unquarantine_k <- function(state, age_cats = length(state)/11, k = 1:age_cats,
                         quarantine.type = "Susceptible", unquarantine.proportion = 1){
  
  #Check if quarantine proportion is vector or numner
  if (length(unquarantine.proportion) == 1 & length(k) > 1){
    unquarantine.proportion <- rep(unquarantine.proportion, length(k))
  }
  
  p.col.name <- switch(quarantine.type,
                       Susceptible  = "S",
                       Exposed      = "E",
                       Asymptomatic = "A",
                       Mild = "I1")
  
  q.col.name <- switch(quarantine.type,
                       Susceptible  = "Q",
                       Exposed      = "QE",
                       Asymptomatic = "QA",
                       Mild = "QI")
  
  #Get column names
  q.col.name <- paste0(q.col.name,k)
  p.col.name <- paste0(p.col.name,k)
  
  #Loop assigning p.names to q.nm
  for (i in 1:length(k)){
    state[,p.col.name[i]] <- state[,p.col.name[i]] + unquarantine.proportion[i]*state[,q.col.name[i]]
    state[,q.col.name[i]] <- (1 - unquarantine.proportion[i])*state[,q.col.name[i]]
  }
  
  return(state)
}

#Function to unquarantine all
unquarantine_all <- function(state, age_cats = length(state)/11, unquarantine.proportion = 1){
  
  qtype <- c("Susceptible","Exposed","Asymptomatic","Mild")
  
  if (length(unquarantine.proportion) == 1){
    unquarantine.proportion <- rep(unquarantine.proportion, 3)
  }
  
  for (i in 1:3){
    state <- unquarantine_k(state, quarantine.type = qtype[i], unquarantine.proportion = unquarantine.proportion[i])
  }
  
  return(state)
  
}


#Function to quarantine all k-states 
quarantine.all.k <-function(state, age_cats = length(state)/11, k = 1:age_cats, quarantine.proportion = 1){
  
  qtype <- c("Susceptible","Exposed","Asymptomatic","Mild")
  
  if (length(quarantine.proportion) == 1){
    quarantine.proportion <- rep(quarantine.proportion, 3)
  }
  
  for (i in 1:3){
    state <- quarantine_k(state, quarantine.type = qtype[i], quarantine.proportion = quarantine.proportion[i],k=k)
  }
  
  return(state)
  
}

#Function to unquarantine all k-states 
unquarantine.all.k <- function(state, age_cats = length(state)/11, k = 1:age_cats, unquarantine.proportion = 1){
  
  qtype <- c("Susceptible","Exposed","Asymptomatic","Mild")
  
  if (length(unquarantine.proportion) == 1){
    unquarantine.proportion <- rep(unquarantine.proportion, 3)
  }
  
  for (i in 1:3){
    state <- unquarantine_k(state, quarantine.type = qtype[i], unquarantine.proportion = unquarantine.proportion[i],k=k)
  }
  
  return(state)
  
}


rnought.change.gamma.1 <- function(R0, params){
  params$gamma.1 <- R0/params$gamma.E
  return(params)
}

rnought.change.gamma.2 <- function(R0, params){
  params$gamma.2 <- R0/params$gamma.E
  return(params)
}

rnought.change.gamma.3 <- function(R0, params){
  params$gamma.3 <- R0/params$gamma.E
  return(params)
}

run.model.periodic <- function(params, state, init.time = 0, end.time = 30, 
                                   periodicity = 7, days = 2, 
                                   quarantine.proportion = 1,
                                   t.step = 0.1, method = "lsoda"){
  
  #Warning if too small periodicity
  if (periodicity > end.time - init.time){
    stop("Invalid periodicity for said time. This is no periodic quarantine.")
  }
  
  #Calculate moments of quarantine
  q.moments <- seq(init.time, end.time, by = periodicity)
  dats      <- as.data.frame(state)
  dats$time <- init.time
  
  for (i in 1:(length(q.moments) - 1)){
    
    #Run model without quarantine
    no.quarantine.model <- run.model.continuous(params, state, init.time = q.moments[i], 
                                    end.time = q.moments[i+1] - days, t.step = t.step, method = "lsoda")
    
    #Get new state after quarantine
    state     <- quarantine_all(no.quarantine.model$state, quarantine.proportion = quarantine.proportion)
    
    #Run model with quarantine
    quarantine.model <- run.model.continuous(params, state, init.time = q.moments[i+1] - days, 
                                              end.time = q.moments[i+1], t.step = t.step, method = "lsoda")
    
    #Get new state after quarantine
    state     <- unquarantine_all(quarantine.model$state, unquarantine.proportion = 1)
    
    #Bind to data
    if (!exists("dats.temp")){
      dats.temp  <- rbind(no.quarantine.model$dats, quarantine.model$dats)
    } else {
      dats.temp  <- rbind(dats.temp, no.quarantine.model$dats, quarantine.model$dats)
    }
    
  }
  
  #Check if we haven't arrived
  if (end.time != q.moments[length(q.moments)]){
    no.quarantine.end <- run.model.continuous(params, state, init.time = q.moments[length(q.moments)],  
                                            end.time =end.time, t.step = t.step, method = "lsoda")
    dats.temp         <- rbind(dats.temp, no.quarantine.end$dats)
  }
  
  #Get ending state
  end.state      <- dats.temp[nrow(dats.temp),] %>% 
    select(-time, -S, -E, -A, -I1, -I2, -I3, -Q, -QA, -QI, -M)
  
  
  return(list(dats = dats.temp, state = end.state))
}



#Model with periodic quarantine of the k-class.    
run.model.periodic.k <- function(params, state, init.time = 0, end.time = 30,age_cats = length(state)/11, k = 1:age_cats,   
                                 periodicity = 7, days = 2, quarantine.proportion = 1,
                                 t.step = 0.1, method = "lsoda"){
  #Warning if too small periodicity
  if (periodicity > end.time - init.time){
    stop("Invalid periodicity for said time. This is no periodic quarantine.")
  }
  
  
  #Calculate moments of quarantine
  q.moments <- seq(init.time, end.time, by = periodicity)
  dats      <- as.data.frame(state)
  dats$time <- init.time
  
  for (i in 1:(length(q.moments) - 1)){
    
    #Run model without quarantine
    no.quarantine.model <- run.model.continuous(params, state, init.time = q.moments[i], 
                                                end.time = q.moments[i+1] - days, t.step = t.step, method = "lsoda")
    
    #Get new state after quarantine
    state     <- quarantine.all.k(no.quarantine.model$state, quarantine.proportion = quarantine.proportion, k=k )
    
    #Run model with quarantine
    quarantine.model <- run.model.continuous(params, state, init.time = q.moments[i+1] - days, 
                                             end.time = q.moments[i+1], t.step = t.step, method = "lsoda")
    
    #Get new state after quarantine
    state     <- unquarantine.all.k(quarantine.model$state, unquarantine.proportion = 1,  age_cats = age_cats, k=k )
    
    #Bind to data
    if (!exists("dats.temp")){
      dats.temp  <- rbind(no.quarantine.model$dats, quarantine.model$dats)
    } else {
      dats.temp  <- rbind(dats.temp, no.quarantine.model$dats, quarantine.model$dats)
    }
    
  }
  
  #Check if we haven't arrived
  if (end.time != q.moments[length(q.moments)]){
    no.quarantine.end <- run.model.continuous(params, state, init.time = q.moments[length(q.moments)],  
                                              end.time =end.time, t.step = t.step, method = "lsoda")
    dats.temp         <- rbind(dats.temp, no.quarantine.end$dats)
  }
  
  #Get ending state
  end.state      <- dats.temp[nrow(dats.temp),] %>% 
    select(-time, -S, -E, -A, -I1, -I2, -I3, -Q, -QA, -QI, -M)
  
  
  return(list(dats = dats.temp, state = end.state))
}
#Model with quarentine of tke s-class and periodic quarantine of the k-class .         
run.model.hybrid.ks<- function(params, state, init.time = 0, end.time = 30,age_cats = length(state)/11, 
                               s=min(ceiling(age_cats/2),age_cats),
                               k = seq(1,age_cats)[-s],   
                               periodicity = 7, days = 2, quarantine.proportion_k = 1,
                               quarantine.proportion_s = 1,t.step = 0.1, method = "lsoda"){
  
  #Warning if too small periodicity
  if (periodicity > end.time - init.time){
    stop("Invalid periodicity for said time. This is no periodic quarantine.")
  }
  
  
  #Calculate moments of quarantine
  q.moments <- seq(init.time, end.time, by = periodicity)
  dats      <- as.data.frame(state)
  dats$time <- init.time
  
  #Get state after quarantine of the s-class
  state     <- quarantine.all.k(state, quarantine.proportion = quarantine.proportion_s, k=s )        
  
  for (i in 1:(length(q.moments) - 1)){
    
    #Run model without quarantine of the k-class
    no.quarantine.model <- run.model.continuous(params, state, init.time = q.moments[i], 
                                                end.time = q.moments[i+1] - days, t.step = t.step, method = "lsoda")
    
    #Get new state after quarantine of the k-class 
    state     <- quarantine.all.k(no.quarantine.model$state, quarantine.proportion = quarantine.proportion_k, k=k )
    #Run model with quarantine of the k-class
    quarantine.model <- run.model.continuous(params, state, init.time = q.moments[i+1] - days, 
                                             end.time = q.moments[i+1], t.step = t.step, method = "lsoda")
    
    #Get new state after quarantine of the k-class
    state     <- unquarantine.all.k(quarantine.model$state, unquarantine.proportion = 1, k=k )
    
    #Bind to data
    if (!exists("dats.temp")){
      dats.temp  <- rbind(no.quarantine.model$dats, quarantine.model$dats)
    } else {
      dats.temp  <- rbind(dats.temp, no.quarantine.model$dats, quarantine.model$dats)
    }
    
  }
  #Get the new state after quarantine the s-class. 
  state     <- unquarantine.all.k(quarantine.model$state, unquarantine.proportion = 1,  k=s )
  #Check if we haven't arrived
  if (end.time != q.moments[length(q.moments)]){
    no.quarantine.end <- run.model.continuous(params, state, init.time = q.moments[length(q.moments)],  
                                              end.time =end.time, t.step = t.step, method = "lsoda")
    dats.temp         <- rbind(dats.temp, no.quarantine.end$dats)
  }
  
  #Get ending state
  end.state      <- dats.temp[nrow(dats.temp),] %>% 
    select(-time, -S, -E, -A, -I1, -I2, -I3, -Q, -QA, -QI, -M)
  
  
  return(list(dats = dats.temp, state = end.state))
}



ggplot.epidemiological.lines.infected <- function(modelo, title = "Evolución de COVID-19", 
                                         xlab = "Día", 
                                         ylab = "Proporción de casos", 
                                         date.init = NULL){
  
  if (!is.null(date.init)){
    modelo$dats$time <- as.Date(modelo$dats$time, origin = date.init)
  }
  
  ggplot(modelo$dats) + 
    geom_line(aes(x = time, y = I1 + QI, color = "Mild infection")) +
    geom_line(aes(x = time, y = I2, color = "Severe infection")) +
    geom_line(aes(x = time, y = I3, color = "Critical infection")) +
    theme_classic()
  
}


ggplot.epidemiological.lines.infected.cat <- function(modelo, title = "Evolución de COVID-19", 
                                                  xlab = "Día", 
                                                  ylab = "Proporción de casos",
                                                  states = ncol(modelo$dats)/11 - 1,
                                                  date.init = NULL){
  
  if (!is.null(date.init)){
    modelo$dats$time <- as.Date(modelo$dats$time, origin = date.init)
  }
  
  plot <- ggplot(modelo$dats) + theme_classic()
  
  for (i in 1:states){
    
    plot <- plot + 
      geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("I1",i)], color = as.character(i), linetype = "Mild infection")) +
      geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("I2",i)], color = as.character(i), linetype = "Severe infection")) +
      geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("I3",i)], color = as.character(i), linetype = "Critical infection")) 
  }
  
  #plot <- plot + scale_color_brewer("Categoría", palette = 1)
  return(plot)
  
}



