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
model.aggregate.data <- function(dats, age_cats = ncol(dats)/11 - 1, 
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
  end.state      <- dats[nrow(dats),] %>% dplyr::select(-time)
  
  #Loop agregating data
  age_cats <- ncol(dats)/11
  for (var in c("S","E","A","I1","I2","I3","QS","QE","QA","QI","M")){
    dats <- model.aggregate.data(dats, variable = var, age_cats = age_cats)
  }
  
  return(list(dats = dats, state = end.state, params = params))
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
                       Susceptible  = "QS",
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
                       Susceptible  = "QS",
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

rnought.change.gamma.1 <- function(R0, params){
  params$gamma.1 <- R0*params$gamma.E
  return(params)
}

rnought.change.gamma.A <- function(R0, params){
  params$gamma.A <- R0*params$gamma.E
  return(params)
}

rnought.change.gamma.2 <- function(R0, params){
  params$gamma.2 <- R0*params$gamma.E
  return(params)
}

rnought.change.gamma.3 <- function(R0, params){
  params$gamma.3 <- R0*params$gamma.E
  return(params)
}

rnought.change.gamma.all <- function(R0, params){
  params <- rnought.change.gamma.1(R0, params)
  params <- rnought.change.gamma.2(R0, params)
  params <- rnought.change.gamma.3(R0, params)
  params <- rnought.change.gamma.A(R0, params)
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
    dplyr::select(-time, -S, -E, -A, -I1, -I2, -I3, -QS, -QE, -QA, -QI, -M)
  
  
  return(list(dats = dats.temp, state = end.state, params = params))
}



cummulative.cases.I1.k <- function(modelo, k){
      
  params <- modelo$params
  dats   <- modelo$dats
  I1.cummulative.k <- with(params,{
    (1 - p.A[k])*gamma.E[k]*(dats[,paste0("E",k)] + dats[,paste0("QE",k)])
  })
  I1.cummulative.k <- integrate.cummulative(modelo$dats$time, I1.cummulative.k)
  return(I1.cummulative.k)
  
}

cummulative.cases.I1 <- function(modelo, age_cats = ncol(modelo$dats)/11 - 1){
  
  I1.cummulative <- rep(0, nrow(modelo$dats))
  
  for (k in 1:age_cats){
    I1.cummulative <- I1.cummulative + cummulative.cases.I1.k(modelo, k)
  }
  
  return(I1.cummulative)
  
}

integrate.cummulative <- function(time, variable){
  integration <- rep(0, length(time))
  if (length(time) > 1){
    for (i in 2:length(time)){
      integration[i] <- integrate.xy(time[1:i], variable[1:i])
    }
  }
  return(integration)
}

cummulative.cases.I2.k <- function(modelo, k){
  
  params <- modelo$params
  dats   <- modelo$dats
  I2.cummulative.k <- with(params,{
    p.I2[k]*beta.1to2[k]*(dats[,paste0("QI",k)] + dats[,paste0("I1",k)])
  })
  I2.cummulative.k <- integrate.cummulative(modelo$dats$time, I2.cummulative.k)
  return(I2.cummulative.k)
  
}

cummulative.cases.I2 <- function(modelo, age_cats = ncol(modelo$dats)/11 - 1){
  
  I2.cummulative <- rep(0, nrow(modelo$dats))
  
  for (k in 1:age_cats){
    I2.cummulative <- I2.cummulative + cummulative.cases.I2.k(modelo, k)
  }
  
  return(I2.cummulative)
  
}

cummulative.cases.I3.k <- function(modelo, k){
  
  params <- modelo$params
  dats   <- modelo$dats
  
  I3.cummulative.k <- with(params,{
    
    I3.cummulative.k <- 0
    
    for (j in 1:nrow(dats)){
      
      #Check if UCI are saturated
      if (dats[j,"I2"] > saturationI2){
        p.I3      <- p.I3.saturated
        beta.2to3 <- beta.2to3.saturated
      } else {
        p.I3      <- p.I3.not.saturated
        beta.2to3 <- beta.2to3.not.saturated
      }
      
      I3.cummulative.k[j] <- p.I3[k]*beta.2to3[k]*dats[j,paste0("I2",k)]
    }
    
    I3.cummulative.k <- integrate.cummulative(modelo$dats$time, I3.cummulative.k)
  })
  
  return(I3.cummulative.k)
  
}

cummulative.cases.I3 <- function(modelo, age_cats = ncol(modelo$dats)/11 - 1){
  
  I3.cummulative <- rep(0, nrow(modelo$dats))
  
  for (k in 1:age_cats){
    I3.cummulative <- I3.cummulative + cummulative.cases.I3.k(modelo, k)
  }
  
  return(I3.cummulative)
  
}

cummulative.cases.A.k <- function(modelo, k){
  
  params <- modelo$params
  dats   <- modelo$dats
  A.cummulative.k <- with(params,{
    p.A[k]*gamma.E[k]*dats[,paste0("E",k)]
  })
  A.cummulative.k <- integrate.cummulative(modelo$dats$time, A.cummulative.k)
  return(A.cummulative.k)
  
}

cummulative.cases.A <- function(modelo, age_cats = ncol(modelo$dats)/11 - 1){
  
  A.cummulative <- rep(0, nrow(modelo$dats))
  
  for (k in 1:age_cats){
    A.cummulative <- A.cummulative + cummulative.cases.A.k(modelo, k)
  }
  
  return(A.cummulative)
  
}

cummulative.cases.E.k <- function(modelo, k){
  
  params <- modelo$params
  dats   <- modelo$dats
  
  E.cummulative.k <- with(params,{
    #Loopeamos para la tasa de contagio
    Isum <- 0
    for (j in 1:age_cats){
      Isum <- Isum + gamma.1[j]*dats[,paste0("I1",j)] + 
        gamma.2[j]*dats[,paste0("I2",j)] + 
        gamma.3[j]*dats[,paste0("I3",j)] + 
        gamma.A[j]*dats[,paste0("A",j)]
    }
    Isum*dats[,paste0("S",k)]
  })
  E.cummulative.k <- integrate.cummulative(modelo$dats$time, E.cummulative.k)
  return(E.cummulative.k)
  
}

cummulative.cases.E <- function(modelo, age_cats = ncol(modelo$dats)/11 - 1){
  
  E.cummulative <- rep(0, nrow(modelo$dats))
  
  for (k in 1:age_cats){
    E.cummulative <- E.cummulative + cummulative.cases.E.k(modelo, k)
  }
  
  return(E.cummulative)
  
}

cummulative.cases.M.k <- function(modelo, k){
  modelo$dats[,paste0("M",k)]
}

cummulative.cases.M <- function(modelo, age_cats = ncol(modelo$dats)/11 - 1){
  
  M.cummulative <- rep(0, nrow(modelo$dats))
  
  for (k in 1:age_cats){
    M.cummulative <- M.cummulative + cummulative.cases.M.k(modelo, k)
  }
  
  return(M.cummulative)
  
}

model.rescale <- function(modelo, scale.factor = 10000){
  
  #Get colnmaes
  other.cols   <- colnames(modelo$dats)[which(colnames(modelo$dats) != "time")]
  scale.factor <- rep(scale.factor, length(other.cols)/length(scale.factor))
  for (i in 1:length(other.cols)){
    modelo$dats[,other.cols[i]] <- scale.factor[i]*modelo$dats[,other.cols[i]]  
  }
  
  
  return(modelo)
  
}

cummulative.rescale <- function(dats.cummulative, scale.factor = 10000){
  
  #Get colnmaes
  other.cols   <- colnames(dats.cummulative)[which(colnames(dats.cummulative) != "time")]
  scale.factor <- rep(scale.factor, length(other.cols)/length(scale.factor))
  for (i in 1:length(other.cols)){
    dats.cummulative[,other.cols[i]] <- scale.factor[i]*dats.cummulative[,other.cols[i]]
  }
  return(dats.cummulative)
}

cummulative.cases <- function(modelo, age_cats = ncol(modelo$dats)/11 - 1){
  
  dats.cummulative <- data.frame(
    time = modelo$dats[,"time"]
  )
  
  for (k in 1:age_cats){
    dats.cummulative[,paste0("E",k)]  <- cummulative.cases.E.k(modelo, k)
    dats.cummulative[,paste0("A",k)]  <- cummulative.cases.A.k(modelo, k)
    dats.cummulative[,paste0("I1",k)] <- cummulative.cases.I1.k(modelo, k)
    dats.cummulative[,paste0("I2",k)] <- cummulative.cases.I2.k(modelo, k)
    dats.cummulative[,paste0("I3",k)] <- cummulative.cases.I3.k(modelo, k)
    dats.cummulative[,paste0("M",k)]  <- cummulative.cases.M.k(modelo, k)
  }
  
  dats.cummulative <- model.aggregate.cummulative(dats.cummulative)
  
  return(dats.cummulative)
  
}



model.aggregate.cummulative <- function(dats.cummulative,  age_cats = (ncol(dats.cummulative)-1)/6){
  
  for (variable in c("E","A","I1","I2","I3","M")){
    dats.cummulative <- model.aggregate.data(dats.cummulative, 
                                             age_cats = age_cats, variable = variable)
  }
  return(dats.cummulative)
}

ggplot.epidemiological.lines.infected <- function(modelo, title = "Evolución de COVID-19", 
                                         xlab = "Día", 
                                         ylab = "Cantidad de casos", 
                                         date.init = NULL,
                                         infect.cats = c("I1","I2","I3"),
                                         scale = 1){
  
  if (!is.null(date.init)){
    modelo$dats$time <- as.Date(modelo$dats$time, origin = date.init)
  }
  
  plot <- ggplot(modelo$dats) 
  if ("I1" %in% infect.cats){
    plot <- plot + geom_line(aes(x = time, y = scale*(I1 + QI), color = "Infección leve")) 
  }
  
  if ("I2" %in% infect.cats){
    plot <- plot + geom_line(aes(x = time, y = scale*I2, color = "Hospitalización grave")) 
  }
  
  if ("I3" %in% infect.cats){
    plot <- plot + geom_line(aes(x = time, y = scale*I3, color = "Hospitalización en UCI")) 
  }
    
  plot + 
    theme_classic() + xlab(xlab) + ylab(ylab) + ggtitle(title) +
    scale_color_brewer("Etapas de la enfermedad", palette = "Dark2") +
    scale_y_continuous(labels = function(x) format(x, scientific = F))
  
}

ggplot.epidemiological.lines.all <- function(modelo, title = "Evolución de COVID-19", 
                                                     xlab = "Día", 
                                                     ylab = "Proporción de casos", 
                                                     date.init = NULL){
  
  if (!is.null(date.init)){
    modelo$dats$time <- as.Date(modelo$dats$time, origin = date.init)
  }
  
  ggplot(modelo$dats) + 
    geom_line(aes(x = time, y = S  + QS,  color = "Susceptibles")) +
    geom_line(aes(x = time, y = E  + QE, color = "Exposed (latent)")) +
    geom_line(aes(x = time, y = A  + QA, color = "Asymptomatic")) +
    geom_line(aes(x = time, y = I1 + QI, color = "Mild infection")) +
    geom_line(aes(x = time, y = I2, color = "Severe infection")) +
    geom_line(aes(x = time, y = I3, color = "Critical infection")) +
    geom_line(aes(x = time, y = M,  color = "Death")) +
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
      geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("I1",i)] + modelo$dats[,paste0("QI",i)], color = as.character(i), linetype = "Mild infection")) +
      geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("I2",i)], color = as.character(i), linetype = "Severe infection")) +
      geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("I3",i)], color = as.character(i), linetype = "Critical infection")) 
  }
  
  #plot <- plot + scale_color_brewer("Categoría", palette = 1)
  return(plot)
  
}

geom_line.variable.from.list <- function(model.list, 
                                      plot =  ggplot(), varname = "S",
                                      linetype = "solid",
                                      color = names(model.list)){
  
  
  k <- 1
  for (this.model in model.list){
    plot <- plot + geom_line(aes_string(x = this.model$dats[,"time"], 
                                        y = this.model$dats[, varname], 
                                        color = as.factor(color[k])),
                             linetype = linetype, data = this.model$dats)
    k <- k + 1
  }
  return(plot)
}

variable.bind <- function(model.list, varname = "dats"){
  dats <- model.list[[1]][varname]
  if (length(model.list) > 1){
    for (i in 2:length(model.list)){
      dats <- dats %>% bind_rows(model.list[[i]][varname])
    }
  }
  return(dats)
}

params.bind <- function(model.list){
  params <- list(model.list[[1]]["params"])
  if (length(model.list) > 1){
    for (i in 2:length(model.list)){
      params <- append(params, model.list[[i]]["params"])
    }
  }
  names(params) <- names(model.list)
  return(params)
}

model.bind <- function(model.list){
  list(dats  = variable.bind(model.list, "dats"),
       state = variable.bind(model.list, "state"),
       params = params.bind(model.list))
}

cummulative.bind <- function(model.list){
  dats <- cummulative.cases(model.list[[1]])
  if (length(model.list) > 1){
    for (i in 2:length(model.list)){
      acumulados         <- cummulative.cases(model.list[[i]])
      for (k in 2:ncol(acumulados))
      acumulados[,k] <- acumulados[,k] + dats[nrow(dats),k]
      dats <- dats %>% bind_rows(acumulados)
    }
  }
  return(dats)
}


ggplot.epidemiological.lines.all.cat <- function(modelo, title = "Evolución de COVID-19", 
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
              geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("S",i)] + modelo$dats[,paste0("QS",i)],  
                            color =  as.character(i), linetype = "Susceptible")) +
              geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("E",i)] + modelo$dats[,paste0("QE",i)], 
                            color = as.character(i), linetype = "Exposed (latent)")) +
              geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("A",i)] + modelo$dats[,paste0("QA",i)], 
                            color = as.character(i), linetype = "Asymptomatic")) +
              geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("I1",i)] + modelo$dats[,paste0("QI",i)], color = as.character(i), linetype = "Mild infection")) +
              geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("I2",i)], 
                             color = as.character(i), linetype = "Severe infection")) +
              geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("I3",i)], 
                             color = as.character(i), linetype = "Critical infection")) 
              geom_line(aes_(x = modelo$dats$time, y = modelo$dats[,paste0("M",i)],  
                             color = as.character(i), linetype = "Death"))
  }
  
  return(plot)
  
}

ggplot.cummulative.lines.all <- function(dats, title = "Evolución de COVID-19", 
                                         xlab = "Día", 
                                         ylab = "Incidencia acumulada de casos", 
                                         states = ncol(dats)/11 - 1,
                                         vars = c("E","A","I1","I2","I3","M"),
                                         date.init = NULL){
  
  if (!is.null(date.init)){
    dats$time <- as.Date(dats$time, origin = date.init)
  }
  
  plot <- ggplot(dats) + theme_classic()
  
  for (i in 1:states){
    
    if ("E" %in% vars){
    plot <- plot + 
      geom_line(aes_(x = dats$time, y = dats[,paste0("E",i)],
                     linetype = as.character(i),
                     color = "Exposed (latent)"))
    }
    
    if ("A" %in% vars){
      plot <- plot + 
        geom_line(aes_(x = dats$time, y = dats[,paste0("A",i)],
                       linetype = as.character(i),
                       color = "Asymptomatic")) 
    }
    
    if ("I1" %in% vars){
      plot <- plot + 
        geom_line(aes_(x = dats$time, y = dats[,paste0("I1",i)],
                     linetype = as.character(i),
                     color = "Mild infection")) 
    }
    
    if ("I2" %in% vars){
      plot <- plot + 
        geom_line(aes_(x = dats$time, y = dats[,paste0("I2",i)],
                     linetype = as.character(i),
                     color = "Severe infection"))
    }
    
    if ("I3" %in% vars){
      plot <- plot + 
        geom_line(aes_(x = dats$time, y = dats[,paste0("I3",i)], 
                     linetype = as.character(i),
                     color = "Critical infection")) 
    }
    
    if ("M" %in% vars){
      geom_line(aes_(x = dats$time, y = dats[,paste0("M",i)], 
                     linetype = as.character(i),
                     color = "Death"))
    }
  }
  
  return(plot)
  
}

time.to.saturation.variable <- function(model, varname = "I2", 
                                        saturation.value = model$params$saturationI2){
  
  #Check if saturation occurs
  saturation <- which(model$dats[varname] > saturation.value)
  
  if (length(saturation) == 0){
    saturation.time <- Inf
  } else {
    saturation.time <- model$dats$time[min(saturation)]
  }
  
  return(saturation.time)
  
}

time.of.peak.variable <- function(model, varname = "I1"){
  
  #Get peak value
  peak.value   <- max(model$dats[varname])
  
  #Get where peak value occurs
  peak.time <- model$dats$time[which(model$dats[varname] == peak.value)]
  
  #Vector
  peak.vars    <- c(peak.value, peak.time)
  names(peak.vars) <- c(paste("Peak of", varname), paste("Time to peak of", varname))
  
  return(peak.vars)
  
}

time.of.peak <- function(model){
  
  #Calcular infectados totales
  list(
    "Mild"     = time.of.peak.variable(model, "I1"),
    "Severe"   = time.of.peak.variable(model, "I2"),
    "Critical" = time.of.peak.variable(model, "I3")
  )
  
} 

time.to.saturation <- function(model, 
                               saturation.values = c("I2" = model$params$saturationI2,
                                                     "I3" = model$params$saturationI3)){
  
  c("I2" = time.to.saturation.variable(model, "I2", saturation.values["I2"]),
    "I3" = time.to.saturation.variable(model, "I3", saturation.values["I3"]))
  
}


ggplot.data.to.fitted <- function(cummulative.model, dats){
  
  cols.model <- colnames(dats)[which(colnames(dats) != "time")]
  
  plot <- ggplot(cummulative.model, aes(x = time)) + theme_bw() +
    xlab("Incidencia acumulada por cada 10,000 habitantes") +
    ylab("Días desde el inicio de la epidemia")
  
  for (var in cols.model){
    plot <- plot + geom_line(aes_string(y = var)) +
      geom_point(aes_string(y = var), size = 2, color = "red", data = dats)
  }
    
  return(plot)
}


