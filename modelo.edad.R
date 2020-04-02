#Ecuacion diferencial del modelo a resolverse
#con deSolve

model <- function(t, y, parameters){
  with(as.list(c(y, parameters)),{
    
    #Variables are all stored in a vector here are the
    #key pairs:
    #S  <- y[1:4]
    #E  <- y[5:8]
    #A  <- y[9:12]
    #I1 <- y[13:16]
    #I2 <- y[17:20]
    #I3 <- y[21:24]
    #M  <- y[25:28]
    #Q  <- y[29:32]
    #QA <- y[33:36]
    #QI <- y[37:40]
    
    #Position for each vector
    Spos   <- 0
    Epos   <- 4
    Apos   <- 8
    I1pos  <- 12
    I2pos  <- 16
    I3pos  <- 20
    Mpos   <- 24
    QSpos  <- 28
    QEpos  <- 32
    QApos  <- 36
    QIpos  <- 40
    
    
    #Instantiate variables for age
    dy  <- rep(0, 10*age_cats)
    
    #Get total infected for catefory
    I2.total <- sum(y[(I2pos + 1):(I2pos + age_cats)])
    I3.total <- sum(y[(I3pos + 1):(I3pos + age_cats)])
    
    #Check if UCI are saturated
    if (I3.total > saturationI3){
      p.M <- p.M.saturated
    } else {
      p.M <- p.M.not.saturated
    }
    
    if (I2.total > saturationI2){
      p.I3      <- p.I3.saturated
      beta.2to3 <- beta.2to3.saturated
      beta.2toR <- beta.2toR.saturated
    } else {
      p.I3      <- p.I3.not.saturated
      beta.2to3 <- beta.2to3.not.saturated
      beta.2toR <- beta.2toR.not.saturated
    }
    
    #Loopeamos para la tasa de contagio
    Isum <- 0
    for (k in 1:age_cats){
      Isum <- Isum + 
        gamma.1[k]*y[I1pos + k] + gamma.2[k]*y[I2pos + k] + 
        gamma.3[k]*y[I3pos + k] + gamma.A[k]*y[Apos + k]
    }
    
    
   #Calculate for each age group
    for (k in 1:age_cats){
      
      #Susceptibles
      #dS/dt[k] = -Isum*S[k]
      dy[Spos + k]  <- -Isum*y[Spos + k]
      
      #Latentes
      #dE/dt[k] = Isum*S[k] - gamma.E[k]*E[k]
      dy[Epos + k]  <- Isum*y[Spos + k] - gamma.E[k]*y[Epos + k]
      
      #Asintomáticos
      #dA/dt[k] = p.A[k]*gamma.E[k]*E[k] - beta.AtoR[k]*A[k]
      dy[Apos + k]  <- p.A[k]*gamma.E[k]*y[Epos + k] - beta.AtoR[k]*y[Apos + k]
      
      #Infectados 1 (leves)
      #dI1/dt[k] =  (1 - p.A[k])*gamma.E[k]*E[k] - (p.I2[k]*beta.1to2[k] + (1 - p.I2[k])*beta.1toR[k])*I1[k]
      dy[I1pos + k] <- (1 - p.A[k])*gamma.E[k]*y[Epos + k] - (p.I2[k]*beta.1to2[k] + 
        (1 - p.I2[k])*beta.1toR[k] - q2[k])*y[I1pos + k]
      
      #Infectados 2 (graves)
      #dI2/dt[k] = p.I2[k]*beta.1to2[k]*(Q2[k] + I1[k]) - (p.I3[k]*beta.2to3[k] + (1 - p.I3[k])*beta.2toR[k])*I2[k]
      dy[I2pos + k] <- p.I2[k]*beta.1to2[k]*(y[QIpos + k] + y[I1pos + k]) - 
        (p.I3[k]*beta.2to3[k] + (1 - p.I3[k])*beta.2toR[k])*y[I2pos + k]
      
      #Infectados 3 (críticos)
      #dI3/dt[k] = p.I3[k]*beta.2to3[k]*I2[k] - ((1 - p.M[k])*beta.3toR[k] + p.M[k]*beta.3toM[k])*I3[k]
      dy[I3pos + k] <- p.I3[k]*beta.2to3[k]*y[I2pos + k] - 
        ((1 - p.M[k])*beta.3toR[k] + p.M[k]*beta.3toM[k])*y[I3pos + k]
      
      #Muertos
      #dM/dt[k] = p.M[k]*beta.3toM[k]*I3[k]
      dy[Mpos + k]  <- p.M[k]*beta.3toM[k]*y[I3pos + k]
      
      #Cuarentena susceptibles
      dy[QSpos + k] <- 0
      
      #Cuarentena de expuestos
      dy[QEpos + k] <- - gamma.E[k]*y[QEpos + k]
      
      #Cuarentena de asintomáticos
      dy[QApos + k] <- p.A[k]*gamma.E[k]*y[QEpos + k] - beta.AtoR[k]*y[QApos + k]
      
      #Cuarentena de infectados
      dy[QIpos + k] <- (1 - p.A[k])*gamma.E[k]*y[QEpos + k] - 
        (p.I2[k]*beta.1to2[k] + (1 - p.I2[k])*beta.1toR[k])*y[QIpos + k]
      
    }
    
    list(dy)
  })
}