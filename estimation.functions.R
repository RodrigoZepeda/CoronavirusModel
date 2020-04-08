#https://stackoverflow.com/questions/17922637/prediction-intervals-for-poisson-regression-on-r
bootSimFun <- function(preddata,fit,data) {
  bdat  <- data[sample(seq(nrow(data)),size=nrow(data),replace=TRUE),]
  bfit  <- update(fit,data=bdat)
  bpred <- predict(bfit,type="response",newdata=preddata)
  rpois(length(bpred),lambda=bpred)
}

estimate.R0 <- function(cummulative.cases, date.cases){
  
  #Calcular el R0
  GT.flu <- generation.time("weibull", c(1,1))
  res.R  <- estimate.R(cummulative.cases, GT=GT.flu, methods=c("TD"),
                       S0 = 1, t = date.cases,
                       date.first.obs = min(date.cases))
  
  R0     <- data.frame(time = as.Date(names(res.R$estimates$TD$R)),
                       y    = res.R$estimates$TD$R, 
                       ymin = res.R$estimates$TD$conf.int$lower, 
                       ymax =res.R$estimates$TD$conf.int$upper)
  
  R0 <- R0[1:(nrow(R0)-1),]
  
  last_date <- max(R0$time) -1
  
  if (R0$y[nrow(R0)] < 1){warning("DATA RESULTS IN R0 < 1")}
  
  colnames(R0) <- c("Fecha","R0", "Lower bound R0","Upper bound R0")
  rownames(R0) <- c()
  
  return(R0)
}

plot.R0 <- function(R0,
                    title = "Número básico de reproducción (R0) de COVID-19",
                    xlab  = "Días desde el inicio de la epidemia",
                    ylab = "R0"){
    ggplot(R0) +
      theme_classic() +
      geom_line(aes(Fecha, R0), color = "deepskyblue4") +
      scale_y_continuous(breaks = seq(0, 5, 0.5)) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "firebrick") +
      geom_ribbon(aes(ymin = `Lower bound R0`, ymax= `Upper bound R0`, x = Fecha), alpha = 0.2,
                  fill = "deepskyblue") +
      ggtitle(title) +
      xlab(xlab) +
      ylab(ylab) 

}

estimate.daily.cases <- function(dats, dias.predict = 10, confidence = 0.95,
                                 nsim = 1000){
  
  intervals <- NA
  if ("n" %in% colnames(dats) & "Fecha" %in% colnames(dats)){
    
    dats$x <- as.Date(dats$Fecha) - as.Date(dats$Fecha[1])
    
    #Modelo predictivo
    total.cases  <- glm(data=dats, n~as.numeric(x), family="poisson")
    newdata      <- data.frame(x = 0:(max(dats$x) + dias.predict))
    simvals      <- raply(nsim, bootSimFun(preddata=newdata, fit = total.cases,
                                           data=dats))
    
    #Intervalos del modelo
    intervals <- t(apply(simvals, 2, quantile, 
                         c((1 - confidence)/2, 0.5, 1 - (1 - confidence)/2)))
    intervals           <- as.data.frame(intervals)
    colnames(intervals) <- c("Lower","Median","Upper")
    intervals$Fecha     <- as.Date(newdata$x, 
                                   origin = as.Date(dats$Fecha[1]))
    
    intervals <- intervals %>% left_join(dats %>% dplyr::select(-x), by = "Fecha")
    
  }
  
  return(intervals)
  
}

#Gráfica de casos incidentes
estimate.cummulative.cases <- function(dats, dias.predict = 10, nsim = 100){
  
  dats <- dats[,c("Fecha","Cummulative")]
  
  dats.predicted <- dats
  if ("Cummulative" %in% colnames(dats) & "Fecha" %in% colnames(dats)){
    
    #Ajuste del modelo
    modelo.fit    <- gam(Cummulative ~ s(as.numeric(Fecha)), data = dats)
    predict.dates <- seq(as.numeric(dats$Fecha[1]), 
                         as.numeric(dats$Fecha[nrow(dats)]) + dias.predict, by = 1)
    dats.pred     <- data.frame(Fecha = predict.dates)
    predichos     <- predict.gam(modelo.fit, dats.pred, se.fit = T) 
    
    #Generate prediction intervals
    beta  <- coef(modelo.fit)
    V     <- vcov(modelo.fit)
    Cv    <- chol(V)
    nus   <- rnorm(nsim * length(beta))
    beta_sims   <- beta + t(Cv) %*% matrix(nus, nrow = length(beta), ncol = nsim)
    sim_idx     <- sample.int(nrow(dats.pred), size = nsim, replace = TRUE)
    sim_dat     <- dats.pred[sim_idx, c("fecha")]
    covar_sim   <- predict(modelo.fit, newdata = dats.pred, type = "lpmatrix")
    linpred_sim <- covar_sim %*% beta_sims
    invlink     <- function(x){x}
    exp_val_sim <- invlink(linpred_sim)
    
    y_sim <- matrix(rnorm(n    = prod(dim(exp_val_sim)), 
                          mean = exp_val_sim, 
                          sd   = summary(modelo.fit)$scale), 
                    nrow = nrow(exp_val_sim), 
                    ncol = ncol(exp_val_sim))
    pred_int_sim  <- apply(y_sim, 1, quantile, prob = c(.025, 0.5, 0.975))
    sim_dat_x0ord <- order(dats.pred$Fecha)
    ymin          <- pred_int_sim[1L, sim_dat_x0ord] 
    y             <- pred_int_sim[2L, sim_dat_x0ord]
    ymax          <- pred_int_sim[3L, sim_dat_x0ord]
    
    #Filter data < 0
    y[which(y < 0)]       <- 0
    ymin[which(ymin < 0)] <- 0
    ymax[which(ymax < 0)] <- 0
    
    #Save to data frame
    dats.predicted <- data.frame(Fecha = as.Date(dats.pred$Fecha - as.numeric(dats$Fecha[1]), origin = dats$Fecha[1]),
                                 ymin = ymin,
                                 y    = y,
                                 ymax = ymax)
    
    #Cummulative dataset
    dats.predicted <- left_join(dats.predicted, dats, by = "Fecha")
    
    colnames(dats.predicted) <- c("Fecha","Intervalo bajo", 
                                  "Casos predichos",
                                  "Intervalo superior", "Observados")
    
    for (i in 2:nrow(dats)){
      if (is.na(dats.predicted$Observados[i])){
        dats.predicted$Observados[i] <- dats.predicted$Observados[i-1]
      }
    }
  }
  return(dats.predicted)
}

plot.cummulative.cases <- function(cummulative.cases,
                                   title = "Casos confirmados acumulados de COVID-19 por fecha de reporte",
                                   xlab = "Fecha de reporte",
                                   ylab = "Número de casos nuevos confirmados"){
  ggplot(cummulative.cases) + 
    geom_ribbon(aes(x    = Fecha, 
                    ymin = `Intervalo bajo`,
                    ymax = `Intervalo superior`), alpha = 0.1) +
    geom_line(aes(x = Fecha, y = `Casos predichos`), color = "deepskyblue4") +
    geom_point(aes(x = Fecha, y = Observados, color = "Observados"), size = 2) +
    theme_bw() +
    ggtitle(title) +
    xlab(xlab) +
    ylab(ylab)+
    scale_color_manual("Datos", values = c("red", "deepskyblue4")) 
}

plot.predicted.cases <- function(predicted.cases,
                                 title = "Casos confirmados de COVID-19 por fecha de reporte",
                                 xlab  = "Fecha de reporte",
                                 ylab  = "Número de casos nuevos confirmados"){
  
  ggplot(predicted.cases, aes(x = Fecha)) +
    geom_col(aes(y = Median, fill = "Predichos")) +
    geom_line(aes(y = Lower, color = "Intervalo"), size = 1, linetype = "dotted") +
    geom_line(aes(y = Upper, color = "Intervalo"), size = 1, linetype = "dotted") +
    geom_point(aes(y = n, color = "Observados"), size = 2) +
    scale_color_manual("Datos", values = c("Observados" = "red","Predichos" = "deepskyblue4", "Intervalo" = "orange")) +
    scale_fill_manual("Datos", values = c("Observados" = "red","Predichos" = "deepskyblue4", "Intervalo" = "orange")) +
    theme_bw() +
    ggtitle(title) + xlab(xlab) + ylab(ylab)
    
}

model.from.cummulative.cases <- function(cummulative.cases, date.cases,
                                         params.file, initial.state = NULL){
  
  #Estimate R0
  R0 <- estimate.R0(cummulative.cases, date.cases)#[nrow(R0)-1,"R0"]
  R0 <- R0[nrow(R0),"R0"]
  
  #Get model parameters
  initial.list <- read.parameters(params.file)
  params       <- initial.list$params
  
  if (!is.null(initial.state)){
    state <- initial.state
  } else {
    state <- initial.list$state  
  }
  
  #Change parameters from R0
  params       <- rnought.change.gamma.1(R0, params)
  params       <- rnought.change.gamma.A(R0, params)
  
  return(list(params = params, state = state))
  
}

