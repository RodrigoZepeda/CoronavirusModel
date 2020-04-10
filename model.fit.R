#Sum of squares for minimization
ssq <- function(params, 
                init.state = state, 
                dats,
                t.extra = 10,
                t.step = 1,
                scale.factor = 10000,
                age_cats = length(init.state)/11,
                model.vars = c(paste0("I1", 1:age_cats), 
                               paste0("I2", 1:age_cats), 
                               paste0("I3", 1:age_cats), 
                               paste0("M", 1:age_cats))){
  
  
  #save dats as data frame
  dats <- as.data.frame(dats)
  
  # solve ODE for a given set of parameters
  ode.model <- run.model.continuous(params    = params, 
                                    state     = init.state,  
                                    init.time = as.numeric(min(dats$time)), 
                                    end.time  = as.numeric(max(dats$time)) + t.extra,
                                    t.step    = t.step)
  
  
  # Filter data that contains time points where data is available
  cummulative.model <- cummulative.cases(ode.model, age_cats = age_cats)
  cummulative.model <- cummulative.rescale(cummulative.cases(ode.model), 
                                           scale.factor = scale.factor)
  model.data        <- cummulative.model %>% filter(time %in% dats$time)
  
  # Evaluate predicted vs experimental residual
  ssqres <-  as.vector(model.data[,model.vars] - dats[,model.vars])
  
  
  # return predicted vs experimental residual
  return(ssqres)
  
}

ssq.fun <- function(params.optimize){
  
  #Optimización de parámetros de contacto
  params$gamma.1 <-  rep(params.optimize[1], 4)
  params$gamma.A <-  rep(params.optimize[1], 4)
  
  ssq(params, state, dats, t.extra = 10, t.step = 1,
      scale.factor = params.optimize[2], 
      length(state)/11, model.vars = c("I1"))
  
}

fit.model <- function(dats, state, 
                      init.pars = c(1/2, 1), 
                      maxiter = 100,
                      lower = rep(0, length(init.pars)), 
                      upper = rep(Inf, length(init.pars))){
  dats <- dats; state <- state
  fitval <- nls.lm(par = init.pars, fn = ssq.fun, 
                   lower = lower, upper = upper,
                   control = nls.lm.control(maxiter = maxiter))
  return(fitval)
}

run.fitted.model <- function(fitval, state, params, age_cats = length(state)/11, 
                             maxiter = 100, init.time = 0, end.time = 100, t.step = 0.5){
  
  #Fitted parameter values
  fitted.params <- coef(fitval)
  
  #Optimización de parámetros de contacto
  params$gamma.1 <- rep(fitted.params[1],4) #c(fitted.params[1], fitted.params[2], fitted.params[3], fitted.params[4])
  params$gamma.A <- rep(fitted.params[1],4) #c(fitted.params[1], fitted.params[2], fitted.params[3], fitted.params[4])
  
  # solve ODE for a given set of parameters
  ode.model <- run.model.continuous(params    = params, 
                                    state     = state,  
                                    init.time = init.time, 
                                    end.time  = end.time,
                                    t.step    = t.step)
  
  
  # Filter data that contains time points where data is available
  cummulative.model <- cummulative.cases(ode.model, age_cats = age_cats)
  
  cummulative.model <- cummulative.rescale(cummulative.cases(ode.model), 
                                           scale.factor = fitted.params[2])
  ode.model         <- model.rescale(ode.model, scale.factor = fitted.params[2])
  
  return(list(cummulative = cummulative.model, incidence = ode.model))
}

# fitval       <- fit.model()
#model.info.1 <- run.fitted.model(fitval, state, params, end.time = 30)
#model.info.2 <- run.fitted.model(fitval, state, params, end.time = 150)

# ggplot.epidemiological.lines.infected(model.info.2$incidence)
# ggplot.data.to.fitted(model.info$cummulative, dats)

