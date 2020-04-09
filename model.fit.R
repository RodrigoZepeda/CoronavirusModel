#Latin hypercube sampling to fit model
#https://www.r-bloggers.com/learning-r-parameter-fitting-for-models-involving-differential-equations/
require(lhs)
require(minpack.lm)

ssq=function(parms){
  
  # inital concentration
  cinit=c(A=1,B=0,C=0)
  # time points for which conc is reported
  # include the points where data is available
  t=c(seq(0,5,0.1),df$time)
  t=sort(unique(t))
  # parameters from the parameter estimation routine
  k1=parms[1]
  k2=parms[2]
  # solve ODE for a given set of parameters
  out=ode(y=cinit,times=t,func=rxnrate,parms=list(k1=k1,k2=k2))
  
  # Filter data that contains time points where data is available
  outdf=data.frame(out)
  outdf=outdf[outdf$time %in% df$time,]
  # Evaluate predicted vs experimental residual
  preddf=melt(outdf,id.var="time",variable.name="species",value.name="conc")
  expdf=melt(df,id.var="time",variable.name="species",value.name="conc")
  ssqres=preddf$conc-expdf$conc
  
  # return predicted vs experimental residual
  return(ssqres)
  
}