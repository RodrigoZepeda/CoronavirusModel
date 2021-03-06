rm(list = ls())

#Required libraries to run the model
library(deSolve)
library(tidyverse)
library(sfsmisc)
library(kableExtra)
library(cowplot)
library(R0)
library(plyr)
library(mgcv)

#Working directory
dir <- "~/Dropbox/Coronavirus/ModelbyCategoriesMarch30"
params.file   <- "parameters.csv"
resource.file <- "resources.csv"
setwd(dir)

source("modelo.edad.R")
source("read.parameters.R")
source("auxiliary.functions.R")

#Parameter list
initial.list <- read.parameters(params.file)
params       <- initial.list$params
state        <- initial.list$state

#You can run model with run model.instance
#microbenchmark::microbenchmark(
params$saturationI2 <- 0.001
params$saturationI3 <- 0.0001
model.1 <- run.model.continuous(params, state,  init.time = 0, end.time = 100)


#And then add quarantine
state   <- quarantine_all(model.1$state)
model.2 <- run.model.continuous(params, state, init.time = 20, end.time = 30)

#Then continue model wth periodic quarantine
state   <- unquarantine_all(model.2$state)
model.3 <- run.model.periodic(params, state, init.time = 30, end.time = 50,
                                  periodicity = 7, days = 2)

#Change R0
params  <- rnought.change.gamma.1(1.2, params)
model.4 <- run.model.continuous(params, model.3$state, init.time = 50, end.time = 60)

#And then add quarantine
state   <- quarantine_all(model.4$state)
params  <- rnought.change.gamma.1(2.1, params)
model.5 <- run.model.continuous(params, state, init.time = 60, end.time = 120)



ggplot() +
  geom_line(aes(x = time, y = I1 + QI, 
                color = "Evolución usual"), data = model.1$dats) + 
  geom_line(aes(x = time, y = I1 + QI,
                color = "Cuarentena poblacional"), 
            data = model.2$dats) + 
  geom_line(aes(x = time, y = I1 + QI,
                color = "Cuarentena periódica"), data = model.3$dats) + 
  geom_line(aes(x = time, y = I1 + QI,
                color = "Evolución usual con R0 = 1.2"), data = model.4$dats) +
  geom_line(aes(x = time, y = I1 + QI,
                color = "Cuarentena poblacional  con R0 = 2.1"), 
            data = model.5$dats) +
  theme_classic()



