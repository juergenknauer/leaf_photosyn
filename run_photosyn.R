rm(list=ls())
setwd("/Users/30044953/Documents/CABLE/_gm_acclim_coord/paper/version_ScienceAdvances/Revision1/Toy_Model")

## Load functions
source("photosyn_functions.R")

## Meteorological conditions
Tair     <- 25
Pressure <- 100
VPD      <- 1.5
Ca       <- seq(100,1500)
PPFD     <- 1500
Tgrowth  <- 15   # ignored if temp_acclim == FALSE
Thome    <- 20   # ignored if temp_acclim == FALSE
beta_b = beta_m = beta_s  <- 1    # water stress scalars (0-1; 1 = no water stress)


## Photosynthetic parameters
Vcmax25 <- 40
bjv25   <- 1.8
Jmax25  <- bjv25*Vcmax25
g0      <- 0.005
g1      <- 3.4
gm25    <- 0.2  # ignored if explicit_gm == FALSE


## Settings
gs_model    <- "USO"
explicit_gm <- FALSE
temp_acclim <- FALSE

out <- photosyn(Tair = Tair, Pressure = Pressure, VPD = VPD, Ca = Ca, PPFD = PPFD,
                Vcmax25 = Vcmax25, Jmax25 = Jmax25, g0 = g0, g1 = g1,  
                explicit_gm = explicit_gm, temp_acclim = temp_acclim,
                gm25 = gm25, beta_s = beta_s, beta_m = beta_m, beta_b = beta_b,
                Tgrowth = Tgrowth, Thome = Thome, 
                gs_model=gs_model,constants=photosyn.constants())




plot(out$An ~ out$Ci)
