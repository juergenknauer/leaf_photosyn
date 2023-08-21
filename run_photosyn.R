rm(list=ls())
setwd("/Users/30044953/Documents/CABLE/_gm_acclim_coord/leaf_photosyn")

## Load functions
source("photosyn_functions.R")

## Meteorological conditions
Tair     <- 25     # Air temperature (degC)
Pressure <- 100    # Atmospheric Pressure (kPa)
VPD      <- 1.5    # Vapour pressure deficit (kPa)
Ca       <- 410    # Atmospheric CO2 concentration (umol mol-1)
PPFD     <- 1500   # Photosynthetically active photon flux density (umol m-2 s-1)
Tgrowth  <- 15     # Growth temperature (ignored if temp_acclim == FALSE)
Thome    <- 20     # Home temperature (ignored if temp_acclim == FALSE)
beta_b = beta_m = beta_s  <- 1  # water stress scalars (0-1; 1 = no water stress)

## Photosynthetic parameters
Vcmax25 <- 40      # Maximum carboxylation rate at 25 degC (umol m-2 s-1)
bjv25   <- 1.8     # Jmax25/Vcmax25 ratio
Jmax25  <- bjv25*Vcmax25 # Maximum electron transport rate (umol m-2 s-1)
g0      <- 0.005   # residual stomatal conductance (mol m-2 s-1)
g1      <- 3.4     # stomatal slope parameter (kPa^0.5)
gm25    <- 0.2     # mesophyll conductance at 25 degC (ignored if explicit_gm == FALSE)

## Settings
gs_model    <- "USO"
explicit_gm <- FALSE
temp_acclim <- TRUE

out <- photosyn(Tair = Tair, Pressure = Pressure, VPD = VPD, Ca = Ca, PPFD = PPFD,
                Vcmax25 = Vcmax25, Jmax25 = Jmax25, g0 = g0, g1 = g1,  
                explicit_gm = explicit_gm, temp_acclim = temp_acclim,
                gm25 = gm25, beta_s = beta_s, beta_m = beta_m, beta_b = beta_b,
                Tgrowth = Tgrowth, Thome = Thome, 
                gs_model=gs_model,constants=photosyn.constants())



## --------------------------------------------------------------------------------
## ---------------------------   Example applications   ---------------------------
## --------------------------------------------------------------------------------

### 1) Compare An with and without thermal acclimation
Tair <- seq(10,35,0.1)
acclim_off <- photosyn(Tair = Tair, Pressure = Pressure, VPD = VPD, Ca = Ca, PPFD = PPFD,
                      Vcmax25 = Vcmax25, Jmax25 = Jmax25, g0 = g0, g1 = g1,  
                      explicit_gm = FALSE, temp_acclim = FALSE,
                      gm25 = gm25, beta_s = beta_s, beta_m = beta_m, beta_b = beta_b,
                      Tgrowth = NULL, Thome = NULL, 
                      gs_model=gs_model,constants=photosyn.constants())

acclim_on <- photosyn(Tair = Tair, Pressure = Pressure, VPD = VPD, Ca = Ca, PPFD = PPFD,
                       Vcmax25 = Vcmax25, Jmax25 = Jmax25, g0 = g0, g1 = g1,  
                       explicit_gm = FALSE, temp_acclim = TRUE,
                       gm25 = gm25, beta_s = beta_s, beta_m = beta_m, beta_b = beta_b,
                       Tgrowth = 25, Thome = 30, 
                       gs_model=gs_model,constants=photosyn.constants())

## plot total An
plot(acclim_off$An ~ Tair,ylab="Net photosynthesis (umol m-2 s-1)",xlab="Leaf temperature (degC)",ylim=c(4,12),type="l")
points(acclim_on$An ~ Tair,col="red",type="l")


## plot sensitivity of An to Tleaf (dAn/An)/(dTair/Tair)
gamma_on  <- diff(acclim_on$An)/acclim_on$An[-length(acclim_on$An)]/(median(diff(Tair))/Tair[-length(Tair)])
gamma_off <- diff(acclim_off$An)/acclim_off$An[-length(acclim_off$An)]/(median(diff(Tair))/Tair[-length(Tair)])

## Note that the following plot shows the temperature sensitivity of An rather than Vcmax or Jmax as in the paper
## An shows a lower temperature optimum compared to Vcmax and Jmax.
plot(gamma_on ~ Tair[-length(Tair)],ylab="beta (dAn/An)/(dTair/Tair)",xlab="Air temperature (degC)")
points(gamma_off ~ Tair[-length(Tair)],col="red")






### 2) Compare An with and without explicit mesophyll conductance
Tair <- 25
Ca   <- seq(100,1500)
gm_off <- photosyn(Tair = Tair, Pressure = Pressure, VPD = VPD, Ca = Ca, PPFD = PPFD,
                    Vcmax25 = Vcmax25, Jmax25 = Jmax25, g0 = g0, g1 = g1,  
                    explicit_gm = FALSE, temp_acclim = FALSE,
                    gm25 = gm25, beta_s = beta_s, beta_m = beta_m, beta_b = beta_b,
                    Tgrowth = 25, Thome = 25, 
                    gs_model=gs_model,constants=photosyn.constants())

gm_on <- photosyn(Tair = Tair, Pressure = Pressure, VPD = VPD, Ca = Ca, PPFD = PPFD,
                  Vcmax25 = Vcmax25, Jmax25 = Jmax25, g0 = g0, g1 = g1,  
                  explicit_gm = TRUE, temp_acclim = FALSE,
                  gm25 = gm25, beta_s = beta_s, beta_m = beta_m, beta_b = beta_b,
                  Tgrowth = 25, Thome = 25, 
                  gs_model=gs_model,constants=photosyn.constants())


# plot net photosynthesis (An) against intercellular CO2 concentration (Ci)
plot(gm_off$An ~ gm_off$Ci,type="l",ylab="Net photosynthesis (umol m-2 s-1)",xlab="Atm. CO2 concentration (umol mol-1)")
points(gm_on$An ~ gm_on$Ci,col="red",type="l")

# plot sensitivty of An to changes in Ci (dAn/An)/(dCa/Ca)
beta_on  <- diff(gm_on$An)/gm_on$An[-length(gm_on$An)]/(median(diff(Ca))/Ca[-length(Ca)])
beta_off <- diff(gm_off$An)/gm_off$An[-length(gm_off$An)]/(median(diff(Ca))/Ca[-length(Ca)])

plot(beta_on ~ Ca[-length(Ca)],ylab="beta (dAn/An)/(dCa/Ca)",xlab="Atm. CO2 concentration (umol mol-1)")
points(beta_off ~ Ca[-length(Ca)],col="red")





### 3) Compare An with different values of bjv (keeping Neff constant)
Tair    <- 25
Ca      <- seq(100,1500)

# photosynthetic parameters
bjv1    <- 1.8
Vcmax25 <- 40
Jmax25  <- bjv1 * Vcmax25 

NcostJV <- 2
bjv2    <- 1.4
Neff    <- Vcmax25 + NcostJV * bjv1 * Vcmax25
Vcmax25_2   <- Neff / (1.0 + NcostJV * bjv2)
Jmax25_2    <- bjv2 * Vcmax25_2

bjv  <- 1.8
bjv_def <- photosyn(Tair = Tair, Pressure = Pressure, VPD = VPD, Ca = Ca, PPFD = PPFD,
                    Vcmax25 = Vcmax25, Jmax25 = Jmax25, g0 = g0, g1 = g1,  
                    explicit_gm = FALSE, temp_acclim = FALSE,
                    gm25 = gm25, beta_s = beta_s, beta_m = beta_m, beta_b = beta_b,
                    Tgrowth = NULL, Thome = NULL, 
                    gs_model=gs_model,constants=photosyn.constants())

bjv <- 1.6
bjv_low <- photosyn(Tair = Tair, Pressure = Pressure, VPD = VPD, Ca = Ca, PPFD = PPFD,
                    Vcmax25 = Vcmax25_2, Jmax25 = Jmax25_2, g0 = g0, g1 = g1,  
                    explicit_gm = FALSE, temp_acclim = FALSE,
                    gm25 = gm25, beta_s = beta_s, beta_m = beta_m, beta_b = beta_b,
                    Tgrowth = NULL, Thome = NULL, 
                    gs_model=gs_model,constants=photosyn.constants())


plot(bjv_def$An ~ Ca,type="l")
points(bjv_low$An ~ Ca,col="red",type="l")


## We can also calculate the fraction of An that is Rubisco-limited (analogous to fGPPC)
fGPPC_def <- sum(bjv_def$Rubisco_lim)/nrow(bjv_def)
fGPPC_low <- sum(bjv_low$Rubisco_lim)/nrow(bjv_low)

print(paste("default_bjv:", round(fGPPC_def,3)))
print(paste("low_bjv:",round(fGPPC_low,3)))
## --> fGPPC is significantly lower with a lower bjv!









