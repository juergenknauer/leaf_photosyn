
#################
### Constants ###
#################
photosyn.constants <- function(){
  
  list(
    
    R_gas    = 8.314,         # universal gas constant
    Kelvin   = 273.15,        # temperature at 0 degC in Kelvin
    Tref     = 25 + 273.15,   # 25 degC in Kelvin
    umol2mol = 1e-06,         # conversion umol to mol
    mol2umol = 1e06,          # conversion mol to umol
    kJ2J     = 1000.0,        # conversion kJ to J
    kPa2Pa   = 1000.0         # conversion kPa to Pa
    
  )
}


#################################
### Photosynthetic parameters ###
#################################

photosyn_params <- function(Ci_based){
  
  params1 <- list(
    fl    = 0.15,        # spectral quality correction (Evans 1987)
    Theta = 0.7,         # curvature factor for electron transport rate
    Ox    = 0.21,        # oxygen concentration (mol mol-1)
    frdc3 = 0.011,       # Ratio of dark respiration to PVM at 25degC
    minVc = 1e-12,       # minimum carboxylation rate
    gm_slope = 0.025,
    CcDiffMax = 0.5e-06, # maximum value of "Diff" in photosyn()
    nit_max = 20,        # maximum number of iterations within photosyn calculations
    Ev    = 59.70,       # activation energy for Vcmax (kJ mol-1) (calculated from Kumarathunge et al. 2019 for Tgrowth=15degC)
    Ej    = 40.71,       # activation energy for Jmax (kJ mol-1) (calculated from Kumarathunge et al. 2019 for Tgrowth=15degC and Thome=25egC)
    Er    = 46.39,       # activation energy for dark respiration (kJ mol-1) # Farquhar_1980: 45000
    Hd    = 200,         # deactivation energy (kJ mol-1)
    dSv   = 639.43,      # entropy term for Vcmax (J mol-1) (calculated from Kumarathunge et al. 2019 for Tgrowth=15degC)
    dSj   = 642.97,      # entropy term for Jmax (J mol-1) (calculated from Kumarathunge et al. 2019 for Tgrowth=15degC and Thome=25egC)
    dSr   = 490,         # entropy factor for dark respiration (J mol-1 K-1) Oleson 2013
    bjv   = 1.78,        # ratio Jmax25 / Vcmax25 (-)
    
    ### Acclimation parameters (used in Eq. 2, values as in Table S2)
    a_dSv = 645.13,
    b_dSv = -0.38,
    c_dSv = 0.0,
    d_dSv = 0.0,
    
    a_dSj = 658.77,
    b_dSj = 0.0,
    c_dSj = -0.84,
    d_dSj = -0.52,
    
    a_Haj = 40.71,
    b_Haj = 0.0,
    c_Haj = 0.0,
    d_Haj = 0.0,
    
    a_Hav = 42.6,
    b_Hav = 1.14,
    c_Hav = 0.0,
    d_Hav = 0.0,
    
    a_bjv = 2.56,
    b_bjv = 0.0,     
    c_bjv = -0.0375,
    d_bjv = -0.0202
  ) 
  
  ### implicit gm (Ci-based)
  if (Ci_based){ # Bernacchi et al. 2001, PCE 24 253-259
    
    params2 <-  list(
      Kc0   = 404.9e-06,  # Michaelis-Menten-constant for C02 at 25?C [mol (CO2) mol-1]
      Ko0   = 278.4e-03,  # Michaelis-Menten-constant for 02 at 25?C [mol (O2) mol-1]
      Gam0  = 42.75e-06,  # Gamma* at 25?C [mol(CO2) mol-1]
      Ec    = 79.430,      # activation energy for Kc [J mol-1]
      Eo    = 36.380,      # activation energy for Ko [J mol-1]
      Eg    = 37.830       # activation energy for Gamma* [J mol-1]
    )
    
  ### explicit gm (Cc-based)
  } else { # Bernacchi et al. 2002, Plant Physiology 130 1992-1998
    
    params2 <- list(
      Kc0   = 272.38e-06,  # Michaelis-Menten-constant for C02 at 25degC [mol (CO2) mol-1]
      Ko0   = 165.82e-03,  # Michaelis-Menten-constant for 02 at 25degC [mol (O2) mol-1]
      Gam0  = 37.43e-06,   # Gamma* at 25degC [mol(CO2) mol-1]
      Ec    = 80.99,       # activation energy for Kc [kJ mol-1]
      Eo    = 23.72,       # activation energy for Ko [kJ mol-1]
      Eg    = 24.46,       # activation energy for Gamma* [kJ mol-1]
      
      Em     = 49.6,       # activation energy for gm [kJ mol-1]
      Hdm    = 437.4,      # deactivation energy for gm [kJ mol-1]
      dSm    = 1400,       # entropy term for gm [J mol-1 K-1]
      fgmmin = 0.15        # fraction of minimum gm on gmmax25 [-]
    )
    
  }
  
  params <- c(params1,params2)
  return(params)
  
}



##############################
### Temperature dependency ###
##############################
# k25 = process rate at 25degC
# Tleaf = leaf temperature in degC
# Ha = activation energy in kJ mol-1
# Hd = deactivation energy in kH mol-1
# dS = entropy term (J mol-1 K-1)

## optional parameters (only considered if Temp_acclim = TRUE)
# Tgrowth = growth temperature (degC)
# Thome = home temperature (degC)
# a_Ha, b_Ha, c_Ha, d_Ha = acclimation constants for Ha (see Table S2 in Knauer et al. 2023)
# a_dS, b_dS, c_dS, d_dS = acclimation constants for dS (see Table S2 in Knauer et al. 2023)

# constants = constants as supplied by function photosyn_constants()
Arrhenius_temp_response <- function(k25, Tleaf, Ha, Hd, dS,
                                    Tgrowth = NULL, Thome = NULL,
                                    a_Ha = NULL, b_Ha = NULL, c_Ha = NULL, d_Ha = NULL, 
                                    a_dS = NULL, b_dS = NULL, c_dS = NULL, d_dS = NULL,
                                    constants = photosyn.constants()){
  
  R_gas <- constants$R_gas
  Tref  <- constants$Tref
  
  Tleaf <- Tleaf + constants$Kelvin
  
  # with acclimation
  if (!is.null(Tgrowth)){
    Ha <- (a_Ha + b_Ha * Tgrowth + c_Ha * Thome + d_Ha * (Tgrowth - Thome)) * constants$kJ2J
    dS <- (a_dS + b_dS * Tgrowth + c_dS * Thome + d_dS * (Tgrowth - Thome))
    Hd <- Hd * constants$kJ2J
  } else {
    Ha    <- Ha * constants$kJ2J
    Hd    <- Hd * constants$kJ2J
  }
  
  kTl = k25 * exp(Ha * (Tleaf - Tref)/(Tref * R_gas * Tleaf)) * 
              (1.0 + exp((Tref * dS - Hd)/(Tref * R_gas))) / 
              (1.0 + exp((Tleaf * dS - Hd)/(Tleaf * R_gas)));
  
  return(kTl)
  
}


## Function that converts Ci-based photosynthetic capacity (Vcmax and Jmax) to
## their corresponding Cc-based values
gm.param.conversion <- function(Ci,Vcmax25Ci,Jmax25Ci,gm_max25,Rd,Ox=0.21,
                                Kc_ci=404.9e-06,Ko_ci=278.4e-03,Gam_ci=42.75e-06,
                                Kc_cc=272.38e-06,Ko_cc=165.82e-03,Gam_cc=37.43e-06,
                                Ci_constants=F,algorithm=c("port","default","plinear")){ 
  
  algorithm <- match.arg(algorithm)
  
  ## convert input to the right units
  Ci          <- Ci * 1e-06
  Vcmax25     <- Vcmax25Ci * 1e-06
  J1          <- Jmax25Ci * 1e-06
  Rd          <- Rd * 1e-06
  
  ## 1) calculate photosynthesis (An-Ci response)
  Je <- J1 * (Ci - Gam_ci)/4/(Ci + 2*Gam_ci)
  Jc <- Vcmax25 * (Ci - Gam_ci)/(Ci + Kc_ci * (1 + Ox/Ko_ci))
  
  Ag <- pmin(Je,Jc)
  An <- Ag - Rd
  
  Rubisco_lim <- which(Jc <= Je)
  RuBP_lim    <- which(Je < Jc)
  
  An[An < 0] <- NA
  
  
  ## 2) calculate Cc based on gm
  gm <- gm_max25
  
  Cc <- Ci - An / gm
  
  ### constrain Cc > 0
  if (any(Cc < 0,na.rm=T)){
    warning("Cc falls below 0!!",immediate.=TRUE)
    Cc <- pmax(0,Cc)
  }
  
  
  ## 3) refit Vcmax and Jmax
  
  if (Ci_constants){  # take ci-based Rubisco kinetic constants also for the Cc-based case
    Gam_cc <- Gam_ci
    Kc_cc  <- Kc_ci
    Ko_cc  <- Ko_ci
  }
  
  
  mod <- nls(An ~ c(pmin(J1 * (Cc - Gam_cc)/4/(Cc + 2*Gam_cc),Vcmax * (Cc - Gam_cc)/(Cc + Kc_cc * (1 + Ox/Ko_cc))) - Rd),
             start=c(Vcmax=Vcmax25,J1=J1),
             algorithm=algorithm)
  
  Vcmax25Cc <- summary(mod)$coef[1,1] * 1e06
  Jmax25Cc  <- summary(mod)$coef[2,1] * 1e06
  Ratio     <- Jmax25Cc/Vcmax25Cc   
  
  
  return(list("Vcmax25Cc"=Vcmax25Cc,"Jmax25Cc"=Jmax25Cc,"bjv25"=Ratio))
  
}




#####################################################
### light inhibition of mitochondrial respiration ###
#####################################################
lightinhib_Rd <- function(PPDF){
  inhib <- 0.5 + 0.5 * exp(-1e05 * pmax(PPFD,0)) 
  return(inhib)
}





####################################################################################
# -------------------------------------------------------------------------------- # 
# ----------------------------- Main function ------------------------------------ #
# -------------------------------------------------------------------------------- #
####################################################################################

photosyn <- function(Tair = 25,           # air temperature (degC)
                     Pressure = 101,      # atmospheric pressure (kPa)
                     rH = NULL,           # relative humidity (-); only needed if gs_model == "Ball&Berry"
                     VPD = 1.5,           # vapour pressure deficit (kPa)
                     Ca = 410,            # atmospheric CO2 concentration (umol mol-1)
                     PPFD = 1500,         # photosynthetically active photon flux density (umol m-2 s-1)
                     Vcmax25 = 50.0,      # (Ci-based) maximum carboxylation rate (umol m-2 s-1)
                     Jmax25 = 90.0,       # (Ci-based) maximum electron transport rate (umol m-2 s-1)
                     g0 = 0.005,          # residual stomatal conductance (mol m-2 s-1)
                     g1 = 3.4,            # stomatal slope parameter (kPa^0.5)
                     explicit_gm = TRUE,  # explicit gm considered?
                     temp_acclim = TRUE,  # thermal acclimation considered? 
                     gm25   = 0.2,        # mesophyll conductance at 25 degC
                     beta_s = 1.0,        # soil water stress scalar for stomatal conductance (0-1)
                     beta_m = 1.0,        # soil water stress scalar for mesophyll conductance (0-1)
                     beta_b = 1.0,        # soil water stress scalar for Jmax and Vcmax (0-1)
                     Tgrowth = NULL,      # growth temperature (degC)
                     Thome = NULL,        # home temperature (degC)
                     gs_model=c("USO","Ball&Berry"),  # stomatal conductance model
                     constants=photosyn.constants()){
  
  gs_model  <- match.arg(gs_model)

  params <- photosyn_params(Ci_based=!explicit_gm)
  sapply(1:length(params),function(x) assign(names(params[x]),params[[x]],pos=sys.frame(-3)))
  sapply(1:length(constants),function(x) assign(names(constants[x]),constants[[x]],pos=sys.frame(-4)))
  
  if (temp_acclim & is.null(Tgrowth)){
    stop("argument Tgrowth must be provided if temp_acclim = TRUE!")
  }

  
  ## determine input length
  met <- data.frame(cbind(Tair,Ca,Pressure,PPFD,VPD,rH))
  ts  <- nrow(met)
  sapply(1:ncol(met),function(x) assign(colnames(met)[x],met[,x],pos=sys.frame(-3)))

  if (explicit_gm){
    exp_params <- gm.param.conversion(Ci=seq(0,1500,1),Vcmax25=Vcmax25,Jmax25=Jmax25,
                                      gm_max25=gm25,Rd=frdc3*Vcmax25,Ox=Ox,
                                      Kc_ci=404.9e-06,Ko_ci=278.4e-03,Gam_ci=42.75e-06,
                                      Kc_cc=Kc0,Ko_cc=Ko0,Gam_cc=Gam0,Ci_constants=FALSE)
    Vcmax25_ci <- Vcmax25
    Vcmax25    <- exp_params[["Vcmax25Cc"]]
    Jmax25     <- exp_params[["Jmax25Cc"]]
  } else {
    Vcmax25_ci <- Vcmax25
  }

  
  
  ## unit conversions
  Tleaf      <- Tair
  TleafK     <- Tleaf + Kelvin
  PPFD       <- PPFD       * umol2mol
  Ca         <- Ca         * umol2mol
  Pressure   <- Pressure   * kPa2Pa
  Vcmax25_ci <- Vcmax25_ci * umol2mol
  

  ## create vectors of length ts
  Rubisco_lim = Jc = Je = Cc = gm_gs = Ci_Ca = Cc_Ca = iWUE = An = Ag = gs = Ci = gm <- numeric(length = ts)
  Kc = Ko = Gam = Vcmax = Jmax = APAR_PSII = J1 = Rd <- numeric(length = ts)
  
  ## loop over inputs
  for (i in 1:ts){
    
    Cc[i]  <- 0.5 * Ca[i]
    Diff   <- CcDiffMax * 2
    nit    <- 1
  
    ## Temperature dependencies
    Kc[i]  <- Arrhenius_temp_response(Kc0,Tleaf[i],Ec,0,0)
    Ko[i]  <- Arrhenius_temp_response(Ko0,Tleaf[i],Eo,0,0)
    Gam[i] <- Arrhenius_temp_response(Gam0,Tleaf[i],Eg,0,0)
    
    if (temp_acclim){
      Vcmax[i] <- Arrhenius_temp_response(Vcmax25,Tleaf[i],NULL,Hd,NULL,
                                          Tgrowth=Tgrowth,Thome=Thome,
                                          a_Hav, b_Hav, c_Hav, d_Hav,
                                          a_dSv, b_dSv, c_dSv, d_dSv) 
      Jmax[i] <- Arrhenius_temp_response(Jmax25,Tleaf[i],NULL,Hd,NULL,
                                         Tgrowth=Tgrowth,Thome=Thome,
                                         a_Haj, b_Haj, c_Haj, d_Haj,
                                         a_dSj, b_dSj, c_dSj, d_dSj) 
    } else {
      Vcmax[i] <- Arrhenius_temp_response(Vcmax25,Tleaf[i],Ev,Hd,dSv)
      Jmax[i]  <- Arrhenius_temp_response(Jmax25,Tleaf[i],Ej,Hd,dSj)
    }
    
   
    # apply water stress scalar and convert rates to right unit
    Vcmax[i] <- Vcmax[i] * beta_b * umol2mol
    Jmax[i]  <- Jmax[i] * beta_b * umol2mol
    
    
    ## electron transport rate
    APAR_PSII[i] <- PPFD[i] * (1 - fl) / 2  # absorbed PPFD by photosystem II
    J1[i] <- (APAR_PSII[i] + Jmax[i] - sqrt((APAR_PSII[i] + Jmax[i])^2 - 4 * Theta * APAR_PSII[i] * Jmax[i])) / (2 * Theta)
  
    ## day respiration
    Rd25  <- frdc3 * Vcmax25_ci
    Rd[i] <- Arrhenius_temp_response(Rd25,Tleaf[i],Er,0.0,0.0) * lightinhib_Rd(PPFD[i])
  
    ## mesophyll conductance
    if (!explicit_gm){ # implicit gm
      gm[i] <- 9e30
    } else {  # explicit gm
      gm[i] <- Arrhenius_temp_response(gm25,Tleaf[i],Em,Hdm,dSm) * beta_m 
      gm[i] <- max(c(fgmmin * gm25), gm[i])
    }
  
    
    ## calculate Photosynthesis, stomatal conductance, intercellular and chloroplastic CO2 concentration
    ## and others for each input step. Find a solution for those through iteration.

    while (abs(Diff) > CcDiffMax & !is.na(Diff) & nit < nit_max){
      
      Je[i] <- J1[i] * (Cc[i] - Gam[i])/4/(Cc[i] + 2*Gam[i])
      Jc[i] <- Vcmax[i] * (Cc[i] - Gam[i])/(Cc[i] + Kc[i] * (1 + Ox/Ko[i]))
      
      Rubisco_lim[i] <- ifelse(Jc[i] <= Je[i],TRUE,FALSE)      
      Ag[i] <- pmin(Je[i],Jc[i]) * 1
      An[i] <- Ag[i] - Rd[i]
      
      if (gs_model == "Ball&Berry"){
        gs[i] <- max(g0*(R_gas*TleafK[i]/Pressure[i]),
                     ((g0 + g1 * beta_s * (An[i] * rH[i])/Ca[i]) * R_gas*TleafK[i]/Pressure[i]))  # gs in ms-1 for water
      } else if (gs_model == "USO"){
        gs[i] <- max(g0*(R_gas*TleafK[i]/Pressure[i]),
                     ((g0 + 1.6 * (1.0 + (g1 * beta_s)/sqrt(max(VPD[i],0.1))) * (An[i]/Ca[i])) * R_gas*TleafK[i]/Pressure[i]))  # gs in ms-1 for water
      }
      
      Ci[i] <- Ca[i] - An[i]/((gs[i]/1.6)*Pressure[i]/(R_gas*TleafK[i]))   # gs in mol m-2 s-1 for CO2
      
      Cc_new  <- Ci[i] - An[i]/gm[i]
      
      Diff  <- Cc[i] - Cc_new
      Cc[i] <- Cc[i] - Diff/2
      
      iWUE[i]  <- An[i]/((gs[i]/1.6)*Pressure[i]/(R_gas*TleafK[i]))
      gm_gs[i] <- gm[i]/((gs[i]/1.6)*Pressure[i]/(R_gas*TleafK[i]))
      Ci_Ca[i] <- Ci[i]/Ca[i]
      Cc_Ca[i] <- Cc[i]/Ca[i]
      
      nit <- nit + 1
    } 
    Cc[i] <- Cc_new # ensure that Cc of output is the one taken for the calculations
  } # end time loop  

  out <- data.frame(Ag,An,Jc,Je,Rd,Vcmax,Jmax,J1,Gam,Kc,Ko,gs,gm,Ci,Cc,iWUE,
                    gm_gs,Ci_Ca,Cc_Ca,Rubisco_lim)
  non_convert <- c("gm_gs","gm","gs","Ci_Ca","Cc_Ca","Rubisco_lim")
  out[,!(colnames(out) %in% non_convert)] <- out[,!(colnames(out) %in% non_convert)] * mol2umol
  
  return(out)
  
}
