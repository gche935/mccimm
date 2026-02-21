###############################################################################################
# Author: Prof. Gordon Cheung (University of Auckland)                                        #
# email: gordon.cheung@auckland.ac.nz                                                         #
# Purpose: Monte Carlo simulation of conditional indirect effects from modsem & Mplus outputs #
# Model: Allows 3-way interaction                                                             #
# Allows simulated function (Sfunction = "    "                                               #
# Version: Beta 1.3.1                                                                         #
# WARNING: This is a beta version. Please report all bugs to the author                       #
###############################################################################################

## == Required packages (load when loading mccimm) == ##
#library(MASS)
#library(ggplot2)
#library(tidyr)
#library(MplusAutomation)


## ====== Function "mccimm_modsem" Monte Carlo Simulation for Confidence Intervals of Moderated Mediation (modsem) ====== ##
#' Monte Carlo Simulation for Confidence Intervals of Moderated Mediation (modsem)
#'
#' Generate confidence intervals of moderated-mediating effects from modsem results using Monte Carlo simulation.
#'
#' \if{html}{
#' \figure{Figure.png}{options: width="75\%" alt="Description of my figure"}
#' }
#' \if{latex}{
#' \figure{Figure.pdf}{options: width=15cm}
#' }
#'
#' @param object modsem object (output from modsem).
#' @param Z name of moderate Z.
#' @param W name of moderator W.
#' @param a1 parameter name of a1 path (main effect).
#' @param a2 parameter name of a2 path (main effect).
#' @param a3 parameter name of a3 path (main effect).
#' @param a4 parameter name of a4 path (main effect).
#' @param z1 parameter name of z1 path (interaction effect).
#' @param z2 parameter name of z2 path (interaction effect).
#' @param z3 parameter name of z3 path (interaction effect).
#' @param z4 parameter name of z4 path (interaction effect).
#' @param w1 parameter name of w1 path (interaction effect).
#' @param w2 parameter name of w2 path (interaction effect).
#' @param w3 parameter name of w3 path (interaction effect).
#' @param w4 parameter name of w4 path (interaction effect).
#' @param zw1 parameter name of zw1 path (3-way interaction effect).
#' @param zw2 parameter name of zw2 path (3-way interaction effect).
#' @param zw3 parameter name of zw3 path (3-way interaction effect).
#' @param zw4 parameter name of zw4 path (3-way interaction effect).
#' @param R number of Monte Carlo simulation samples (in millions). For example, R=5 (default) generates 5,000,000 simulated samples.
#'
#' @return mccimm output for plotting Johnson-Neyman Figure.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' # modsem object is "est_lms" & output mccimm object is mcObject
#'
#' mcObject <- mccimm_modsem(est_lms, a1="a1", a2="a2", a3="a3a", a1="z1", Z="Autonomy")
#'
#'
#' # Change 2-Way Figure Title and/or Axis Labels Afterwards
#' p_int <- p_int + ggplot2::labs(title = "Replace with your Figure Title",
#'                                   x = "Replace with your X-axis Label",
#'                                   y = "Replace with your Y-axis Label")
#'
#' # Change Figure Legend Labels Afterwards
#' p_int <- p_int + ggplot2::scale_color_manual(name = "Replace with your Legend Title",
#'                                           values = c("line1" = "black", "line2" = "grey"),
#'                                           labels = c("Z at mean - 1sd", "Z at mean + 1sd"))
#'
#' # Save the New Figure
#' ggplot2::ggsave("New Standardized Interaction Figure.png", width = 22.86, height = 16.51, units = "cm")
#'

  mccimm_modsem <- function(object, Z="NA", W="NA",
                   a1="NA", z1="NA", w1="NA", zw1="NA",
                   a2="NA", z2="NA", w2="NA", zw2="NA",
                   a3="NA", z3="NA", w3="NA", zw3="NA",
                   a4="NA", z4="NA", w4="NA", zw4="NA",
                   R=5) {

  ## --- Initial Inputs for programming --- ##

#  object <- est_lms
#  Z <- "Autonomy"          ## Specify moderator Z
#  W <- "NA"
#  a1 <- "a1"
#  a2 <- "a2"
#  a3 <- "a3a"
#  a4 <- "NA"
#  z1 <- "z1"
#  z2 <- "NA"
#  z3 <- "NA"
#  z4 <- "NA"
#  w1 <- "NA"
#  w2 <- "NA"
#  w3 <- "NA"
#  w4 <- "NA"
#  zw1 <- "NA"
#  zw2 <- "NA"
#  zw3 <- "NA"
#  zw4 <- "NA"

#  R <- 1 # Number of simulated samples = R*1e6 (default: R = 5)
  ## ------------------------------- ##


  ## -- Extract defined parameters and vcov -- ##
  varZ <- "NA"
  varW <- "NA"
  if (Z != "NA") varZ <- paste0(Z, "~~", Z)
  if (W != "NA") varW <- paste0(W, "~~", W)
  dp <- c(a1, a2, a3, a4, z1, z2, z3, z4, w1, w2, w3, w4, zw1, zw2, zw3, zw4, varZ, varW)
  dp <- dp[dp != "NA"]


  temp <- modsem_coef(object)
  estcoeff <- temp[dp]
  Temp3 <- modsem_vcov(object)
  Tech3 <- Temp3[dp, dp]

  dd <- modsem::standardized_estimates(object, correction=TRUE)
  stdyx.temp <- dd[, "est"]
  names(stdyx.temp) <- names(temp)
  stdyx.estcoeff <- stdyx.temp[dp]

#  stdyx.temp <- modsem_coef(object, standardized=TRUE)
#  stdyx.estcoeff <- stdyx.temp[dp]

  return_mccimm <- mccimm(estcoeff, stdyx.estcoeff, Tech3,
                        Z, W,
                        varZ, varW,
                        a1, z1, w1, zw1,
                        a2, z2, w2, zw2,
                        a3, z3, w3, zw3,
                        a4, z4, w4, zw4,
                        R=5)

  return(return_mccimm)

}  ## end (Function "mccimm_modsem") ##




## ====== Function "mccimm_lavaan" Monte Carlo Simulation for Confidence Intervals of Moderated Mediation (lavaan) ====== ##
#' Monte Carlo Simulation for Confidence Intervals of Moderated Mediation (lavaan)
#'
#' Generate confidence intervals of moderated-mediating effects from lavaan results using Monte Carlo simulation.
#'
#' \if{html}{
#' \figure{Figure.png}{options: width="75\%" alt="Description of my figure"}
#' }
#' \if{latex}{
#' \figure{Figure.pdf}{options: width=15cm}
#' }
#'
#' @param object lavaan object (output from lavaan).
#' @param Z name of moderate Z.
#' @param W name of moderator W.
#' @param a1 parameter name of a1 path (main effect).
#' @param a2 parameter name of a2 path (main effect).
#' @param a3 parameter name of a3 path (main effect).
#' @param a4 parameter name of a4 path (main effect).
#' @param z1 parameter name of z1 path (interaction effect).
#' @param z2 parameter name of z2 path (interaction effect).
#' @param z3 parameter name of z3 path (interaction effect).
#' @param z4 parameter name of z4 path (interaction effect).
#' @param w1 parameter name of w1 path (interaction effect).
#' @param w2 parameter name of w2 path (interaction effect).
#' @param w3 parameter name of w3 path (interaction effect).
#' @param w4 parameter name of w4 path (interaction effect).
#' @param zw1 parameter name of zw1 path (3-way interaction effect).
#' @param zw2 parameter name of zw2 path (3-way interaction effect).
#' @param zw3 parameter name of zw3 path (3-way interaction effect).
#' @param zw4 parameter name of zw4 path (3-way interaction effect).
#' @param R number of Monte Carlo simulation samples (in millions). For example, R=5 (default) generates 5,000,000 simulated samples.
#'
#' @return mccimm output for plotting Johnson-Neyman Figure.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' # lavaan object is "est_lms" & output mccimm object is mcObject
#'
#' mcObject <- mccimm_lavaan(est_lms, a1="a1", a2="a2", a3="a3a", a1="z1", Z="Autonomy")
#'
#'
#' # Change 2-Way Figure Title and/or Axis Labels Afterwards
#' p_int <- p_int + ggplot2::labs(title = "Replace with your Figure Title",
#'                                   x = "Replace with your X-axis Label",
#'                                   y = "Replace with your Y-axis Label")
#'
#' # Change Figure Legend Labels Afterwards
#' p_int <- p_int + ggplot2::scale_color_manual(name = "Replace with your Legend Title",
#'                                           values = c("line1" = "black", "line2" = "grey"),
#'                                           labels = c("Z at mean - 1sd", "Z at mean + 1sd"))
#'
#' # Save the New Figure
#' ggplot2::ggsave("New Standardized Interaction Figure.png", width = 22.86, height = 16.51, units = "cm")
#'

  mccimm_lavaan <- function(object, Z="NA", W="NA",
                   a1="NA", z1="NA", w1="NA", zw1="NA",
                   a2="NA", z2="NA", w2="NA", zw2="NA",
                   a3="NA", z3="NA", w3="NA", zw3="NA",
                   a4="NA", z4="NA", w4="NA", zw4="NA",
                   R=5) {

  ## --- Initial Inputs for programming --- ##

#  R <- 1 # Number of simulated samples = R*1e6 (default: R = 5)
  ## ------------------------------- ##


  ## -- Extract defined parameters and vcov -- ##
  varZ <- "NA"
  varW <- "NA"
  if (Z != "NA") varZ <- paste0(Z, "~~", Z)
  if (W != "NA") varW <- paste0(W, "~~", W)
  dp <- c(a1, a2, a3, a4, z1, z2, z3, z4, w1, w2, w3, w4, zw1, zw2, zw3, zw4, varZ, varW)
  dp <- dp[dp != "NA"]


  temp <- lavaan::coef(object)
  estcoeff <- temp[dp]
  Temp3 <- lavaan::vcov(object)
  Tech3 <- Temp3[dp, dp]

  dd <- lavaan::parameterEstimates(object, standardized=TRUE, remove.nonfree=TRUE, remove.def=TRUE)
  stdyx.temp <- dd[, "std.all"]
  names(stdyx.temp) <- names(temp)
  stdyx.estcoeff <- stdyx.temp[dp]

  return_mccimm <- mccimm(estcoeff, stdyx.estcoeff, Tech3,
                        Z, W,
                        varZ, varW,
                        a1, z1, w1, zw1,
                        a2, z2, w2, zw2,
                        a3, z3, w3, zw3,
                        a4, z4, w4, zw4,
                        R=5)

  return(return_mccimm)

}  ## end (Function "mccimm_lavaan") ##





## ====== Function "mccimm_mplus" Monte Carlo Confidence Intervals for Moderated Mediation (mplus) ====== ##
#' Monte Carlo Simulation for Confidence Intervals of Moderated Mediation (mplus)
#'
#' Generate confidence intervals of moderated-mediating effects from modsem results using Monte Carlo simulation.
#' Location of estimated parameters can be found in mccimm::TECH1().
#'
#' \if{html}{
#' \figure{Figure.png}{options: width="75\%" alt="Description of my figure"}
#' }
#' \if{latex}{
#' \figure{Figure.pdf}{options: width=15cm}
#' }
#'
#' @param mplus_output_file Mplus output (.out) file (output from Mplus).
#' @param results_file Mplus a text file (.txt) that saves the Mplus results (RESULTS IS "filename.txt" in Mplus SAVEDATA:).
#' @param Z name of moderate Z.
#' @param W name of moderator W.
#' @param varZ location of parameter varZ in Mplus Tech1 outputs.
#' @param varW location of parameter varW in Mplus Tech1 outputs.
#' @param a1 location of parameter a1 in Mplus Tech1 outputs.
#' @param a2 location of parameter a2 in Mplus Tech1 outputs.
#' @param a3 location of parameter a3 in Mplus Tech1 outputs.
#' @param a4 location of parameter a4 in Mplus Tech1 outputs.
#' @param z1 location of parameter z1 in Mplus Tech1 outputs.
#' @param z2 location of parameter z2 in Mplus Tech1 outputs.
#' @param z3 location of parameter z3 in Mplus Tech1 outputs.
#' @param z4 location of parameter z4 in Mplus Tech1 outputs.
#' @param w1 location of parameter w1 in Mplus Tech1 outputs.
#' @param w2 location of parameter w2 in Mplus Tech1 outputs.
#' @param w3 location of parameter w3 in Mplus Tech1 outputs.
#' @param w4 location of parameter w4 in Mplus Tech1 outputs.
#' @param zw1 location of parameter zw1 in Mplus Tech1 outputs.
#' @param zw2 location of parameter zw2 in Mplus Tech1 outputs.
#' @param zw3 location of parameter zw3 in Mplus Tech1 outputs.
#' @param zw4 location of parameter zw4 in Mplus Tech1 outputs.
#' @param R number of Monte Carlo simulation samples (in millions). For example, R=5 (default) generates 5,000,000 simulated samples.
#'
#' @return mccimm output for plotting Johnson-Neyman Figure.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' # mplus_output_file is "model cc4.out", results_file is "Model_CC4.txt" & moderator Z is "AUTO"
#' # output mccimm object is mcObject
#'
#' mcObject <- mccimm_mplus("model cc4.out", "Model_CC4.txt", Z = "AUTO", varZ = "72",
#'             a1 = "60", a2 = "65", z1 = "62")#'
#'
#' # Change 2-Way Figure Title and/or Axis Labels Afterwards
#' p_int <- p_int + ggplot2::labs(title = "Replace with your Figure Title",
#'                                   x = "Replace with your X-axis Label",
#'                                   y = "Replace with your Y-axis Label")
#'
#' # Change Figure Legend Labels Afterwards
#' p_int <- p_int + ggplot2::scale_color_manual(name = "Replace with your Legend Title",
#'                                           values = c("line1" = "black", "line2" = "grey"),
#'                                           labels = c("Z at mean - 1sd", "Z at mean + 1sd"))
#'
#' # Save the New Figure
#' ggplot2::ggsave("New Standardized Interaction Figure.png", width = 22.86, height = 16.51, units = "cm")
#'

mccimm_mplus <- function(mplus_output_file = "mplus_output.out",
                   results_file = "results.txt",
                   Z="NA", W="NA",
                   varZ="NA", varW="NA",
                   a1="NA", a2="NA", a3="NA", a4="NA",
                   z1="NA", z2="NA", z3="NA", z4="NA",
                   w1="NA", w2="NA", w3="NA", w4="NA",
                   zw1="NA", zw2="NA", zw3="NA", zw4="NA",
                   R=5) {

  mplus_output <- readModels(mplus_output_file)
  results <- mplus_output$parameters$unstandardized
  temp <- scan(results_file, sep="")
  stdyx.results <- mplus_output$parameters$stdyx.standardized
  no.parameters <- mplus_output$summaries$Parameters
  stdyx.temp <- temp[-(1:(2*no.parameters))]

  Temp3 <- mplus_output$tech3$paramCov
  Temp3[upper.tri(Temp3, diag = FALSE)] <- 0

  Tech3 <- Temp3 + t(Temp3)
  Tech3 <- Tech3 - diag(diag(Temp3))

  ## -- Extract defined parameters and vcov -- ##
  dp <- c(varZ, varW, a1, a2, a3, a4, z1, z2, z3, z4, w1, w2, w3, w4, zw1, zw2, zw3, zw4)

  dp_no <- suppressWarnings(as.numeric(dp))
  dp_list <- c("varZ", "varW", "a1", "a2", "a3", "a4", "z1", "z2", "z3", "z4", "w1", "w2", "w3", "w4", "zw1", "zw2", "zw3", "zw4")
  non_na_list <- dp_list[which(!is.na(dp_no))]
  dp <- dp[dp != "NA"]
  dp <- as.numeric(dp)

  estcoeff <- temp[dp]
  names(estcoeff) <- non_na_list
  stdyx.estcoeff <- stdyx.temp[dp]
  names(stdyx.estcoeff) <- non_na_list

  Tech3 <- Tech3[dp, dp]
  rownames(Tech3) <- non_na_list
  colnames(Tech3) <- non_na_list

  # -- Reassigning variable names -- #
  if (varZ != "NA") varZ <- "varZ"
  if (varW != "NA") varW <- "varW"
  if (a1 != "NA") a1 <- "a1"
  if (a2 != "NA") a2 <- "a2"
  if (a3 != "NA") a3 <- "a3"
  if (a4 != "NA") a4 <- "a4"
  if (z1 != "NA") z1 <- "z1"
  if (z2 != "NA") z2 <- "z2"
  if (z3 != "NA") z3 <- "z3"
  if (z4 != "NA") z4 <- "z4"
  if (w1 != "NA") w1 <- "w1"
  if (w2 != "NA") w2 <- "w2"
  if (w3 != "NA") w3 <- "w3"
  if (w4 != "NA") w4 <- "w4"
  if (zw1 != "NA") zw1 <- "zw1"
  if (zw2 != "NA") zw2 <- "zw2"
  if (zw3 != "NA") zw3 <- "zw3"
  if (zw4 != "NA") zw4 <- "zw4"

  return_mccimm <- mccimm(estcoeff, stdyx.estcoeff, Tech3,
                        Z, W,
                        varZ, varW,
                        a1, z1, w1, zw1,
                        a2, z2, w2, zw2,
                        a3, z3, w3, zw3,
                        a4, z4, w4, zw4,
                        R=5)

  return(return_mccimm)


  ## ------------------------------- ##

}  ## end (Function "mccimm_mplus") ##



## ====== Sub-Function "mccimm" Monte Carlo Confidence Intervals for Moderated Mediation ====== ##

mccimm <- function(estcoeff, stdyx.estcoeff, Tech3,
                   Z="NA", W="NA",
                   varZ="NA", varW="NA",
                   a1="NA", z1="NA", w1="NA", zw1="NA",
                   a2="NA", z2="NA", w2="NA", zw2="NA",
                   a3="NA", z3="NA", w3="NA", zw3="NA",
                   a4="NA", z4="NA", w4="NA", zw4="NA",
                   R=5) {

  ## -- Number of Moderating Effects -- ##
  NoModz <- 0
  NoModw <- 0
  if (z1 != "NA") NoModz <- NoModz + 1
  if (z2 != "NA") NoModz <- NoModz + 1
  if (z3 != "NA") NoModz <- NoModz + 1
  if (z4 != "NA") NoModz <- NoModz + 1
  if (w1 != "NA") NoModw <- NoModw + 1
  if (w2 != "NA") NoModw <- NoModw + 1
  if (w3 != "NA") NoModw <- NoModw + 1
  if (w4 != "NA") NoModw <- NoModw + 1
  NoMod <- NoModz + NoModw
  ## ----- ##


  ## -- Location of Moderator for calculation of Index MM -- ##
  if (NoMod == 1) {
    if (z1 != "NA") PoMod <- 1
    if (z2 != "NA") PoMod <- 2
    if (z3 != "NA") PoMod <- 3
    if (z4 != "NA") PoMod <- 4
    if (w1 != "NA") PoMod <- 1
    if (w2 != "NA") PoMod <- 2
    if (w3 != "NA") PoMod <- 3
    if (w4 != "NA") PoMod <- 4
  } # end (if (NoMod == 1))
  ## ------------------------------ ##


  ## -- Initialize a-paths to 1 -- ##
  U7Xa1 <- 1
  U7Xa2 <- 1
  U7Xa3 <- 1
  U7Xa4 <- 1
  U7XSa1 <- matrix(1, R*1e6)
  U7XSa2 <- matrix(1, R*1e6)
  U7XSa3 <- matrix(1, R*1e6)
  U7XSa4 <- matrix(1, R*1e6)
  ## ------------------------------ ##


  ## -- Initialize interaction paths to 0 -- ##
  U7Xz1 <- 0
  U7Xz2 <- 0
  U7Xz3 <- 0
  U7Xz4 <- 0
  U7XSz1 <- matrix(0, R*1e6)
  U7XSz2 <- matrix(0, R*1e6)
  U7XSz3 <- matrix(0, R*1e6)
  U7XSz4 <- matrix(0, R*1e6)

  U7Xw1 <- 0
  U7Xw2 <- 0
  U7Xw3 <- 0
  U7Xw4 <- 0
  U7XSw1 <- matrix(0, R*1e6)
  U7XSw2 <- matrix(0, R*1e6)
  U7XSw3 <- matrix(0, R*1e6)
  U7XSw4 <- matrix(0, R*1e6)

  U7Xzw1 <- 0
  U7Xzw2 <- 0
  U7Xzw3 <- 0
  U7Xzw4 <- 0
  U7XSzw1 <- matrix(0, R*1e6)
  U7XSzw2 <- matrix(0, R*1e6)
  U7XSzw3 <- matrix(0, R*1e6)
  U7XSzw4 <- matrix(0, R*1e6)
  ## ------------------------------ ##


  ## -- Initialize stdZ and stdW -- ##
  stdZ <- 1
  stdW <- 1
  SstdZ <- matrix(1, R*1e6)
  SstdW <- matrix(1, R*1e6)
  ## ------------------------------ ##


  ## -- Initialize standardized a-paths to 1 -- ##
  Z7Xa1 <- 1
  Z7Xa2 <- 1
  Z7Xa3 <- 1
  Z7Xa4 <- 1
  ## ------------------------------ ##


  ## -- Initialize standardized interaction paths to 0 -- ##
  Z7Xz1 <- 0
  Z7Xz2 <- 0
  Z7Xz3 <- 0
  Z7Xz4 <- 0

  Z7Xw1 <- 0
  Z7Xw2 <- 0
  Z7Xw3 <- 0
  Z7Xw4 <- 0

  Z7Xzw1 <- 0
  Z7Xzw2 <- 0
  Z7Xzw3 <- 0
  Z7Xzw4 <- 0
  ## ------------------------------ ##


  ## -- Monte Carlo Simulation of R*1e6 samples, default: R = 5 -- ##
  mcmc <- MASS::mvrnorm(n=R*1e6, mu=estcoeff, Sigma=Tech3, tol = 1e-6)

  # -- Retain simulated samples with variance larger than or equal to 0 -- #
  if (NoMod == 1) {
    if (NoModz != 0){
      mcmc <- mcmc[which(mcmc[, varZ] >= 0),]
    } else {
      mcmc <- mcmc[which(mcmc[, varW] >= 0),]
    }
  } else if (NoMod == 2) {
    mcmc <- mcmc[which(mcmc[, varZ] >= 0),]
    mcmc <- mcmc[which(mcmc[, varW] >= 0),]
  } # end (if NoMod)

  b.no <- nrow(mcmc)
  R.no <- format(R*1e6, scientific = FALSE)

  # -- Print number of bootstrap samples -- #
  cat("\n", "   Number of requested simulated samples = ", R.no)
  cat("\n", "   Number of completed simulated samples = ", b.no, rep("\n",2))

  ## ------------------------------ ##


  ### --- No Moderating Effect --- ###
  if (NoMod == 0) {

    # -- Define estimated parameters for calculating indirect effects -- #
    if (any(names(estcoeff) %in% a1)) U7Xa1 <- estcoeff[a1]
    if (any(names(estcoeff) %in% a2)) U7Xa2 <- estcoeff[a2]
    if (any(names(estcoeff) %in% a3)) U7Xa3 <- estcoeff[a3]
    if (any(names(estcoeff) %in% a4)) U7Xa4 <- estcoeff[a4]

    # -- Calculate Estimated Indirect Effect -- #
    estM  <- U7Xa1*U7Xa2*U7Xa3*U7Xa4

    # -- Capture simulated parameters for calculating indirect effects -- #
    if (any(names(estcoeff) %in% a1)) U7XSa1 <- mcmc[, a1]
    if (any(names(estcoeff) %in% a2)) U7XSa2 <- mcmc[, a2]
    if (any(names(estcoeff) %in% a3)) U7XSa3 <- mcmc[, a3]
    if (any(names(estcoeff) %in% a4)) U7XSa4 <- mcmc[, a4]

    # -- Calculate Simulated Indirect Effect -- #
    abM <- U7XSa1*U7XSa2*U7XSa3*U7XSa4

    # == Confidence Intervals and p-value == #

    # -- Calculate Percentile Probability -- #
    if (quantile(abM,probs=0.5)>0) {
      pM = 2*(sum(abM<0)/b.no)
    } else {
      pM = 2*(sum(abM>0)/b.no)
    }


    #### Percentile Confidence Intervals of Conditional Indirect Effects ####

    PCI <- matrix(1:8, nrow = 1, dimnames = list(c("        "),
                                             c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    PCI[1,1] <- format(round(quantile(abM,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[1,2] <- format(round(quantile(abM,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[1,3] <- format(round(quantile(abM,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[1,4] <- format(round(estM, digits = 4), nsmall = 4, scientific = FALSE)
    PCI[1,5] <- format(round(quantile(abM,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[1,6] <- format(round(quantile(abM,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[1,7] <- format(round(quantile(abM,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[1,8] <- format(round(pM, digits = 4), nsmall = 4, scientific = FALSE)

    # Bias-Corrected Factor
    zM = qnorm(sum(abM<estM)/b.no)

    # Calculate Bias-Corrected Probability

    if ((estM>0 & min(abM)>0) | (estM<0 & max(abM)<0)) {
      pbM = 0
    } else if (qnorm(sum(abM>0)/b.no)+2*zM<0) {
      pbM = 2*pnorm((qnorm(sum(abM>0)/b.no)+2*zM))
    } else {
      pbM = 2*pnorm(-1*(qnorm(sum(abM>0)/b.no)+2*zM))
    }

    #### Bias-Corrected Confidence Intervals ####
    BCCI <- matrix(1:8, nrow = 1, dimnames = list(c("        "),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    BCCI[1,1] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[1,2] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[1,3] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[1,4] <- format(round(estM, digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[1,5] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[1,6] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[1,7] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[1,8] <- format(round(pbM, digits = 4), nsmall = 4, scientific = FALSE)

    cat("\n")
    cat("Percentile Confidence Intervals for Indirect Effect", rep("\n", 2))
    rownames(PCI) <- rep("    ", nrow(PCI))
    print(PCI, quote=FALSE, right=TRUE)
    cat("\n")

    cat("\n")
    cat("Bias-Corrected Confidence Intervals for Indirect Effect", rep("\n", 2))
    rownames(BCCI) <- rep("    ", nrow(BCCI))
    print(BCCI, quote=FALSE, right=TRUE)
    cat("\n")


    # - Define estimated parameters for calculating standardized indirect effects - #
    if (any(names(stdyx.estcoeff) %in% a1)) Z7Xa1 <- stdyx.estcoeff[a1]
    if (any(names(stdyx.estcoeff) %in% a2)) Z7Xa2 <- stdyx.estcoeff[a2]
    if (any(names(stdyx.estcoeff) %in% a3)) Z7Xa3 <- stdyx.estcoeff[a3]
    if (any(names(stdyx.estcoeff) %in% a4)) Z7Xa4 <- stdyx.estcoeff[a4]

    # - Calculate Estimated Indirect Effect - #
    estM  <- Z7Xa1*Z7Xa2*Z7Xa3*Z7Xa4

    # -- Print standardized indirect effects -- #
    cat("\n", "   Standardized indirect effects = ", estM, rep("\n",2))

  } # end (if (NoMod == 0))
  ### --- End No Moderating Effect --- ###


  ### --- Only One Moderating Effect Z --- ###
  if ((NoMod == 1) & (NoModz == 1)) {

    # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% a1)) U7Xa1 <- estcoeff[a1]
    if (any(names(estcoeff) %in% a2)) U7Xa2 <- estcoeff[a2]
    if (any(names(estcoeff) %in% a3)) U7Xa3 <- estcoeff[a3]
    if (any(names(estcoeff) %in% a4)) U7Xa4 <- estcoeff[a4]
    if (any(names(estcoeff) %in% z1)) U7Xz1 <- estcoeff[z1]
    if (any(names(estcoeff) %in% z2)) U7Xz2 <- estcoeff[z2]
    if (any(names(estcoeff) %in% z3)) U7Xz3 <- estcoeff[z3]
    if (any(names(estcoeff) %in% z4)) U7Xz4 <- estcoeff[z4]
    stdZ <- sqrt(estcoeff[varZ])

    # Calculate Estimated Index MM #
    if (PoMod == 1) estIMM <- U7Xz1*U7Xa2*U7Xa3*U7Xa4
    if (PoMod == 2) estIMM <- U7Xa1*U7Xz2*U7Xa3*U7Xa4
    if (PoMod == 3) estIMM <- U7Xa1*U7Xa2*U7Xz3*U7Xa4
    if (PoMod == 4) estIMM <- U7Xa1*U7Xa2*U7Xa3*U7Xz4

    # Capture simulated parameters for calculating indirect effects #
    if (any(names(estcoeff) %in% a1)) U7XSa1 <- mcmc[, a1]
    if (any(names(estcoeff) %in% a2)) U7XSa2 <- mcmc[, a2]
    if (any(names(estcoeff) %in% a3)) U7XSa3 <- mcmc[, a3]
    if (any(names(estcoeff) %in% a4)) U7XSa4 <- mcmc[, a4]
    if (any(names(estcoeff) %in% z1)) U7XSz1 <- mcmc[, z1]
    if (any(names(estcoeff) %in% z2)) U7XSz2 <- mcmc[, z2]
    if (any(names(estcoeff) %in% z3)) U7XSz3 <- mcmc[, z3]
    if (any(names(estcoeff) %in% z4)) U7XSz4 <- mcmc[, z4]
    SstdZ <- sqrt(mcmc[,varZ])

    # Calculate Simulated Index MM
    if (any(names(estcoeff) %in% z1)) IndexMM <- U7XSz1*U7XSa2*U7XSa3*U7XSa4
    if (any(names(estcoeff) %in% z2)) IndexMM <- U7XSa1*U7XSz2*U7XSa3*U7XSa4
    if (any(names(estcoeff) %in% z3)) IndexMM <- U7XSa1*U7XSa2*U7XSz3*U7XSa4
    if (any(names(estcoeff) %in% z4)) IndexMM <- U7XSa1*U7XSa2*U7XSa3*U7XSz4

    #### Percentile and Bias-Corrected Confidence Intervals of Conditional Indirect Effects ####

    PCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))
    BCCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                               c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    for (X in 1:5) {
      level <- X - 3
      estX <- (U7Xa1+U7Xz1*level*stdZ)*(U7Xa2+U7Xz2*level*stdZ)*(U7Xa3+U7Xz3*level*stdZ)*(U7Xa4+U7Xz4*level*stdZ)
      abX <- (U7XSa1+U7XSz1*level*SstdZ)*(U7XSa2+U7XSz2*level*SstdZ)*(U7XSa3+U7XSz3*level*SstdZ)*(U7XSa4+U7XSz4*level*SstdZ)
      zX = qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

      ## Percentile Confidence Intervals ##
      PCI[X,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

      # Percentile p-value #
      if (quantile(abX,probs=0.5)>0) {
        PCI[X,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
      } else {
        PCI[X,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
      }

      ## Bias-Corrected Confidence Intervals ##
      BCCI[X,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

      # Bias-Corrected p-value
      if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
        BCCI[X,8] = 0
      } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
        BCCI[X,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
      } else {
        BCCI[X,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
      }
    }

    # Percentile Confidence Intervals for Index MM
    PCI[6,1] <- format(round(quantile(IndexMM,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,2] <- format(round(quantile(IndexMM,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,3] <- format(round(quantile(IndexMM,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,4] <- format(round(estIMM, digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,5] <- format(round(quantile(IndexMM,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,6] <- format(round(quantile(IndexMM,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,7] <- format(round(quantile(IndexMM,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

    # Percentile Probability for Index MM
    if (quantile(IndexMM,probs=0.5)>0) {
      PCI[6,8] <- format(round(2*(sum(IndexMM<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)  # Probability
    } else {
      PCI[6,8] <- format(round(2*(sum(IndexMM>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)  # Probability
    }

    # Bias-Corrected Confidence Intervals for Index MM
    zIMM = qnorm(sum(IndexMM<estIMM)/b.no)
    BCCI[6,1] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,2] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,3] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,4] <- format(round(estIMM, digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,5] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,6] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,7] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

    # Bias-Corrected Probability for Index MM
    if ((estIMM>0 & min(IndexMM)>0) | (estIMM<0 & max(IndexMM)<0)) {
      BCCI[6,8] = 0
    } else if (qnorm(sum(IndexMM>0)/b.no)+2*zIMM<0) {
      BCCI[6,8] = format(round(2*pnorm((qnorm(sum(IndexMM>0)/b.no)+2*zIMM)), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      BCCI[6,8] = format(round(2*pnorm(-1*(qnorm(sum(IndexMM>0)/b.no)+2*zIMM)), digits = 4), nsmall = 4, scientific = FALSE)
    }

    cat("\n")
    cat("Percentile Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(PCI, quote=FALSE, right=TRUE)
    cat("\n")

    cat("\n")
    cat("Bias-Corrected Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(BCCI, quote=FALSE, right=TRUE)
    cat("\n")

    # -- Plot Standardized 2-Way Interaction Effects -- #
    # Define estimated parameters for calculating indirect effects
    if (any(names(stdyx.estcoeff) %in% a1)) Z7Xa1 <- stdyx.estcoeff[a1]
    if (any(names(stdyx.estcoeff) %in% a2)) Z7Xa2 <- stdyx.estcoeff[a2]
    if (any(names(stdyx.estcoeff) %in% a3)) Z7Xa3 <- stdyx.estcoeff[a3]
    if (any(names(stdyx.estcoeff) %in% a4)) Z7Xa4 <- stdyx.estcoeff[a4]
    if (any(names(stdyx.estcoeff) %in% z1)) Z7Xz1 <- stdyx.estcoeff[z1]
    if (any(names(stdyx.estcoeff) %in% z2)) Z7Xz2 <- stdyx.estcoeff[z2]
    if (any(names(stdyx.estcoeff) %in% z3)) Z7Xz3 <- stdyx.estcoeff[z3]
    if (any(names(stdyx.estcoeff) %in% z4)) Z7Xz4 <- stdyx.estcoeff[z4]

    estX.lo <- (Z7Xa1-Z7Xz1)*(Z7Xa2-Z7Xz2)*(Z7Xa3-Z7Xz3)*(Z7Xa4-Z7Xz4)
    estX.hi <- (Z7Xa1+Z7Xz1)*(Z7Xa2+Z7Xz2)*(Z7Xa3+Z7Xz3)*(Z7Xa4+Z7Xz4)

    # -- Print standardized indirect effects -- #
    cat("\n", "Standardized indirect effects when Z is low (Mean - 1 S.D.) = ", estX.lo)
    cat("\n", "Standardized indirect effects when Z is high (Mean + 1 S.D.) = ", estX.hi, "\n")

    Two_Way_Figure(estX.lo, estX.hi) # Plot standardized Figure with Sub-Function

  } # End one moderating effect for Z


  ## --- Only One Moderating Effect W --- ##
  if ((NoMod == 1) & (NoModw == 1)) {

    # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% a1)) U7Xa1 <- estcoeff[a1]
    if (any(names(estcoeff) %in% a2)) U7Xa2 <- estcoeff[a2]
    if (any(names(estcoeff) %in% a3)) U7Xa3 <- estcoeff[a3]
    if (any(names(estcoeff) %in% a4)) U7Xa4 <- estcoeff[a4]
    if (any(names(estcoeff) %in% w1)) U7Xw1 <- estcoeff[w1]
    if (any(names(estcoeff) %in% w2)) U7Xw2 <- estcoeff[w2]
    if (any(names(estcoeff) %in% w3)) U7Xw3 <- estcoeff[w3]
    if (any(names(estcoeff) %in% w4)) U7Xw4 <- estcoeff[w4]
    stdW <- sqrt(estcoeff[varW])


    # Calculate Estimated Index MM #
    if (PoMod == 1) estIMM <- U7Xw1*U7Xa2*U7Xa3*U7Xa4
    if (PoMod == 2) estIMM <- U7Xa1*U7Xw2*U7Xa3*U7Xa4
    if (PoMod == 3) estIMM <- U7Xa1*U7Xa2*U7Xw3*U7Xa4
    if (PoMod == 4) estIMM <- U7Xa1*U7Xa2*U7Xa3*U7Xw4

    # Capture simulated parameters for calculating indirect effects #
    if (any(names(estcoeff) %in% a1)) U7XSa1 <- mcmc[, a1]
    if (any(names(estcoeff) %in% a2)) U7XSa2 <- mcmc[, a2]
    if (any(names(estcoeff) %in% a3)) U7XSa3 <- mcmc[, a3]
    if (any(names(estcoeff) %in% a4)) U7XSa4 <- mcmc[, a4]
    if (any(names(estcoeff) %in% w1)) U7XSw1 <- mcmc[, w1]
    if (any(names(estcoeff) %in% w2)) U7XSw2 <- mcmc[, w2]
    if (any(names(estcoeff) %in% w3)) U7XSw3 <- mcmc[, w3]
    if (any(names(estcoeff) %in% w4)) U7XSw4 <- mcmc[, w4]
    SstdW <- sqrt(mcmc[,varW])

    # Calculate Simulated Index MM
    if (any(names(estcoeff) %in% w1)) IndexMM <- U7XSw1*U7XSa2*U7XSa3*U7XSa4
    if (any(names(estcoeff) %in% w2)) IndexMM <- U7XSa1*U7XSw2*U7XSa3*U7XSa4
    if (any(names(estcoeff) %in% w3)) IndexMM <- U7XSa1*U7XSa2*U7XSw3*U7XSa4
    if (any(names(estcoeff) %in% w4)) IndexMM <- U7XSa1*U7XSa2*U7XSa3*U7XSw4


    #### Percentile and Bias-Corrected Confidence Intervals of Conditional Indirect Effects ####

    PCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))
    BCCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                               c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    for (X in 1:5) {
      level <- X - 3
      estX <- (U7Xa1+U7Xw1*level*stdW)*(U7Xa2+U7Xw2*level*stdW)*(U7Xa3+U7Xw3*level*stdW)*(U7Xa4+U7Xw4*level*stdW)
      abX <- (U7XSa1+U7XSw1*level*SstdW)*(U7XSa2+U7XSw2*level*SstdW)*(U7XSa3+U7XSw3*level*SstdW)*(U7XSa4+U7XSw4*level*SstdW)
      zX = qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

      ## Percentile Confidence Intervals ##
      PCI[X,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

      # Percentile p-value #
      if (quantile(abX,probs=0.5)>0) {
        PCI[X,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
      } else {
        PCI[X,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
      }

      ## Bias-Corrected Confidence Intervals ##
      BCCI[X,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

      # Bias-Corrected p-value
      if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
        BCCI[X,8] = 0
      } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
        BCCI[X,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
      } else {
        BCCI[X,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
      }
    }

    # Percentile Confidence Intervals for Index MM
    PCI[6,1] <- format(round(quantile(IndexMM,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,2] <- format(round(quantile(IndexMM,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,3] <- format(round(quantile(IndexMM,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,4] <- format(round(estIMM, digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,5] <- format(round(quantile(IndexMM,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,6] <- format(round(quantile(IndexMM,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
    PCI[6,7] <- format(round(quantile(IndexMM,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

    # Percentile Probability for Index MM
    if (quantile(IndexMM,probs=0.5)>0) {
      PCI[6,8] <- format(round(2*(sum(IndexMM<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)  # Probability
    } else {
      PCI[6,8] <- format(round(2*(sum(IndexMM>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)  # Probability
    }


    # Bias-Corrected Confidence Intervals for Index MM
    zIMM = qnorm(sum(IndexMM<estIMM)/b.no)
    BCCI[6,1] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,2] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,3] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,4] <- format(round(estIMM, digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,5] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,6] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCI[6,7] <- format(round(quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

    # Bias-Corrected Probability for Index MM
    if ((estIMM>0 & min(IndexMM)>0) | (estIMM<0 & max(IndexMM)<0)) {
      BCCI[6,8] = 0
    } else if (qnorm(sum(IndexMM>0)/b.no)+2*zIMM<0) {
      BCCI[6,8] = format(round(2*pnorm((qnorm(sum(IndexMM>0)/b.no)+2*zIMM)), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      BCCI[6,8] = format(round(2*pnorm(-1*(qnorm(sum(IndexMM>0)/b.no)+2*zIMM)), digits = 4), nsmall = 4, scientific = FALSE)
    }

    cat("\n")
    cat("Percentile Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(PCI, quote=FALSE, right=TRUE)
    cat("\n")

    cat("\n")
    cat("Bias-Corrected Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(BCCI, quote=FALSE, right=TRUE)
    cat("\n")


    # -- Plot Standardized 2-Way Interaction Effects -- #
    # Define estimated parameters for calculating indirect effects
    if (any(names(stdyx.estcoeff) %in% a1)) Z7Xa1 <- stdyx.estcoeff[a1]
    if (any(names(stdyx.estcoeff) %in% a2)) Z7Xa2 <- stdyx.estcoeff[a2]
    if (any(names(stdyx.estcoeff) %in% a3)) Z7Xa3 <- stdyx.estcoeff[a3]
    if (any(names(stdyx.estcoeff) %in% a4)) Z7Xa4 <- stdyx.estcoeff[a4]
    if (any(names(stdyx.estcoeff) %in% w1)) Z7Xw1 <- stdyx.estcoeff[w1]
    if (any(names(stdyx.estcoeff) %in% w2)) Z7Xw2 <- stdyx.estcoeff[w2]
    if (any(names(stdyx.estcoeff) %in% w3)) Z7Xw3 <- stdyx.estcoeff[w3]
    if (any(names(stdyx.estcoeff) %in% w4)) Z7Xw4 <- stdyx.estcoeff[w4]

    estX.lo <- (Z7Xa1-Z7Xw1)*(Z7Xa2-Z7Xw2)*(Z7Xa3-Z7Xw3)*(Z7Xa4-Z7Xw4)
    estX.hi <- (Z7Xa1+Z7Xw1)*(Z7Xa2+Z7Xw2)*(Z7Xa3+Z7Xw3)*(Z7Xa4+Z7Xw4)

    # -- Print standardized indirect effects -- #
    cat("\n", "Standardized indirect effects when W is low (Mean - 1 S.D.) = ", estX.lo)
    cat("\n", "Standardized indirect effects when W is high (Mean + 1 S.D.) = ", estX.hi, "\n")

    Two_Way_Figure(estX.lo, estX.hi) # Plot standardized Figure with Sub-Function

  } # End one moderating effect for W



  ## --- More Than One Moderating Effect for Z only --- ##
  if ((NoMod > 1) & (NoModw == 0)) {

    # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% a1)) U7Xa1 <- estcoeff[a1]
    if (any(names(estcoeff) %in% a2)) U7Xa2 <- estcoeff[a2]
    if (any(names(estcoeff) %in% a3)) U7Xa3 <- estcoeff[a3]
    if (any(names(estcoeff) %in% a4)) U7Xa4 <- estcoeff[a4]
    if (any(names(estcoeff) %in% z1)) U7Xz1 <- estcoeff[z1]
    if (any(names(estcoeff) %in% z2)) U7Xz2 <- estcoeff[z2]
    if (any(names(estcoeff) %in% z3)) U7Xz3 <- estcoeff[z3]
    if (any(names(estcoeff) %in% z4)) U7Xz4 <- estcoeff[z4]
    stdZ <- sqrt(estcoeff[varZ])

    # Capture simulated parameters for calculating indirect effects #
    if (any(names(estcoeff) %in% a1)) U7XSa1 <- mcmc[, a1]
    if (any(names(estcoeff) %in% a2)) U7XSa2 <- mcmc[, a2]
    if (any(names(estcoeff) %in% a3)) U7XSa3 <- mcmc[, a3]
    if (any(names(estcoeff) %in% a4)) U7XSa4 <- mcmc[, a4]
    if (any(names(estcoeff) %in% z1)) U7XSz1 <- mcmc[, z1]
    if (any(names(estcoeff) %in% z2)) U7XSz2 <- mcmc[, z2]
    if (any(names(estcoeff) %in% z3)) U7XSz3 <- mcmc[, z3]
    if (any(names(estcoeff) %in% z4)) U7XSz4 <- mcmc[, z4]
    SstdZ <- sqrt(mcmc[,varZ])

    #### Percentile and Bias-Corrected Confidence Intervals of Conditional Indirect Effects ####

    PCI <- matrix(1:40, nrow = 5, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))
    BCCI <- matrix(1:40, nrow = 5, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                               c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    for (X in 1:5) {
      level <- X - 3
      estX <- (U7Xa1+U7Xz1*level*stdZ)*(U7Xa2+U7Xz2*level*stdZ)*(U7Xa3+U7Xz3*level*stdZ)*(U7Xa4+U7Xz4*level*stdZ)
      abX <- (U7XSa1+U7XSz1*level*SstdZ)*(U7XSa2+U7XSz2*level*SstdZ)*(U7XSa3+U7XSz3*level*SstdZ)*(U7XSa4+U7XSz4*level*SstdZ)
      zX = qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

      ## Percentile Confidence Intervals ##
      PCI[X,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

      # Percentile p-value #
      if (quantile(abX,probs=0.5)>0) {
        PCI[X,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
      } else {
        PCI[X,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
      }

      ## Bias-Corrected Confidence Intervals ##
      BCCI[X,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

      # Bias-Corrected p-value
      if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
        BCCI[X,8] = 0
      } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
        BCCI[X,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
      } else {
        BCCI[X,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
      }
    }

    cat("\n")
    cat("Percentile Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(PCI, quote=FALSE, right=TRUE)
    cat("\n")
    cat("Index of Moderated Mediation is not defined for two or more moderating effects", rep("\n",3))

    cat("\n")
    cat("Bias-Corrected Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(BCCI, quote=FALSE, right=TRUE)
    cat("\n")
    cat("Index of Moderated Mediation is not defined for two or more moderating effects","\n")

    # -- Plot Standardized 2-Way Interaction Effects -- #
    # Define estimated parameters for calculating indirect effects
    if (any(names(stdyx.estcoeff) %in% a1)) Z7Xa1 <- stdyx.estcoeff[a1]
    if (any(names(stdyx.estcoeff) %in% a2)) Z7Xa2 <- stdyx.estcoeff[a2]
    if (any(names(stdyx.estcoeff) %in% a3)) Z7Xa3 <- stdyx.estcoeff[a3]
    if (any(names(stdyx.estcoeff) %in% a4)) Z7Xa4 <- stdyx.estcoeff[a4]
    if (any(names(stdyx.estcoeff) %in% z1)) Z7Xz1 <- stdyx.estcoeff[z1]
    if (any(names(stdyx.estcoeff) %in% z2)) Z7Xz2 <- stdyx.estcoeff[z2]
    if (any(names(stdyx.estcoeff) %in% z3)) Z7Xz3 <- stdyx.estcoeff[z3]
    if (any(names(stdyx.estcoeff) %in% z4)) Z7Xz4 <- stdyx.estcoeff[z4]

    estX.lo <- (Z7Xa1-Z7Xz1)*(Z7Xa2-Z7Xz2)*(Z7Xa3-Z7Xz3)*(Z7Xa4-Z7Xz4)
    estX.hi <- (Z7Xa1+Z7Xz1)*(Z7Xa2+Z7Xz2)*(Z7Xa3+Z7Xz3)*(Z7Xa4+Z7Xz4)

    # -- Print standardized indirect effects -- #
    cat("\n", "Standardized indirect effects when Z is low (Mean - 1 S.D.) = ", estX.lo)
    cat("\n", "Standardized indirect effects when Z is high (Mean + 1 S.D.) = ", estX.hi, "\n")

    Two_Way_Figure(estX.lo, estX.hi) # Plot standardized Figure with Sub-Function

  } # End more than one moderating effect for Z


  ## --- More Than One Moderating Effect for W only --- ##
  if ((NoMod > 1) & (NoModz == 0)) {

    # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% a1)) U7Xa1 <- estcoeff[a1]
    if (any(names(estcoeff) %in% a2)) U7Xa2 <- estcoeff[a2]
    if (any(names(estcoeff) %in% a3)) U7Xa3 <- estcoeff[a3]
    if (any(names(estcoeff) %in% a4)) U7Xa4 <- estcoeff[a4]
    if (any(names(estcoeff) %in% w1)) U7Xw1 <- estcoeff[w1]
    if (any(names(estcoeff) %in% w2)) U7Xw2 <- estcoeff[w2]
    if (any(names(estcoeff) %in% w3)) U7Xw3 <- estcoeff[w3]
    if (any(names(estcoeff) %in% w4)) U7Xw4 <- estcoeff[w4]
    stdW <- sqrt(estcoeff[varW])

    # Capture simulated parameters for calculating indirect effects #
    if (any(names(estcoeff) %in% a1)) U7XSa1 <- mcmc[, a1]
    if (any(names(estcoeff) %in% a2)) U7XSa2 <- mcmc[, a2]
    if (any(names(estcoeff) %in% a3)) U7XSa3 <- mcmc[, a3]
    if (any(names(estcoeff) %in% a4)) U7XSa4 <- mcmc[, a4]
    if (any(names(estcoeff) %in% w1)) U7XSw1 <- mcmc[, w1]
    if (any(names(estcoeff) %in% w2)) U7XSw2 <- mcmc[, w2]
    if (any(names(estcoeff) %in% w3)) U7XSw3 <- mcmc[, w3]
    if (any(names(estcoeff) %in% w4)) U7XSw4 <- mcmc[, w4]
    SstdW <- sqrt(mcmc[,varW])

    #### Percentile and Bias-Corrected Confidence Intervals of Conditional Indirect Effects ####

    PCI <- matrix(1:40, nrow = 5, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))
    BCCI <- matrix(1:40, nrow = 5, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    for (X in 1:5) {
      level <- X - 3
      estX <- (U7Xa1+U7Xw1*level*stdW)*(U7Xa2+U7Xw2*level*stdW)*(U7Xa3+U7Xw3*level*stdW)*(U7Xa4+U7Xw4*level*stdW)
      abX <- (U7XSa1+U7XSw1*level*SstdW)*(U7XSa2+U7XSw2*level*SstdW)*(U7XSa3+U7XSw3*level*SstdW)*(U7XSa4+U7XSw4*level*SstdW)
      zX = qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

      ## Percentile Confidence Intervals ##
      PCI[X,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
      PCI[X,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

      # Percentile p-value #
      if (quantile(abX,probs=0.5)>0) {
        PCI[X,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
      } else {
        PCI[X,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
      }

      ## Bias-Corrected Confidence Intervals ##
      BCCI[X,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
      BCCI[X,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

      # Bias-Corrected p-value
      if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
        BCCI[X,8] = 0
      } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
        BCCI[X,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
      } else {
        BCCI[X,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
      }
    }


    cat("\n")
    cat("Percentile Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(PCI, quote=FALSE, right=TRUE)
    cat("\n")
    cat("Index of Moderated Mediation is not defined for two or more moderating effects", rep("\n",3))

    cat("\n")
    cat("Bias-Corrected Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(BCCI, quote=FALSE, right=TRUE)
    cat("\n")
    cat("Index of Moderated Mediation is not defined for two or more moderating effects","\n")


    # -- Plot Standardized 2-Way Interaction Effects -- #
    # Define estimated parameters for calculating standardized indirect effects
    if (any(names(stdyx.estcoeff) %in% a1)) Z7Xa1 <- stdyx.estcoeff[a1]
    if (any(names(stdyx.estcoeff) %in% a2)) Z7Xa2 <- stdyx.estcoeff[a2]
    if (any(names(stdyx.estcoeff) %in% a3)) Z7Xa3 <- stdyx.estcoeff[a3]
    if (any(names(stdyx.estcoeff) %in% a4)) Z7Xa4 <- stdyx.estcoeff[a4]
    if (any(names(stdyx.estcoeff) %in% w1)) Z7Xw1 <- stdyx.estcoeff[w1]
    if (any(names(stdyx.estcoeff) %in% w2)) Z7Xw2 <- stdyx.estcoeff[w2]
    if (any(names(stdyx.estcoeff) %in% w3)) Z7Xw3 <- stdyx.estcoeff[w3]
    if (any(names(stdyx.estcoeff) %in% w4)) Z7Xw4 <- stdyx.estcoeff[w4]

    estX.lo <- (Z7Xa1-Z7Xw1)*(Z7Xa2-Z7Xw2)*(Z7Xa3-Z7Xw3)*(Z7Xa4-Z7Xw4)
    estX.hi <- (Z7Xa1+Z7Xw1)*(Z7Xa2+Z7Xw2)*(Z7Xa3+Z7Xw3)*(Z7Xa4+Z7Xw4)

    # -- Print standardized indirect effects -- #
    cat("\n", "Standardized indirect effects when W is low (Mean - 1 S.D.) = ", estX.lo)
    cat("\n", "Standardized indirect effects when W is high (Mean + 1 S.D.) = ", estX.hi, "\n")

    Two_Way_Figure(estX.lo, estX.hi) # Plot standardized Figure with Sub-Function

  } # End more than one moderating effect for W



  ## --- More Than One Moderating Effect for Z & W --- ##
  if ((NoMod > 1) & (NoModz != 0) & (NoModw != 0)) {

  # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% a1)) U7Xa1 <- estcoeff[a1]
    if (any(names(estcoeff) %in% a2)) U7Xa2 <- estcoeff[a2]
    if (any(names(estcoeff) %in% a3)) U7Xa3 <- estcoeff[a3]
    if (any(names(estcoeff) %in% a4)) U7Xa4 <- estcoeff[a4]
    if (any(names(estcoeff) %in% z1)) U7Xz1 <- estcoeff[z1]
    if (any(names(estcoeff) %in% z2)) U7Xz2 <- estcoeff[z2]
    if (any(names(estcoeff) %in% z3)) U7Xz3 <- estcoeff[z3]
    if (any(names(estcoeff) %in% z4)) U7Xz4 <- estcoeff[z4]
    if (any(names(estcoeff) %in% w1)) U7Xw1 <- estcoeff[w1]
    if (any(names(estcoeff) %in% w2)) U7Xw2 <- estcoeff[w2]
    if (any(names(estcoeff) %in% w3)) U7Xw3 <- estcoeff[w3]
    if (any(names(estcoeff) %in% w4)) U7Xw4 <- estcoeff[w4]
    if (any(names(estcoeff) %in% zw1)) U7Xzw1 <- estcoeff[zw1]
    if (any(names(estcoeff) %in% zw2)) U7Xzw2 <- estcoeff[zw2]
    if (any(names(estcoeff) %in% zw3)) U7Xzw3 <- estcoeff[zw3]
    if (any(names(estcoeff) %in% zw4)) U7Xzw4 <- estcoeff[zw4]
    stdZ <- sqrt(estcoeff[varZ])
    stdW <- sqrt(estcoeff[varW])


  # Capture simulated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% a1)) U7XSa1 <- mcmc[, a1]
    if (any(names(estcoeff) %in% a2)) U7XSa2 <- mcmc[, a2]
    if (any(names(estcoeff) %in% a3)) U7XSa3 <- mcmc[, a3]
    if (any(names(estcoeff) %in% a4)) U7XSa4 <- mcmc[, a4]
    if (any(names(estcoeff) %in% z1)) U7XSz1 <- mcmc[, z1]
    if (any(names(estcoeff) %in% z2)) U7XSz2 <- mcmc[, z2]
    if (any(names(estcoeff) %in% z3)) U7XSz3 <- mcmc[, z3]
    if (any(names(estcoeff) %in% z4)) U7XSz4 <- mcmc[, z4]
    if (any(names(estcoeff) %in% w1)) U7XSw1 <- mcmc[, w1]
    if (any(names(estcoeff) %in% w2)) U7XSw2 <- mcmc[, w2]
    if (any(names(estcoeff) %in% w3)) U7XSw3 <- mcmc[, w3]
    if (any(names(estcoeff) %in% w4)) U7XSw4 <- mcmc[, w4]
    if (any(names(estcoeff) %in% zw1)) U7XSzw1 <- mcmc[, zw1]
    if (any(names(estcoeff) %in% zw2)) U7XSzw2 <- mcmc[, zw2]
    if (any(names(estcoeff) %in% zw3)) U7XSzw3 <- mcmc[, zw3]
    if (any(names(estcoeff) %in% zw4)) U7XSzw4 <- mcmc[, zw4]
    SstdZ <- sqrt(mcmc[,varZ])
    SstdW <- sqrt(mcmc[,varW])

    #### Percentile and Bias-Corrected Confidence Intervals of Conditional Indirect Effects ####

    PCI <- array(data = 1:200,
             dim = c(5, 8, 5),
             dimnames = list(c("Z = Mean-2SD","Z = Mean-1SD", "Z = Mean", "Z = Mean+1SD", "Z = Mean+2SD"),
                             c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value"),
                             c("W = Mean-2SD","W = Mean - 1SD", "W = Mean", "W = Mean+1SD", "W = Mean+2SD")))
    BCCI <- array(data = 1:200,
             dim = c(5, 8, 5),
             dimnames = list(c("Z = Mean-2SD","Z = Mean-1SD", "Z = Mean", "Z = Mean+1SD", "Z = Mean+2SD"),
                             c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value"),
                             c("W = Mean-2SD","W = Mean - 1SD", "W = Mean", "W = Mean+1SD", "W = Mean+2SD")))

    for (Z in 1:5) {
      for (W in 1:5) {
        levelZ <- Z - 3
        levelW <- W - 3
        estX <- (U7Xa1+U7Xz1*levelZ*stdZ+U7Xw1*levelW*stdW+U7Xzw1*levelZ*stdZ*levelW*stdW)*
                (U7Xa2+U7Xz2*levelZ*stdZ+U7Xw2*levelW*stdW+U7Xzw2*levelZ*stdZ*levelW*stdW)*
                (U7Xa3+U7Xz3*levelZ*stdZ+U7Xw3*levelW*stdW+U7Xzw3*levelZ*stdZ*levelW*stdW)*
                (U7Xa4+U7Xz4*levelZ*stdZ+U7Xw4*levelW*stdW+U7Xzw4*levelZ*stdZ*levelW*stdW)
        abX <- (U7XSa1+U7XSz1*levelZ*SstdZ+U7XSw1*levelW*SstdW+U7XSzw1*levelZ*SstdZ*levelW*SstdW)*
               (U7XSa2+U7XSz2*levelZ*SstdZ+U7XSw2*levelW*SstdW+U7XSzw2*levelZ*SstdZ*levelW*SstdW)*
               (U7XSa3+U7XSz3*levelZ*SstdZ+U7XSw3*levelW*SstdW+U7XSzw3*levelZ*SstdZ*levelW*SstdW)*
               (U7XSa4+U7XSz4*levelZ*SstdZ+U7XSw4*levelW*SstdW+U7XSzw4*levelZ*SstdZ*levelW*SstdW)
        zX = qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

        ## Percentile Confidence Intervals ##
        PCI[Z,1,W] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
        PCI[Z,2,W] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
        PCI[Z,3,W] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
        PCI[Z,4,W] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
        PCI[Z,5,W] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
        PCI[Z,6,W] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
        PCI[Z,7,W] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

        # Percentile p-value #
        if (quantile(abX,probs=0.5)>0) {
          PCI[Z,8,W] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
        } else {
          PCI[Z,8,W] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
        }

        ## Bias-Corrected Confidence Intervals ##
        BCCI[Z,1,W] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
        BCCI[Z,2,W] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
        BCCI[Z,3,W] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
        BCCI[Z,4,W] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
        BCCI[Z,5,W] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
        BCCI[Z,6,W] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
        BCCI[Z,7,W] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)
        # Bias-Corrected Probability
        if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
          BCCI[Z,8,W] = 0
        } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
          BCCI[Z,8,W] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
        } else {
          BCCI[Z,8,W] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
        }
      }
    }


    cat("\n")
    cat("Percentile Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(PCI, quote=FALSE, right=TRUE)
    cat("\n")
    cat("Note: Index of Moderated Mediation is not defined for two or more moderating effects", rep("\n",3))

    cat("\n")
    cat("Bias-Corrected Confidence Intervals for Conditional Indirect Effects", rep("\n", 2))
    print(BCCI, quote=FALSE, right=TRUE)
    cat("\n")
    cat("Note: Index of Moderated Mediation is not defined for two or more moderating effects","\n")

    ## Slope Difference Tests ##
    PCISDT <- matrix(1:48, nrow = 6, dimnames = list(c("HiZ/HiW - HiZ/LoW", "HiZ/HiW - LoZ/HiW", "HiZ/LoW - LoZ/LoW",
                                                       "LoZ/HiW - LoZ/LoW", "HiZ/HiW - LoZ/LoW", "HiZ/LoW - LoZ/HiW"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))
    BCCISDT <- matrix(1:48, nrow = 6, dimnames = list(c("HiZ/HiW - HiZ/LoW", "HiZ/HiW - LoZ/HiW", "HiZ/LoW - LoZ/LoW",
                                                      "LoZ/HiW - LoZ/LoW", "HiZ/HiW - LoZ/LoW", "HiZ/LoW - LoZ/HiW"),
                                               c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))


    ## Slope Difference Test 1 "HiZ/HiW - HiZ/LoW" ##
    estX <- (U7Xa1+U7Xz1*stdZ+U7Xw1*stdW+U7Xzw1*stdZ*stdW)*(U7Xa2+U7Xz2*stdZ+U7Xw2*stdW+U7Xzw2*stdZ*stdW)*
            (U7Xa3+U7Xz3*stdZ+U7Xw3*stdW+U7Xzw3*stdZ*stdW)*(U7Xa4+U7Xz4*stdZ+U7Xw4*stdW+U7Xzw4*stdZ*stdW) -
            (U7Xa1+U7Xz1*stdZ-U7Xw1*stdW-U7Xzw1*stdZ*stdW)*(U7Xa2+U7Xz2*stdZ-U7Xw2*stdW-U7Xzw2*stdZ*stdW)*
            (U7Xa3+U7Xz3*stdZ-U7Xw3*stdW-U7Xzw3*stdZ*stdW)*(U7Xa4+U7Xz4*stdZ-U7Xw4*stdW-U7Xzw4*stdZ*stdW)
    abX <- (U7XSa1+U7XSz1*SstdZ+U7XSw1*SstdW+U7XSzw1*SstdZ*SstdW)*(U7XSa2+U7XSz2*SstdZ+U7XSw2*SstdW+U7XSzw2*SstdZ*SstdW)*
           (U7XSa3+U7XSz3*SstdZ+U7XSw3*SstdW+U7XSzw3*SstdZ*SstdW)*(U7XSa4+U7XSz4*SstdZ+U7XSw4*SstdW+U7XSzw4*SstdZ*SstdW) -
           (U7XSa1+U7XSz1*SstdZ-U7XSw1*SstdW-U7XSzw1*SstdZ*SstdW)*(U7XSa2+U7XSz2*SstdZ-U7XSw2*SstdW-U7XSzw2*SstdZ*SstdW)*
           (U7XSa3+U7XSz3*SstdZ-U7XSw3*SstdW-U7XSzw3*SstdZ*SstdW)*(U7XSa4+U7XSz4*SstdZ-U7XSw4*SstdW-U7XSzw4*SstdZ*SstdW)
    zX = qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

    ## Percentile Confidence Intervals ##
    PCISDT[1,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[1,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[1,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[1,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[1,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[1,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[1,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

    # Percentile p-value #
    if (quantile(abX,probs=0.5)>0) {
      PCISDT[1,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      PCISDT[1,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    }

    ## Bias-Corrected Confidence Intervals ##
    BCCISDT[1,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[1,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[1,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[1,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[1,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[1,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[1,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

    # Bias-Corrected Probability
    if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
      BCCISDT[1,8] = 0
    } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
      BCCISDT[1,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      BCCISDT[1,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    }

    ## Slope Difference Test 2 "HiZ/HiW - LoZ/HiW" ##
    estX <- (U7Xa1+U7Xz1*stdZ+U7Xw1*stdW+U7Xzw1*stdZ*stdW)*(U7Xa2+U7Xz2*stdZ+U7Xw2*stdW+U7Xzw2*stdZ*stdW)*
            (U7Xa3+U7Xz3*stdZ+U7Xw3*stdW+U7Xzw3*stdZ*stdW)*(U7Xa4+U7Xz4*stdZ+U7Xw4*stdW+U7Xzw4*stdZ*stdW) -
            (U7Xa1-U7Xz1*stdZ+U7Xw1*stdW-U7Xzw1*stdZ*stdW)*(U7Xa2-U7Xz2*stdZ+U7Xw2*stdW-U7Xzw2*stdZ*stdW)*
            (U7Xa3-U7Xz3*stdZ+U7Xw3*stdW-U7Xzw3*stdZ*stdW)*(U7Xa4-U7Xz4*stdZ+U7Xw4*stdW-U7Xzw4*stdZ*stdW)
    abX <- (U7XSa1+U7XSz1*SstdZ+U7XSw1*SstdW+U7XSzw1*SstdZ*SstdW)*(U7XSa2+U7XSz2*SstdZ+U7XSw2*SstdW+U7XSzw2*SstdZ*SstdW)*
           (U7XSa3+U7XSz3*SstdZ+U7XSw3*SstdW+U7XSzw3*SstdZ*SstdW)*(U7XSa4+U7XSz4*SstdZ+U7XSw4*SstdW+U7XSzw4*SstdZ*SstdW) -
           (U7XSa1-U7XSz1*SstdZ+U7XSw1*SstdW-U7XSzw1*SstdZ*SstdW)*(U7XSa2-U7XSz2*SstdZ+U7XSw2*SstdW-U7XSzw2*SstdZ*SstdW)*
           (U7XSa3-U7XSz3*SstdZ+U7XSw3*SstdW-U7XSzw3*SstdZ*SstdW)*(U7XSa4-U7XSz4*SstdZ+U7XSw4*SstdW-U7XSzw4*SstdZ*SstdW)
    zX = qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

    ## Percentile Confidence Intervals ##
    PCISDT[2,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[2,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[2,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[2,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[2,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[2,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[2,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

    # Percentile p-value #
    if (quantile(abX,probs=0.5)>0) {
      PCISDT[2,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      PCISDT[2,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    }

    ## Bias-Corrected Confidence Intervals ##
    BCCISDT[2,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[2,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[2,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[2,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[2,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[2,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[2,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

    # Bias-Corrected Probability
    if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
      BCCISDT[2,8] = 0
    } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
      BCCISDT[2,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      BCCISDT[2,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    }

    ## Slope Difference Test 3 "HiZ/LoW - LoZ/LoW" ##
    estX <- (U7Xa1+U7Xz1*stdZ-U7Xw1*stdW-U7Xzw1*stdZ*stdW)*(U7Xa2+U7Xz2*stdZ-U7Xw2*stdW-U7Xzw2*stdZ*stdW)*
            (U7Xa3+U7Xz3*stdZ-U7Xw3*stdW-U7Xzw3*stdZ*stdW)*(U7Xa4+U7Xz4*stdZ-U7Xw4*stdW-U7Xzw4*stdZ*stdW) -
            (U7Xa1-U7Xz1*stdZ-U7Xw1*stdW+U7Xzw1*stdZ*stdW)*(U7Xa2-U7Xz2*stdZ-U7Xw2*stdW+U7Xzw2*stdZ*stdW)*
            (U7Xa3-U7Xz3*stdZ-U7Xw3*stdW+U7Xzw3*stdZ*stdW)*(U7Xa4-U7Xz4*stdZ-U7Xw4*stdW+U7Xzw4*stdZ*stdW)
    abX <- (U7XSa1+U7XSz1*SstdZ-U7XSw1*SstdW-U7XSzw1*SstdZ*SstdW)*(U7XSa2+U7XSz2*SstdZ-U7XSw2*SstdW-U7XSzw2*SstdZ*SstdW)*
           (U7XSa3+U7XSz3*SstdZ-U7XSw3*SstdW-U7XSzw3*SstdZ*SstdW)*(U7XSa4+U7XSz4*SstdZ-U7XSw4*SstdW-U7XSzw4*SstdZ*SstdW) -
           (U7XSa1-U7XSz1*SstdZ-U7XSw1*SstdW+U7XSzw1*SstdZ*SstdW)*(U7XSa2-U7XSz2*SstdZ-U7XSw2*SstdW+U7XSzw2*SstdZ*SstdW)*
           (U7XSa3-U7XSz3*SstdZ-U7XSw3*SstdW+U7XSzw3*SstdZ*SstdW)*(U7XSa4-U7XSz4*SstdZ-U7XSw4*SstdW+U7XSzw4*SstdZ*SstdW)
    zX = qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

    ## Percentile Confidence Intervals ##
    PCISDT[3,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[3,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[3,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[3,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[3,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[3,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[3,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

    # Percentile p-value #
    if (quantile(abX,probs=0.5)>0) {
      PCISDT[3,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      PCISDT[3,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    }

    ## Bias-Corrected Confidence Intervals ##
    BCCISDT[3,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[3,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[3,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[3,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[3,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[3,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[3,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

    # Bias-Corrected Probability
    if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
      BCCISDT[3,8] = 0
    } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
      BCCISDT[3,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      BCCISDT[3,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    }

    ## Slope Difference Test 4 "LoZ/HiW - LoZ/LoW" ##
    estX <- (U7Xa1-U7Xz1*stdZ+U7Xw1*stdW-U7Xzw1*stdZ*stdW)*(U7Xa2-U7Xz2*stdZ+U7Xw2*stdW-U7Xzw2*stdZ*stdW)*
            (U7Xa3-U7Xz3*stdZ+U7Xw3*stdW-U7Xzw3*stdZ*stdW)*(U7Xa4-U7Xz4*stdZ+U7Xw4*stdW-U7Xzw4*stdZ*stdW) -
            (U7Xa1-U7Xz1*stdZ-U7Xw1*stdW+U7Xzw1*stdZ*stdW)*(U7Xa2-U7Xz2*stdZ-U7Xw2*stdW+U7Xzw2*stdZ*stdW)*
            (U7Xa3-U7Xz3*stdZ-U7Xw3*stdW+U7Xzw3*stdZ*stdW)*(U7Xa4-U7Xz4*stdZ-U7Xw4*stdW+U7Xzw4*stdZ*stdW)
    abX <- (U7XSa1-U7XSz1*SstdZ+U7XSw1*SstdW-U7XSzw1*SstdZ*SstdW)*(U7XSa2-U7XSz2*SstdZ+U7XSw2*SstdW-U7XSzw2*SstdZ*SstdW)*
           (U7XSa3-U7XSz3*SstdZ+U7XSw3*SstdW-U7XSzw3*SstdZ*SstdW)*(U7XSa4-U7XSz4*SstdZ+U7XSw4*SstdW-U7XSzw4*SstdZ*SstdW) -
           (U7XSa1-U7XSz1*SstdZ-U7XSw1*SstdW+U7XSzw1*SstdZ*SstdW)*(U7XSa2-U7XSz2*SstdZ-U7XSw2*SstdW+U7XSzw2*SstdZ*SstdW)*
           (U7XSa3-U7XSz3*SstdZ-U7XSw3*SstdW+U7XSzw3*SstdZ*SstdW)*(U7XSa4-U7XSz4*SstdZ-U7XSw4*SstdW+U7XSzw4*SstdZ*SstdW)
    zX <- qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

    ## Percentile Confidence Intervals ##
    PCISDT[4,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[4,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[4,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[4,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[4,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[4,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[4,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

    # Percentile p-value #
    if (quantile(abX,probs=0.5)>0) {
      PCISDT[4,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      PCISDT[4,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    }

    ## Bias-Corrected Confidence Intervals ##
    BCCISDT[4,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[4,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[4,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[4,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[4,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[4,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[4,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

    # Bias-Corrected Probability
    if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
      BCCISDT[4,8] = 0
    } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
    BCCISDT[4,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      BCCISDT[4,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    }


    ## Slope Difference Test 5 "HiZ/HiW - LoZ/LoW" ##
    estX <- (U7Xa1+U7Xz1*stdZ+U7Xw1*stdW+U7Xzw1*stdZ*stdW)*(U7Xa2+U7Xz2*stdZ+U7Xw2*stdW+U7Xzw2*stdZ*stdW)*
            (U7Xa3+U7Xz3*stdZ+U7Xw3*stdW+U7Xzw3*stdZ*stdW)*(U7Xa4+U7Xz4*stdZ+U7Xw4*stdW+U7Xzw4*stdZ*stdW) -
            (U7Xa1-U7Xz1*stdZ-U7Xw1*stdW+U7Xzw1*stdZ*stdW)*(U7Xa2-U7Xz2*stdZ-U7Xw2*stdW+U7Xzw2*stdZ*stdW)*
            (U7Xa3-U7Xz3*stdZ-U7Xw3*stdW+U7Xzw3*stdZ*stdW)*(U7Xa4-U7Xz4*stdZ-U7Xw4*stdW+U7Xzw4*stdZ*stdW)
    abX <- (U7XSa1+U7XSz1*SstdZ+U7XSw1*SstdW+U7XSzw1*SstdZ*SstdW)*(U7XSa2+U7XSz2*SstdZ+U7XSw2*SstdW+U7XSzw2*SstdZ*SstdW)*
           (U7XSa3+U7XSz3*SstdZ+U7XSw3*SstdW+U7XSzw3*SstdZ*SstdW)*(U7XSa4+U7XSz4*SstdZ+U7XSw4*SstdW+U7XSzw4*SstdZ*SstdW) -
           (U7XSa1-U7XSz1*SstdZ-U7XSw1*SstdW+U7XSzw1*SstdZ*SstdW)*(U7XSa2-U7XSz2*SstdZ-U7XSw2*SstdW+U7XSzw2*SstdZ*SstdW)*
           (U7XSa3-U7XSz3*SstdZ-U7XSw3*SstdW+U7XSzw3*SstdZ*SstdW)*(U7XSa4-U7XSz4*SstdZ-U7XSw4*SstdW+U7XSzw4*SstdZ*SstdW)
    zX <- qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

    ## Percentile Confidence Intervals ##
    PCISDT[5,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[5,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[5,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[5,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[5,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[5,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[5,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

    # Percentile p-value #
    if (quantile(abX,probs=0.5)>0) {
      PCISDT[5,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      PCISDT[5,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    }

    ## Bias-Corrected Confidence Intervals ##
    BCCISDT[5,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[5,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[5,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[5,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[5,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[5,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[5,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

    # Bias-Corrected Probability
    if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
      BCCISDT[5,8] = 0
    } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
      BCCISDT[5,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      BCCISDT[5,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    }


    ## Slope Difference Test 6 "HiZ/LoW - LoZ/HiW" ##
    estX <- (U7Xa1+U7Xz1*stdZ-U7Xw1*stdW-U7Xzw1*stdZ*stdW)*(U7Xa2+U7Xz2*stdZ-U7Xw2*stdW-U7Xzw2*stdZ*stdW)*
            (U7Xa3+U7Xz3*stdZ-U7Xw3*stdW-U7Xzw3*stdZ*stdW)*(U7Xa4+U7Xz4*stdZ-U7Xw4*stdW-U7Xzw4*stdZ*stdW) -
            (U7Xa1-U7Xz1*stdZ+U7Xw1*stdW-U7Xzw1*stdZ*stdW)*(U7Xa2-U7Xz2*stdZ+U7Xw2*stdW-U7Xzw2*stdZ*stdW)*
            (U7Xa3-U7Xz3*stdZ+U7Xw3*stdW-U7Xzw3*stdZ*stdW)*(U7Xa4-U7Xz4*stdZ+U7Xw4*stdW-U7Xzw4*stdZ*stdW)
     abX <- (U7XSa1+U7XSz1*SstdZ-U7XSw1*SstdW-U7XSzw1*SstdZ*SstdW)*(U7XSa2+U7XSz2*SstdZ-U7XSw2*SstdW-U7XSzw2*SstdZ*SstdW)*
            (U7XSa3+U7XSz3*SstdZ-U7XSw3*SstdW-U7XSzw3*SstdZ*SstdW)*(U7XSa4+U7XSz4*SstdZ-U7XSw4*SstdW-U7XSzw4*SstdZ*SstdW) -
            (U7XSa1-U7XSz1*SstdZ+U7XSw1*SstdW-U7XSzw1*SstdZ*SstdW)*(U7XSa2-U7XSz2*SstdZ+U7XSw2*SstdW-U7XSzw2*SstdZ*SstdW)*
            (U7XSa3-U7XSz3*SstdZ+U7XSw3*SstdW-U7XSzw3*SstdZ*SstdW)*(U7XSa4-U7XSz4*SstdZ+U7XSw4*SstdW-U7XSzw4*SstdZ*SstdW)
      zX <- qnorm(sum(abX<estX)/b.no)  # Bias-Corrected Factor

    ## Percentile Confidence Intervals ##
    PCISDT[6,1] <- format(round(quantile(abX,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[6,2] <- format(round(quantile(abX,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[6,3] <- format(round(quantile(abX,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[6,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[6,5] <- format(round(quantile(abX,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[6,6] <- format(round(quantile(abX,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
    PCISDT[6,7] <- format(round(quantile(abX,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)

    # Percentile p-value #
    if (quantile(abX,probs=0.5)>0) {
      PCISDT[6,8] <- format(round(2*(sum(abX<0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      PCISDT[6,8] <- format(round(2*(sum(abX>0)/b.no), digits = 4), nsmall = 4, scientific = FALSE)
    }

    ## Bias-Corrected Confidence Intervals ##
    BCCISDT[6,1] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[6,2] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[6,3] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[6,4] <- format(round(estX, digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[6,5] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[6,6] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
    BCCISDT[6,7] <- format(round(quantile(abX,probs=pnorm(2*zX+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)

    # Bias-Corrected Probability
    if ((estX>0 & min(abX)>0) | (estX<0 & max(abX)<0)) {
      BCCISDT[6,8] = 0
    } else if (qnorm(sum(abX>0)/b.no)+2*zX<0) {
      BCCISDT[6,8] = format(round(2*pnorm((qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    } else {
      BCCISDT[6,8] = format(round(2*pnorm(-1*(qnorm(sum(abX>0)/b.no)+2*zX)), digits = 4), nsmall = 4, scientific = FALSE)
    }

    cat("\n")
    cat("\n")
    cat("Percentile Confidence Intervals for Simple Slope Tests: Hi = +1SD; Lo = -1SD", rep("\n", 2))
    print(PCISDT, quote=FALSE, right=TRUE)
    cat("\n")

    cat("\n")
    cat("Bias-Corrected Confidence Intervals for Simple Slope Tests: Hi = +1SD; Lo = -1SD", rep("\n", 2))
    print(BCCISDT, quote=FALSE, right=TRUE)
    cat("\n")


    ## == Create 3-way Standardized Interaction Figure == ##

    # Define estimated parameters for calculating indirect effects
    if (any(names(stdyx.estcoeff) %in% a1)) Z7Xa1 <- stdyx.estcoeff[a1]
    if (any(names(stdyx.estcoeff) %in% a2)) Z7Xa2 <- stdyx.estcoeff[a2]
    if (any(names(stdyx.estcoeff) %in% a3)) Z7Xa3 <- stdyx.estcoeff[a3]
    if (any(names(stdyx.estcoeff) %in% a4)) Z7Xa4 <- stdyx.estcoeff[a4]
    if (any(names(stdyx.estcoeff) %in% z1)) Z7Xz1 <- stdyx.estcoeff[z1]
    if (any(names(stdyx.estcoeff) %in% z2)) Z7Xz2 <- stdyx.estcoeff[z2]
    if (any(names(stdyx.estcoeff) %in% z3)) Z7Xz3 <- stdyx.estcoeff[z3]
    if (any(names(stdyx.estcoeff) %in% z4)) Z7Xz4 <- stdyx.estcoeff[z4]
    if (any(names(stdyx.estcoeff) %in% w1)) Z7Xw1 <- stdyx.estcoeff[w1]
    if (any(names(stdyx.estcoeff) %in% w2)) Z7Xw2 <- stdyx.estcoeff[w2]
    if (any(names(stdyx.estcoeff) %in% w3)) Z7Xw3 <- stdyx.estcoeff[w3]
    if (any(names(stdyx.estcoeff) %in% w4)) Z7Xw4 <- stdyx.estcoeff[w4]
    if (any(names(stdyx.estcoeff) %in% zw1)) Z7Xzw1 <- stdyx.estcoeff[zw1]
    if (any(names(stdyx.estcoeff) %in% zw2)) Z7Xzw2 <- stdyx.estcoeff[zw2]
    if (any(names(stdyx.estcoeff) %in% zw3)) Z7Xzw3 <- stdyx.estcoeff[zw3]
    if (any(names(stdyx.estcoeff) %in% zw4)) Z7Xzw4 <- stdyx.estcoeff[zw4]

    ## Unstandardized Slopes ##
#    TWx1 <- PCI[4,4,4] # Hi-Z, Hi-W
#    TWx2 <- PCI[4,4,2] # Hi-Z, Lo-W
#    TWx3 <- PCI[2,4,4] # Lo-Z, Hi-W
#    TWx4 <- PCI[2,4,2] # Lo-Z, Lo-W

    ## Standardized Slopes ##
    TWx1 <- (Z7Xa1+Z7Xz1+Z7Xw1+Z7Xzw1)*(Z7Xa2+Z7Xz2+Z7Xw2+Z7Xzw2)*(Z7Xa3+Z7Xz3+Z7Xw3+Z7Xzw3)*(Z7Xa4+Z7Xz4+Z7Xw4+Z7Xzw4) # Hi-Z, Hi-W
    TWx2 <- (Z7Xa1+Z7Xz1-Z7Xw1-Z7Xzw1)*(Z7Xa2+Z7Xz2-Z7Xw2-Z7Xzw2)*(Z7Xa3+Z7Xz3-Z7Xw3-Z7Xzw3)*(Z7Xa4+Z7Xz4-Z7Xw4-Z7Xzw4) # Hi-Z, Lo-W
    TWx3 <- (Z7Xa1-Z7Xz1+Z7Xw1-Z7Xzw1)*(Z7Xa2-Z7Xz2+Z7Xw2-Z7Xzw2)*(Z7Xa3-Z7Xz3+Z7Xw3-Z7Xzw3)*(Z7Xa4-Z7Xz4+Z7Xw4-Z7Xzw4) # Lo-Z, Hi-W
    TWx4 <- (Z7Xa1-Z7Xz1-Z7Xw1+Z7Xzw1)*(Z7Xa2-Z7Xz2-Z7Xw2+Z7Xzw2)*(Z7Xa3-Z7Xz3-Z7Xw3+Z7Xzw3)*(Z7Xa4-Z7Xz4-Z7Xw4+Z7Xzw4) # Lo-Z, Lo-W

    # -- Print standardized indirect effects -- #
    cat("\n", "Standardized indirect effects when Z is high (Mean + 1 S.D.) and W is high (Mean + 1 S.D.) = ", TWx1)
    cat("\n", "Standardized indirect effects when Z is high (Mean + 1 S.D.) and W is low (Mean - 1 S.D.) = ", TWx2)
    cat("\n", "Standardized indirect effects when Z is low (Mean - 1 S.D.) and W is high (Mean + 1 S.D.) = ", TWx3)
    cat("\n", "Standardized indirect effects when Z is low (Mean - 1 S.D.) and W is low (Mean - 1 S.D.) = ", TWx4, "\n")


    df_wide <- data.frame(
      x_var = c(-2,2),
      line1 = c(TWx1*-2, TWx1*2),
      line2 = c(TWx2*-2, TWx2*2),
      line3 = c(TWx3*-2, TWx3*2),
      line4 = c(TWx4*-2, TWx4*2)
    )

    # Convert to long format
    df_long <- df_wide %>%
      pivot_longer(
        cols = starts_with("line"),
        names_to = "line_id",
        values_to = "y_value"
      )
    y.upper = 1.5*max(df_long$y_value) - 0.5*mean(df_long$y_value)
    y.lower = 1.5*min(df_long$y_value) - 0.5*mean(df_long$y_value)

    p_int <<- ggplot(data = df_long, aes(x = x_var, y = y_value, color = line_id, linetype = line_id)) +
                    xlim(-2.5, 2.5) +
                    ylim(y.lower, y.upper) +
                    geom_line(linewidth=1) +
                    scale_linetype_manual(name = "Moderator Levels",
                                          values = c("line1" = "solid", "line2" = "twodash", "line3" = "solid", "line4" = "twodash"),
                                          labels = c("Hi-Z, Hi-W", "Hi-Z, Lo-W", "Lo-Z, Hi-W", "Lo-Z, Lo-W")) +
                    scale_color_manual(values = c("line1" = "black", "line2" = "black", "line3" = "grey", "line4"="grey")) +
                    guides(color = "none") +
                    labs(title = "Standardized Three-Way Interaction Effects",
                             x = "X-axis Label",
                             y = "Indirect Effect of X on Y (through M)") +
                    theme_classic() # Optional: Apply a theme

    ggplot2::ggsave("3-Way Interaction Figure.png", width = 22.86, height = 16.51, units = "cm")

    cat("\n")
    cat("Figure p_int is saved in '3-Way Interaction Figure.png'", rep("\n", 2))

    ## == Close 3-Way Interaction Figure

  }


  ## -- Prepare Object for Outputs -- ##
  if (NoModz == 1) {
    sd_z <- stdZ
    Z <- Z
  } else if (NoModw == 1) {
    sd_z <- stdW
    Z <- W
  } else {
    sd_z <- 0
    Z <- "NA"
  }

  return(list(a1=U7Xa1, a2=U7Xa2, a3=U7Xa3, a4=U7Xa4,
              z1=U7Xz1, z2=U7Xz2, z3=U7Xz3, z4=U7Xz4,
              w1=U7Xw1, w2=U7Xw2, w3=U7Xw3, w4=U7Xw4, sd_z=sd_z,
              Sa1=U7XSa1, Sa2=U7XSa2, Sa3=U7XSa3, Sa4=U7XSa4,
              Sz1=U7XSz1, Sz2=U7XSz2, Sz3=U7XSz3, Sz4=U7XSz4,
              Sw1=U7XSw1, Sw2=U7XSw2, Sw3=U7XSw3, Sw4=U7XSw4,
              b.no = b.no, NoModz = NoModz, NoModw = NoModw, Z=Z))

}  ## ===== End (function mccimm) ===== ##



## ===== Simulation of Defined Function (modsem) ===== ##
#' Monte Carlo Simulation for Confidence Intervals of Defined Function (modsem)
#'
#' Generate confidence intervals of defined function from modsem results using Monte Carlo simulation.
#'
#' @param object modsem object (output from modsem).
#' @param Sfunction function of estimated parameters from modsem object.
#' @param R number of Monte Carlo simulation samples (in millions). For example, R=5 (default) generates 5,000,000 simulated samples.
#'
#' @return mccimm confidence intervals.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' # modsem object is "est_lms"
#'
#' mccimm_modsem_fun(est_lms, Sfunction = "a1*a2")
#'

mccimm_modsem_fun <- function(object, Sfunction="NULL", R=5) {

  dp <- all.vars(parse(text = Sfunction))
  var_label <- parameter_estimates(object)$label
  var_label <- var_label[var_label != ""]
  if (all(unlist(dp) %in% unlist(var_label)) == FALSE) {
    stop("Element(s) in the Sfunction not an estimated parameter in modsem object")
  }

  ## Extract defined parameters and vcov ##


  temp <- modsem_coef(object)
  estcoeff <- temp[dp]
  Temp3 <- modsem_vcov(object)
  Tech3 <- Temp3[dp, dp]


  ## -- Monte Carlo Simulation of R*1e6 samples, default: R = 5 -- ##
  mcmc <- MASS::mvrnorm(n=R*1e6, mu=estcoeff, Sigma=Tech3, tol = 1e-6)

  b.no <- nrow(mcmc)
  R.no <- format(R*1e6, scientific = FALSE)

  # ===== Print number of simulated samples
  cat("\n", "   Number of requested simulated samples = ", R.no)
  cat("\n", "   Number of completed simulated samples = ", b.no, rep("\n",2))


  cat("Simulated Defined Function Values", rep("\n", 2))

  # ==================================================================== #

  # Calculate estimated parameter from Dfunction
  list2env(as.list(estcoeff), envir = .GlobalEnv)
  estM  <- eval(parse(text=Sfunction))


  # Calculate Simulated parameter from Dfunction
  mcmc <- as.data.frame(mcmc)
  mcmc <- mcmc %>%
    mutate(abM = eval(parse(text=Sfunction)))
  abM <- mcmc[, "abM"]

  #### Confidence Intervals and p-value ####

  # Calculate Percentile Probability
  if (quantile(abM,probs=0.5)>0) {
    pM = 2*(sum(abM<0)/b.no)
  } else {
    pM = 2*(sum(abM>0)/b.no)
  }

  #### Percentile Confidence Intervals of Conditional Indirect Effects ####
  PCI <- matrix(1:8, nrow = 1, dimnames = list(c("        "),
                                             c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

  PCI[1,1] <- format(round(quantile(abM,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,2] <- format(round(quantile(abM,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,3] <- format(round(quantile(abM,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,4] <- format(round(estM, digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,5] <- format(round(quantile(abM,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,6] <- format(round(quantile(abM,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,7] <- format(round(quantile(abM,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,8] <- format(round(pM, digits = 4), nsmall = 4, scientific = FALSE)

  # Bias-Corrected Factor
  zM = qnorm(sum(abM<estM)/b.no)

  # Calculate Bias-Corrected Probability

  if ((estM>0 & min(abM)>0) | (estM<0 & max(abM)<0)) {
    pbM = 0
  } else if (qnorm(sum(abM>0)/b.no)+2*zM<0) {
    pbM = 2*pnorm((qnorm(sum(abM>0)/b.no)+2*zM))
  } else {
    pbM = 2*pnorm(-1*(qnorm(sum(abM>0)/b.no)+2*zM))
  }

  #### Bias-Corrected Confidence Intervals ####

  BCCI <- matrix(1:8, nrow = 1, dimnames = list(c("        "),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

  BCCI[1,1] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,2] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,3] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,4] <- format(round(estM, digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,5] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,6] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,7] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,8] <- format(round(pbM, digits = 4), nsmall = 4, scientific = FALSE)

  cat("\n")
  cat("Percentile Confidence Intervals for Sfunction", rep("\n", 2))
  rownames(PCI) <- rep("    ", nrow(PCI))
  print(PCI, quote=FALSE, right=TRUE)
  cat("\n")

  cat("\n")
  cat("Bias-Corrected Confidence Intervals for Sfunction", rep("\n", 2))
  rownames(BCCI) <- rep("    ", nrow(BCCI))
  print(BCCI, quote=FALSE, right=TRUE)
  cat("\n")

}  ## ===== End (function mccimm_modsem_fun) ===== ##





## ===== Simulation of Defined Function (lavaan) ===== ##
#' Monte Carlo Simulation for Confidence Intervals of Defined Function (lavaan)
#'
#' Generate confidence intervals of defined function from lavaan results using Monte Carlo simulation.
#'
#' @param object lavaan object (output from lavaan).
#' @param Sfunction function of estimated parameters from lavaan object.
#' @param R number of Monte Carlo simulation samples (in millions). For example, R=5 (default) generates 5,000,000 simulated samples.
#'
#' @return mccimm confidence intervals.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' # lavaan object is "est_lms"
#'
#' mccimm_lavaan_fun(est_lms, Sfunction = "a1*a2")
#'

mccimm_lavaan_fun <- function(object, Sfunction="NULL", R=5) {

  dp <- all.vars(parse(text = Sfunction))
  var_label <- lavaan::parameterEstimates(object)$label
  var_label <- var_label[var_label != ""]
  if (all(unlist(dp) %in% unlist(var_label)) == FALSE) {
    stop("Element(s) in the Sfunction not an estimated parameter in modsem object")
  }

  ## Extract defined parameters and vcov ##
  temp <- lavaan::coef(object)
  estcoeff <- temp[dp]
  Temp3 <- lavaan::vcov(object)
  Tech3 <- Temp3[dp, dp]

  ## -- Monte Carlo Simulation of R*1e6 samples, default: R = 5 -- ##
  mcmc <- MASS::mvrnorm(n=R*1e6, mu=estcoeff, Sigma=Tech3, tol = 1e-6)

  b.no <- nrow(mcmc)
  R.no <- format(R*1e6, scientific = FALSE)

  # ===== Print number of simulated samples
  cat("\n", "   Number of requested simulated samples = ", R.no)
  cat("\n", "   Number of completed simulated samples = ", b.no, rep("\n",2))


  cat("Simulated Defined Function Values", rep("\n", 2))

  # ==================================================================== #

  # Calculate estimated parameter from Dfunction
  list2env(as.list(estcoeff), envir = .GlobalEnv)
  estM  <- eval(parse(text=Sfunction))


  # Calculate Simulated parameter from Dfunction
  mcmc <- as.data.frame(mcmc)
  mcmc <- mcmc %>%
    mutate(abM = eval(parse(text=Sfunction)))
  abM <- mcmc[, "abM"]

  #### Confidence Intervals and p-value ####

  # Calculate Percentile Probability
  if (quantile(abM,probs=0.5)>0) {
    pM = 2*(sum(abM<0)/b.no)
  } else {
    pM = 2*(sum(abM>0)/b.no)
  }

  #### Percentile Confidence Intervals of Conditional Indirect Effects ####
  PCI <- matrix(1:8, nrow = 1, dimnames = list(c("        "),
                                             c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

  PCI[1,1] <- format(round(quantile(abM,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,2] <- format(round(quantile(abM,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,3] <- format(round(quantile(abM,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,4] <- format(round(estM, digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,5] <- format(round(quantile(abM,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,6] <- format(round(quantile(abM,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,7] <- format(round(quantile(abM,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,8] <- format(round(pM, digits = 4), nsmall = 4, scientific = FALSE)

  # Bias-Corrected Factor
  zM = qnorm(sum(abM<estM)/b.no)

  # Calculate Bias-Corrected Probability

  if ((estM>0 & min(abM)>0) | (estM<0 & max(abM)<0)) {
    pbM = 0
  } else if (qnorm(sum(abM>0)/b.no)+2*zM<0) {
    pbM = 2*pnorm((qnorm(sum(abM>0)/b.no)+2*zM))
  } else {
    pbM = 2*pnorm(-1*(qnorm(sum(abM>0)/b.no)+2*zM))
  }

  #### Bias-Corrected Confidence Intervals ####

  BCCI <- matrix(1:8, nrow = 1, dimnames = list(c("        "),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

  BCCI[1,1] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,2] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,3] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,4] <- format(round(estM, digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,5] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,6] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,7] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,8] <- format(round(pbM, digits = 4), nsmall = 4, scientific = FALSE)

  cat("\n")
  cat("Percentile Confidence Intervals for Sfunction", rep("\n", 2))
  rownames(PCI) <- rep("    ", nrow(PCI))
  print(PCI, quote=FALSE, right=TRUE)
  cat("\n")

  cat("\n")
  cat("Bias-Corrected Confidence Intervals for Sfunction", rep("\n", 2))
  rownames(BCCI) <- rep("    ", nrow(BCCI))
  print(BCCI, quote=FALSE, right=TRUE)
  cat("\n")

}  ## ===== End (function mccimm_lavaan_fun) ===== ##





## ====== Function "mccimm_mplus_fun" Monte Carlo Confidence Intervals for Defined Function (mplus) ====== ##
#' Monte Carlo Simulation for Confidence Intervals of Defined Function (mplus)
#'
#' Generate confidence intervals of a defined function from Mplus results using Monte Carlo simulation.
#' Location of estimated parameters can be found in mccimm::TECH1().
#'
#' @param mplus_output_file Mplus output (.out) file (output from Mplus).
#' @param results_file Mplus a text file (.txt) that saves the Mplus results (RESULTS IS "filename.txt" in Mplus SAVEDATA:).
#' @param a1 location of parameter a1 in Mplus Tech1 outputs.
#' @param a2 location of parameter a2 in Mplus Tech1 outputs.
#' @param a3 location of parameter a3 in Mplus Tech1 outputs.
#' @param a4 location of parameter a4 in Mplus Tech1 outputs.
#' @param z1 location of parameter z1 in Mplus Tech1 outputs.
#' @param z2 location of parameter z2 in Mplus Tech1 outputs.
#' @param z3 location of parameter z3 in Mplus Tech1 outputs.
#' @param z4 location of parameter z4 in Mplus Tech1 outputs.
#' @param w1 location of parameter w1 in Mplus Tech1 outputs.
#' @param w2 location of parameter w2 in Mplus Tech1 outputs.
#' @param w3 location of parameter w3 in Mplus Tech1 outputs.
#' @param w4 location of parameter w4 in Mplus Tech1 outputs.
#' @param zw1 location of parameter zw1 in Mplus Tech1 outputs.
#' @param zw2 location of parameter zw2 in Mplus Tech1 outputs.
#' @param zw3 location of parameter zw3 in Mplus Tech1 outputs.
#' @param zw4 location of parameter zw4 in Mplus Tech1 outputs.
#' @param Sfunction function of estimated parameters from Mplus outputs.
#' @param R number of Monte Carlo simulation samples (in millions). For example, R=5 (default) generates 5,000,000 simulated samples.
#'
#' @return confidence intervals of defined function
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' # mplus_output_file is "model cc4.out", results_file is "Model_CC4.txt" & moderator Z is "AUTO"
#' # output mccimm object is mcObject
#'
#' mcObject <- mccimm_mplus("model cc4.out", "Model_CC4.txt",
#'             a1 = "60", a2 = "65", Sfunction="a1*a2")
#'

mccimm_mplus_fun <- function(mplus_output_file = "mplus_output.out",
                             results_file = "results.txt",
                             a1="NA", a2="NA", a3="NA", a4="NA",
                             z1="NA", z2="NA", z3="NA", z4="NA",
                             w1="NA", w2="NA", w3="NA", w4="NA",
                             zw1="NA", zw2="NA", zw3="NA", zw4="NA",
                             Sfunction="NULL", R=5) {

  mplus_output <- readModels(mplus_output_file)
  results <- mplus_output$parameters$unstandardized
  temp <- scan(results_file, sep="")

  Temp3 <- mplus_output$tech3$paramCov
  Temp3[upper.tri(Temp3, diag = FALSE)] <- 0

  Tech3 <- Temp3 + t(Temp3)
  Tech3 <- Tech3 - diag(diag(Temp3))

  ## -- Extract defined parameters and vcov -- ##
  dp <- c(varZ, varW, a1, a2, a3, a4, z1, z2, z3, z4, w1, w2, w3, w4, zw1, zw2, zw3, zw4)

  dp_no <- suppressWarnings(as.numeric(dp))
  dp_list <- c("varZ", "varW", "a1", "a2", "a3", "a4", "z1", "z2", "z3", "z4", "w1", "w2", "w3", "w4", "zw1", "zw2", "zw3", "zw4")
  non_na_list <- dp_list[which(!is.na(dp_no))]
  dp <- dp[dp != "NA"]
  dp <- as.numeric(dp)

  estcoeff <- temp[dp]
  names(estcoeff) <- non_na_list

  Tech3 <- Tech3[dp, dp]
  rownames(Tech3) <- non_na_list
  colnames(Tech3) <- non_na_list

  # -- Reassigning variable names -- #
  if (varZ != "NA") varZ <- "varZ"
  if (varW != "NA") varW <- "varW"
  if (a1 != "NA") a1 <- "a1"
  if (a2 != "NA") a2 <- "a2"
  if (a3 != "NA") a3 <- "a3"
  if (a4 != "NA") a4 <- "a4"
  if (z1 != "NA") z1 <- "z1"
  if (z2 != "NA") z2 <- "z2"
  if (z3 != "NA") z3 <- "z3"
  if (z4 != "NA") z4 <- "z4"
  if (w1 != "NA") w1 <- "w1"
  if (w2 != "NA") w2 <- "w2"
  if (w3 != "NA") w3 <- "w3"
  if (w4 != "NA") w4 <- "w4"
  if (zw1 != "NA") zw1 <- "zw1"
  if (zw2 != "NA") zw2 <- "zw2"
  if (zw3 != "NA") zw3 <- "zw3"
  if (zw4 != "NA") zw4 <- "zw4"

  ## -- Monte Carlo Simulation of R*1e6 samples, default: R = 5 -- ##
  mcmc <- MASS::mvrnorm(n=R*1e6, mu=estcoeff, Sigma=Tech3, tol = 1e-6)

  b.no <- nrow(mcmc)
  R.no <- format(R*1e6, scientific = FALSE)

  # ===== Print number of simulated samples
  cat("\n", "   Number of requested simulated samples = ", R.no)
  cat("\n", "   Number of completed simulated samples = ", b.no, rep("\n",2))


  cat("Simulated Defined Function Values", rep("\n", 2))

  # ==================================================================== #

  # Calculate estimated parameter from Dfunction
  list2env(as.list(estcoeff), envir = .GlobalEnv)
  estM  <- eval(parse(text=Sfunction))


  # Calculate Simulated parameter from Dfunction
  mcmc <- as.data.frame(mcmc)
  mcmc <- mcmc %>%
    mutate(abM = eval(parse(text=Sfunction)))
  abM <- mcmc[, "abM"]

  #### Confidence Intervals and p-value ####

  # Calculate Percentile Probability
  if (quantile(abM,probs=0.5)>0) {
    pM = 2*(sum(abM<0)/b.no)
  } else {
    pM = 2*(sum(abM>0)/b.no)
  }

  #### Percentile Confidence Intervals of Conditional Indirect Effects ####
  PCI <- matrix(1:8, nrow = 1, dimnames = list(c("        "),
                                             c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

  PCI[1,1] <- format(round(quantile(abM,c(0.005)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,2] <- format(round(quantile(abM,c(0.025)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,3] <- format(round(quantile(abM,c(0.05)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,4] <- format(round(estM, digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,5] <- format(round(quantile(abM,c(0.95)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,6] <- format(round(quantile(abM,c(0.975)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,7] <- format(round(quantile(abM,c(0.995)), digits = 4), nsmall = 4, scientific = FALSE)
  PCI[1,8] <- format(round(pM, digits = 4), nsmall = 4, scientific = FALSE)

  # Bias-Corrected Factor
  zM = qnorm(sum(abM<estM)/b.no)

  # Calculate Bias-Corrected Probability

  if ((estM>0 & min(abM)>0) | (estM<0 & max(abM)<0)) {
    pbM = 0
  } else if (qnorm(sum(abM>0)/b.no)+2*zM<0) {
    pbM = 2*pnorm((qnorm(sum(abM>0)/b.no)+2*zM))
  } else {
    pbM = 2*pnorm(-1*(qnorm(sum(abM>0)/b.no)+2*zM))
  }

  #### Bias-Corrected Confidence Intervals ####

  BCCI <- matrix(1:8, nrow = 1, dimnames = list(c("        "),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

  BCCI[1,1] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.005))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,2] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.025))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,3] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.050))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,4] <- format(round(estM, digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,5] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.950))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,6] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.975))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,7] <- format(round(quantile(abM,probs=pnorm(2*zM+qnorm(0.995))), digits = 4), nsmall = 4, scientific = FALSE)
  BCCI[1,8] <- format(round(pbM, digits = 4), nsmall = 4, scientific = FALSE)

  cat("\n")
  cat("Percentile Confidence Intervals for Sfunction", rep("\n", 2))
  rownames(PCI) <- rep("    ", nrow(PCI))
  print(PCI, quote=FALSE, right=TRUE)
  cat("\n")

  cat("\n")
  cat("Bias-Corrected Confidence Intervals for Sfunction", rep("\n", 2))
  rownames(BCCI) <- rep("    ", nrow(BCCI))
  print(BCCI, quote=FALSE, right=TRUE)
  cat("\n")

}  ## ===== End (function mccimm_mplus_fun) ===== ##





## ===== FUNCTION JN_plot to plot Johnson-Neyman Figure ===== ##
#' Generate Johnson-Neyman Figure from mccimm object
#'
#' Generate Johnson-Neyman Figure from mccimm object
#'
#' @param mccimmObject mccimm object (output from mccimm_modsem or mccimm_mplus).
#' @param min_z minimum values of moderator Z on the graph.
#' @param max_z maximum values of moderator Z on the graph.
#' @param detail detail of the Johnson-Neyman graph. The higher the number, the smoother the lines, but it will also take longer to generate.
#' @param lower.quantile lower quantile for the confidence intervals. 0.025 for 95-percent confidence intervals (default). 0.005 for 99-percent confidence intervals.
#' @param upper.quantile upper quantile for the confidence intervals. 0.975 for 95-percent confidence intervals (default). 0.995 for 99-percent confidence intervals.
#' @param alpha alpha = 0.2 (default) Thickness of the lines.
#' @param sd.line horiznotal line that shows mean +/- 2 (default) values of the moderator.
#'
#' @return a figure for changing the titles.
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' # modsem object is mcObject & output JN_plot object is JN_figure
#'
#' JN_figure <- JN_plot(mcObject)
#'

JN_plot <- function (mccimmObject, ci="bc",
                     min_z = -3, max_z = 3, detail = 300,
                     lower.quantile = 0.025, upper.quantile = 0.975,
                     alpha = 0.2, sd.line = 2)
{

  cat("\n")
  cat("Generating Johnson-Neyman Figure ... ", rep("\n", 2))


  ## -- Get values from mccimmObject -- ##
  a1 <- mccimmObject$a1
  a2 <- mccimmObject$a2
  a3 <- mccimmObject$a3
  a4 <- mccimmObject$a4
  z1 <- mccimmObject$z1
  z2 <- mccimmObject$z2
  z3 <- mccimmObject$z3
  z4 <- mccimmObject$z4
  w1 <- mccimmObject$w1
  w2 <- mccimmObject$w2
  w3 <- mccimmObject$w3
  w4 <- mccimmObject$w4

  sd_z <- mccimmObject$sd_z

  Sa1 <- mccimmObject$Sa1
  Sa2 <- mccimmObject$Sa2
  Sa3 <- mccimmObject$Sa3
  Sa4 <- mccimmObject$Sa4
  Sz1 <- mccimmObject$Sz1
  Sz2 <- mccimmObject$Sz2
  Sz3 <- mccimmObject$Sz3
  Sz4 <- mccimmObject$Sz4
  Sw1 <- mccimmObject$Sw1
  Sw2 <- mccimmObject$Sw2
  Sw3 <- mccimmObject$Sw3
  Sw4 <- mccimmObject$Sw4

  b.no <- mccimmObject$b.no
  NoModz <- mccimmObject$NoModz
  NoModw <- mccimmObject$NoModw
  Z <- mccimmObject$Z


  ## -- Calculating Parameters -- ##
  mean_z <- 0 ## estX
  min_z_abs <- min_z + mean_z
  max_z_abs <- max_z + mean_z

  valsz <- seq(min_z_abs, max_z_abs, length.out = detail)
  if (NoModz == 1) {
    zestX <- (a1+z1*valsz)*(a2+z2*valsz)*(a3+z3*valsz)*(a4+z4*valsz)
  } else {
    zestX <- (a1+w1*valsz)*(a2+w2*valsz)*(a3+w3*valsz)*(a4+w4*valsz)
  }

  CI <- data.frame(betax = zestX, betax.lower = NA_real_, betax.upper = NA_real_, z = valsz)  # Initialize CI

  for (i in seq_len(detail)) {
    if (NoModz == 1) {
      zabX <- (Sa1+Sz1*valsz[i])*(Sa2+Sz2*valsz[i])*(Sa3+Sz3*valsz[i])*(Sa4+Sz4*valsz[i])
    } else {
      zabX <- (Sa1+Sw1*valsz[i])*(Sa2+Sw2*valsz[i])*(Sa3+Sw3*valsz[i])*(Sa4+Sw4*valsz[i])
    }

    if (ci == "percent") {
      CI[i, "betax.lower"] <- stats::quantile(zabX, probs = lower.quantile, na.rm = TRUE)  ##
      CI[i, "betax.upper"] <- stats::quantile(zabX, probs = upper.quantile, na.rm = TRUE)  ##
    } else {  ## ci == "bc"
      zX = qnorm(sum(zabX<zestX[i])/b.no)  # Bias-Corrected Factor
      CI[i, "betax.lower"] <- stats::quantile(zabX, probs = pnorm(2*zX+qnorm(lower.quantile)), na.rm = TRUE)  ##
      CI[i, "betax.upper"] <- stats::quantile(zabX, probs = pnorm(2*zX+qnorm(upper.quantile)), na.rm = TRUE)  ##
    }
  }

  CI$sig <- CI$betax.lower > 0 | CI$betax.upper < 0

  df_plot <- data.frame(z = CI$z, slope = as.numeric(CI$betax),
        lower_all = as.numeric(CI$betax.lower), upper_all = as.numeric(CI$betax.upper),
        significant = CI$sig)
  siglabel <- "CI excludes 0"
  significance_chr <- ifelse(df_plot$significant, siglabel, "n.s.")
  df_plot$Significance <- factor(significance_chr, levels = c(siglabel, "n.s."))
  df_plot$run_id <- cumsum(c(0, diff(as.integer(df_plot$significant)) != 0))
  x_start <- mean_z - sd.line * sd_z
  x_end <- mean_z + sd.line * sd_z
  if (x_start < min_z_abs && x_end > max_z_abs) {
    warning("Truncating SD-range on the right and left!")
  } else if (x_start < min_z_abs) {
    warning("Truncating SD-range on the left!")
  } else if (x_end > max_z_abs) {
    warning("Truncating SD-range on the right!")
  }
  x_start <- max(x_start, min_z_abs)
  x_end <- min(x_end, max_z_abs)
  y_start <- y_end <- 0
  hline_label <- sprintf("+/- %s SDs of %s", sd.line, Z)
  data_hline <- data.frame(x_start, x_end, y_start, y_end, hline_label)
  breaks <- c(siglabel, "n.s.", hline_label)
  values <- structure(c("cyan3", "red", "black"), names = breaks)
  y_range <- range(c(df_plot$lower_all, df_plot$upper_all, 0), na.rm = TRUE)
  if (!all(is.finite(y_range))) y_range <- c(-1, 1)
  flip_idx <- which(diff(as.integer(df_plot$significant)) != 0)
  approx_jn <- numeric(0)

  if (length(flip_idx) > 0) {
    for (k in flip_idx) {
      z0 <- df_plot$z[k]
      z1 <- df_plot$z[k + 1]
      lo0 <- df_plot$lower_all[k]
      lo1 <- df_plot$lower_all[k + 1]
      hi0 <- df_plot$upper_all[k]
      hi1 <- df_plot$upper_all[k + 1]
      use_lower <- xor(lo0 > 0, lo1 > 0)
      if (use_lower) {
        t <- (0 - lo0)/(lo1 - lo0)
      } else {
        t <- (0 - hi0)/(hi1 - hi0)
      }
      t <- min(max(t, 0), 1)
      approx_jn <- c(approx_jn, z0 + t * (z1 - z0))
    }
  }
  slope <- NULL
  lower_all <- NULL
  upper_all <- NULL
  Significance <- NULL
  run_id <- NULL


  ## -- Plot with ggplot2 -- ##
  p_jn <- ggplot2::ggplot(df_plot, ggplot2::aes(x = z, y = slope)) +
       ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_all, ymax = upper_all,
                     fill = Significance, group = run_id), alpha = alpha, na.rm = TRUE) +
       ggplot2::geom_line(ggplot2::aes(color = Significance,
                     group = run_id), linewidth = 1, na.rm = TRUE) +
       ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
       suppressWarnings(ggplot2::geom_segment(mapping = ggplot2::aes(x = x_start,
                     xend = x_end, y = y_start, yend = y_end, color = hline_label,
                     fill = hline_label), data = data_hline, linewidth = 1.5)) +
       ggplot2::ggtitle("Johnson-Neyman Figure of Moderated-Mediating Effect") +
       ggplot2::scale_discrete_manual(aesthetics = c("colour",
                    "fill"), name = "", values = values, breaks = breaks, drop = FALSE) +
       ggplot2::scale_y_continuous(limits = y_range) +
       ggplot2::labs(x = Z, y = paste("Indirect Effect")) +
       ggplot2::theme_minimal() +
       ggplot2::theme(plot.title = element_text(hjust = 0.5))
       if (length(approx_jn) > 0) {
         top_y <- suppressWarnings(max(df_plot$slope[is.finite(df_plot$slope)], na.rm = TRUE))
         if (!is.finite(top_y)) top_y <- y_range[2]
         for (zstar in approx_jn) {
           if (is.finite(zstar) && zstar >= min_z_abs && zstar <= max_z_abs) {
             p_jn <- p_jn + ggplot2::geom_vline(xintercept = zstar, linetype = "dashed", color = "red") +
                    ggplot2::annotate("text", x = zstar, y = top_y, label = paste("JN point (~):",
                              round(zstar, 2)), hjust = -0.1, vjust = 1, color = "black")
           }
         }
       }

  ggplot2::ggsave("Johnson-Neyman Figure.png", width = 22.86, height = 16.51, units = "cm")

  ## -- Print the Figure -- ##
  print(p_jn)
  return(p_jn)

  cat("\n")
  cat("Figure is p_jn and is saved as 'Johnson-Neyman Figure.png'", rep("\n", 2))


} ## ===== End (plot Johnson-Neyman Figure) ===== ##


## ===== Function TECH1 to identify position of reference parameters for Mplus outputs ===== ##
#' Present Mplus results with matching numbers on TECH1 outputs from Mplus
#'
#' Present Mplus estimated parameters with matching number on TECH1 outputs (for inputs in mccimm_mplus)
#'
#' @param mplusoutput.file Mplus output file (with TECH 1 on the OUTPUT: command in Mplus input file).
#'
#' @return outputs on screen
#' @export
#' @examples
#'
#' ## -- Example -- ##
#'
#' # Mplus output file "example_d3.out" & output TECH1 object is Tech1
#'
#' Tech1 <- TECH1("example_d3.out")
#'
TECH1 <- function(mplusoutput.file = "example_d3.out") {
  mplus_output <- MplusAutomation::readModels(mplusoutput.file)
  results <- mplus_output$parameters
  spec <<- mplus_output$tech1$parameterSpecification

  ## ----- One-level ----- ##
  if (ncol(results$unstandardized) == 6) {

    Tech1 <- matrix(nrow = 0, ncol = 4)
    colnames(Tech1) <- c("par", "Row_Name", "Col_Name", "Value")
    Tech1 <- data.frame(Tech1)

    nuW <- Tech1_1_level(spec, "nu")
    Tech1 <- rbind(Tech1, nuW)
    lambdaW <- Tech1_1_level(spec, "lambda")
    Tech1 <- rbind(Tech1, lambdaW)
    thetaW <- Tech1_1_level(spec, "theta")
    Tech1 <- rbind(Tech1, thetaW)
    alphaW <- Tech1_1_level(spec, "alpha")
    Tech1 <- rbind(Tech1, alphaW)
    betaW <- Tech1_1_level(spec, "beta")
    Tech1 <- rbind(Tech1, betaW)
    psiW <- Tech1_1_level(spec, "psi")
    Tech1 <- rbind(Tech1, psiW)

    Tech1 <- Tech1[order(Tech1$Value),]

    results <- mplus_output$parameters
    res.tech1 <- data.frame(lapply(results$unstandardized, unlist))
    res.tech1["tech1"] <- 0


    ## -- Match TECH1 information to parameters file -- ##
    for (i in 1:mplus_output$summaries$Parameters) {
      t.par <- Tech1[i, "par"]
      t.row <- Tech1[i, "Row_Name"]
      t.col <- Tech1[i, "Col_Name"]
      if (t.par == "psi" | t.par == "theta") {
        if (t.row != t.col) {
          lhs = paste0(t.row,".WITH")
          ind <- which(res.tech1[,"paramHeader"] == lhs & res.tech1[,"param"] == t.col)
        } else {
          ind <- which((res.tech1[,"paramHeader"] == "Variances" | res.tech1[,"paramHeader"] == "Residual.Variances") & res.tech1[,"param"] == t.col)
        }
      }
      if (t.par == "lambda") {
        lhs = paste0(t.row,".BY")
        ind <- which(res.tech1[,"paramHeader"] == lhs & res.tech1[,"param"] == t.col)
      }
      if (t.par == "alpha" | t.par == "nu") {
        ind <- which((res.tech1[,"paramHeader"] == "Means" | res.tech1[,"paramHeader"] == "Intercepts") & res.tech1[,"param"] == t.col)
      }
      if (t.par == "beta") {
        lhs = paste0(t.row,".ON")
        ind <- which(res.tech1[,"paramHeader"] == lhs & res.tech1[,"param"] == t.col)
      }
      res.tech1[ind,"tech1"] <- i
    }
    ## end (for i) ------- ##
    return(res.tech1)
  }

  ## ----- Two-levels ----- ##
  if (ncol(results$unstandardized) == 7) {
    Tech1 <- matrix(nrow = 0, ncol = 5)
    colnames(Tech1) <- c("level", "par", "Row_Name", "Col_Name", "Value")
    Tech1 <- data.frame(Tech1)

    nuW <- Tech1_2_level(spec, "Within", "nu")
    Tech1 <- rbind(Tech1, nuW)
    lambdaW <- Tech1_2_level(spec, "Within", "lambda")
    Tech1 <- rbind(Tech1, lambdaW)
    thetaW <- Tech1_2_level(spec, "Within", "theta")
    Tech1 <- rbind(Tech1, thetaW)
    alphaW <- Tech1_2_level(spec, "Within", "alpha")
    Tech1 <- rbind(Tech1, alphaW)
    betaW <- Tech1_2_level(spec, "Within", "beta")
    Tech1 <- rbind(Tech1, betaW)
    psiW <- Tech1_2_level(spec, "Within", "psi")
    Tech1 <- rbind(Tech1, psiW)

    nuB <- Tech1_2_level(spec, "Between", "nu")
    Tech1 <- rbind(Tech1, nuB)
    lambdaB <- Tech1_2_level(spec, "Between", "lambda")
    Tech1 <- rbind(Tech1, lambdaB)
    thetaB <- Tech1_2_level(spec, "Between", "theta")
    Tech1 <- rbind(Tech1, thetaB)
    alphaB <- Tech1_2_level(spec, "Between", "alpha")
    Tech1 <- rbind(Tech1, alphaB)
    betaB <- Tech1_2_level(spec, "Between", "beta")
    Tech1 <- rbind(Tech1, betaB)
    psiB <- Tech1_2_level(spec, "Between", "psi")
    Tech1 <- rbind(Tech1, psiB)

    Tech1 <- Tech1[order(Tech1$Value),]

    results <- mplus_output$parameters
    res.tech1 <- data.frame(lapply(results$unstandardized, unlist))
    res.tech1["tech1"] <- 0


    ## -- Match TECH1 information to parameters file -- ##
    for (i in 1:mplus_output$summaries$Parameters) {
      t.level <- Tech1[i,"level"]
      t.par <- Tech1[i, "par"]
      t.row <- Tech1[i, "Row_Name"]
      t.col <- Tech1[i, "Col_Name"]
      if (t.par == "psi" | t.par == "theta") {
        if (t.row != t.col) {
          lhs = paste0(t.row,".WITH")
          ind <- which(res.tech1[,"paramHeader"] == lhs & res.tech1[,"param"] == t.col & res.tech1[,"BetweenWithin"] == t.level)
        } else {
          ind <- which((res.tech1[,"paramHeader"] == "Variances" | res.tech1[,"paramHeader"] == "Residual.Variances") & res.tech1[,"param"] == t.col &
             res.tech1[,"BetweenWithin"] == t.level)
        }
      }
      if (t.par == "lambda") {
        lhs = paste0(t.row,".BY")
        ind <- which(res.tech1[,"paramHeader"] == lhs & res.tech1[,"param"] == t.col & res.tech1[,"BetweenWithin"] == t.level)
      }
      if (t.par == "alpha" | t.par == "nu") {
        ind <- which((res.tech1[,"paramHeader"] == "Means" | res.tech1[,"paramHeader"] == "Intercepts") & res.tech1[,"param"] == t.col &
             res.tech1[,"BetweenWithin"] == t.level)
      }
      if (t.par == "beta") {
        lhs = paste0(t.row,".ON")
        ind <- which(res.tech1[,"paramHeader"] == lhs & res.tech1[,"param"] == t.col & res.tech1[,"BetweenWithin"] == t.level)
      }
      res.tech1[ind,"tech1"] <- i
    }
    ## ------- ##

    return(res.tech1)
  }

  sink("Tech1.txt")
  print(res.tech1)
  sink()
  cat("\n")
  cat("Output is saved as 'Tech1.txt'", rep("\n", 2))

} # end (function TECH1)




## -- Sub-function identifying what the numbers in TECH1 refer to (One-Level) -- ##
Tech1_1_level <- function(spec, para) {
  if (names(spec)[1] != "X") {
    p.tech1 <- paste0("spec$", para)
  } else {
    p.tech1 <- paste0("spec$X$", para)
  }
  parameter <- eval(parse(text = p.tech1))
  parameter_greater_than_0 <- which(parameter > 0, arr.ind = TRUE)
  if (length(parameter_greater_than_0) > 0) {
    result_parameter <- data.frame(
      par = para,
      Row_Name = rownames(parameter)[parameter_greater_than_0[, "row"]],
      Col_Name = colnames(parameter)[parameter_greater_than_0[, "col"]],
      Value = parameter[parameter_greater_than_0]
    )
    return(result_parameter)
  }
} # end (Sub-Function "Tech1_1_level")

## -- Sub-function identifying what the numbers in TECH1 refer to (Two-Level) -- ##
Tech1_2_level <- function(spec, level, para) {
  p.tech1 <- paste0("spec$", toupper(level), "$", para)
  parameter <- eval(parse(text = p.tech1))
  parameter_greater_than_0 <- which(parameter > 0, arr.ind = TRUE)
  if (length(parameter_greater_than_0) > 0) {
    result_parameter <- data.frame(
      level = level,
      par = para,
      Row_Name = rownames(parameter)[parameter_greater_than_0[, "row"]],
      Col_Name = colnames(parameter)[parameter_greater_than_0[, "col"]],
      Value = parameter[parameter_greater_than_0]
    )
    return(result_parameter)
  }
} # end (Sub-Function "Tech1_2_level")




## ===== Sub_Function to Create 2-way Standardized Interaction Figure ===== ##

Two_Way_Figure <- function(estX.lo, estX.hi) {

  # -- Calculate endpoints of two lines-- #
  df_wide <- data.frame(
    x_var = c(-2,2),
    line1 = c(estX.lo*-2, estX.lo*2),
    line2 = c(estX.hi*-2, estX.hi*2)
  )

  # -- Convert to long format -- #
  df_long <- df_wide %>%
    pivot_longer(
      cols = starts_with("line"),
      names_to = "line_id",
      values_to = "y_value"
    )
    y.upper = 1.5*max(df_long$y_value) - 0.5*mean(df_long$y_value)
    y.lower = 1.5*min(df_long$y_value) - 0.5*mean(df_long$y_value)

  # -- Plot Figure p_int -- #
  p_int <<- ggplot2::ggplot(data = df_long, aes(x = x_var, y = y_value, color = line_id, linetype = line_id)) +
                    xlim(-2.5, 2.5) +
                    ylim(y.lower, y.upper) +
           geom_line(linewidth=1) +
           scale_linetype_manual(name = "Moderator Levels",
                                 values = c("line1" = "solid", "line2" = "solid"),
                                 labels = c("Mean - 1 S.D.", "Mean + 1 S.D.")) +
           scale_color_manual(name = "Moderator Levels",
                              values = c("line1" = "black", "line2" = "grey"),
                              labels = c("Mean - 1 S.D.", "Mean + 1 S.D.")) +
           guides(linetype = "none") +
           labs(title = "Standardized Moderated-Mediating Effects",
                    x = "X-axis Label",
                    y = "Mediating Effects of X on Y through M") +
           theme_classic()  # Optional: Apply a theme

  print(p_int) # show p_int

  ggplot2::ggsave("2-Way Standardized Interaction Figure.png", width = 22.86, height = 16.51, units = "cm")

  cat("\n")
  cat("Figure is p_int and is saved as '2-Way standardized Interaction Figure.png'", rep("\n", 2))

} # end (Sub-Function: Two-Way_Figure)

## ===== Close 2-Way Standardized Interaction Figure ===== ##
