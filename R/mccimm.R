#######################################################################################
# Author: Prof. Gordon Cheung (University of Auckland)                                #
# email: gordon.cheung@auckland.ac.nz                                                 #
# Purpose: Monte Carlo simulation of conditional indirect effects from modsem outputs #
# Model: Allows 3-way interaction                                                     #
# Allows simulated function (Sfunction = "    "                                       #
# Version: Beta 1.3.1                                                                 #
# WARNING: This is beta version. Please report all bugs to the author                 #
#######################################################################################

## == Required packages (load when loading mccimm_modsem) == ##
library(MASS)
library(ggplot2)
# library (openxlsx)


# ==================== Creating Function "mccimm_modsem" Monte Carlo Confidence Intervals for Moderated Mediation (modsem) ==================== #
mccimm_modsem <- function(object = est_lms, Z="NA", W="NA",
                   A1="NA", Z1="NA", W1="NA", ZW1="NA",
                   A2="NA", Z2="NA", W2="NA", ZW2="NA",
                   A3="NA", Z3="NA", W3="NA", ZW3="NA",
                   A4="NA", Z4="NA", W4="NA", ZW4="NA",
                   R=5) {

  ## --- Initial Inputs for programming --- ##

  # setwd("C:/Research/mccimm/mccimm_modsem")  ## set working directory

#  object <- est_lms
#  Z <- "Autonomy"          ## Specify moderator Z
#  W <- "NA"
#  A1 <- "a1"
#  A2 <- "a2"
#  A3 <- "a3a"
#  A4 <- "NA"
#  Z1 <- "z1"
#  Z2 <- "NA"
#  Z3 <- "NA"
#  Z4 <- "NA"
#  W1 <- "NA"
#  W2 <- "NA"
#  W3 <- "NA"
#  W4 <- "NA"
#  ZW1 <- "NA"
#  ZW2 <- "NA"
#  ZW3 <- "NA"
#  ZW4 <- "NA"

#  R <- 1 # Number of simulated samples = R*1e6 (default: R = 5)
  ## ------------------------------- ##


  ## Extract defined parameters and vcov ##
  varZ <- "NA"
  varW <- "NA"
  if (Z != "NA") varZ <- paste0(Z, "~~", Z)
  if (W != "NA") varW <- paste0(W, "~~", W)
  dp <- c(A1, A2, A3, A4, Z1, Z2, Z3, Z4, W1, W2, W3, W4, ZW1, ZW2, ZW3, ZW4, varZ, varW)
  dp <- dp[dp != "NA"]


  temp <- modsem_coef(object)
  estcoeff <- temp[dp]
  Temp3 <- modsem_vcov(object)
  Tech3 <- Temp3[dp, dp]


  ## -- Number of Moderating Effects -- ##
  NoModz <- 0
  NoModw <- 0
  if (Z1 != "NA") NoModz <- NoModz + 1
  if (Z2 != "NA") NoModz <- NoModz + 1
  if (Z3 != "NA") NoModz <- NoModz + 1
  if (Z4 != "NA") NoModz <- NoModz + 1
  if (W1 != "NA") NoModw <- NoModw + 1
  if (W2 != "NA") NoModw <- NoModw + 1
  if (W3 != "NA") NoModw <- NoModw + 1
  if (W4 != "NA") NoModw <- NoModw + 1

  NoMod <- NoModz + NoModw
  ## ----- ##

  ## -- Location of Moderator for calculation of Index MM -- ##
  if (NoMod == 1) {
    if (Z1 != "NA") PoMod <- 1
    if (Z2 != "NA") PoMod <- 2
    if (Z3 != "NA") PoMod <- 3
    if (Z4 != "NA") PoMod <- 4
    if (W1 != "NA") PoMod <- 1
    if (W2 != "NA") PoMod <- 2
    if (W3 != "NA") PoMod <- 3
    if (W4 != "NA") PoMod <- 4
  }

  ## -- Initialize a-paths to 1 -- ##
  a1 <- 1
  a2 <- 1
  a3 <- 1
  a4 <- 1
  Sa1 <- matrix(1, R*1e6)
  Sa2 <- matrix(1, R*1e6)
  Sa3 <- matrix(1, R*1e6)
  Sa4 <- matrix(1, R*1e6)
  ## ------------------------------ ##

  ## -- Initialize interaction paths to 0 -- ##
  z1 <- 0
  z2 <- 0
  z3 <- 0
  z4 <- 0
  Sz1 <- matrix(0, R*1e6)
  Sz2 <- matrix(0, R*1e6)
  Sz3 <- matrix(0, R*1e6)
  Sz4 <- matrix(0, R*1e6)

  w1 <- 0
  w2 <- 0
  w3 <- 0
  w4 <- 0
  Sw1 <- matrix(0, R*1e6)
  Sw2 <- matrix(0, R*1e6)
  Sw3 <- matrix(0, R*1e6)
  Sw4 <- matrix(0, R*1e6)

  zw1 <- 0
  zw2 <- 0
  zw3 <- 0
  zw4 <- 0
  Szw1 <- matrix(0, R*1e6)
  Szw2 <- matrix(0, R*1e6)
  Szw3 <- matrix(0, R*1e6)
  Szw4 <- matrix(0, R*1e6)

  ## -- Initialize stdZ and stdW -- ##
  stdZ <- 1
  stdW <- 1
  SstdZ <- matrix(1, R*1e6)
  SstdW <- matrix(1, R*1e6)

  ## ------------------------------ ##

  ## -- Monte Carlo Simulation of R*1e6 samples, default: R = 5 -- ##
  mcmc <- MASS::mvrnorm(n=R*1e6, mu=estcoeff, Sigma=Tech3, tol = 1e-6)


  # ===== Retain simulated samples with variance larger than or equal to 0
  if (NoMod == 1) {
    if (NoModz != 0){
      mcmc <- mcmc[which(mcmc[, varZ] >= 0),]
    } else {
      mcmc <- mcmc[which(mcmc[, varW] >= 0),]
    }
  } else if (NoMod == 2) {
    mcmc <- mcmc[which(mcmc[, varZ] >= 0),]
    mcmc <- mcmc[which(mcmc[, varW] >= 0),]
  }

  b.no <- nrow(mcmc)
  R.no <- format(R*1e6, scientific = FALSE)

  # ===== Print number of bootstrap samples
  cat("\n", "   Number of requested simulated sample = ", R.no)
  cat("\n", "   Number of completed simulated sample = ", b.no, rep("\n",2))

  # ==================================================================== #


  ### --- No Moderating Effect --- ###
  if (NoMod == 0) {

    # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% A1)) a1 <- estcoeff[A1]
    if (any(names(estcoeff) %in% A2)) a2 <- estcoeff[A2]
    if (any(names(estcoeff) %in% A3)) a3 <- estcoeff[A3]
    if (any(names(estcoeff) %in% A4)) a4 <- estcoeff[A4]

    # Calculate Estimated Indirect Effect
    estM  <- a1*a2*a3*a4

    # Capture simulated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% A1)) Sa1 <- mcmc[, A1]
    if (any(names(estcoeff) %in% A2)) Sa2 <- mcmc[, A2]
    if (any(names(estcoeff) %in% A3)) Sa3 <- mcmc[, A3]
    if (any(names(estcoeff) %in% A4)) Sa4 <- mcmc[, A4]

    # Calculate Simulated Indirect Effect
    abM <- Sa1*Sa2*Sa3*Sa4
    ##################################


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
    cat("Percentile Confidence Intervals for Indirect Effect", rep("\n", 2))
    rownames(PCI) <- rep("    ", nrow(PCI))
    print(PCI, quote=FALSE, right=TRUE)
    cat("\n")

    cat("\n")
    cat("Bias-Corrected Confidence Intervals for Indirect Effect", rep("\n", 2))
    rownames(BCCI) <- rep("    ", nrow(BCCI))
    print(BCCI, quote=FALSE, right=TRUE)
    cat("\n")

  }
  ### --- End No Moderating Effect --- ###


  ### --- Only One Moderating Effect Z --- ###
  if ((NoMod == 1) & (NoModz == 1)) {

    # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% A1)) a1 <- estcoeff[A1]
    if (any(names(estcoeff) %in% A2)) a2 <- estcoeff[A2]
    if (any(names(estcoeff) %in% A3)) a3 <- estcoeff[A3]
    if (any(names(estcoeff) %in% A4)) a4 <- estcoeff[A4]
    if (any(names(estcoeff) %in% Z1)) z1 <- estcoeff[Z1]
    if (any(names(estcoeff) %in% Z2)) z2 <- estcoeff[Z2]
    if (any(names(estcoeff) %in% Z3)) z3 <- estcoeff[Z3]
    if (any(names(estcoeff) %in% Z4)) z4 <- estcoeff[Z4]
    stdZ <- sqrt(estcoeff[varZ])

    # Calculate Estimated Index MM #
    if (PoMod == 1) estIMM <- z1*a2*a3*a4
    if (PoMod == 2) estIMM <- a1*z2*a3*a4
    if (PoMod == 3) estIMM <- a1*a2*z3*a4
    if (PoMod == 4) estIMM <- a1*a2*a3*z4

    # Capture simulated parameters for calculating indirect effects #
    if (any(names(estcoeff) %in% A1)) Sa1 <- mcmc[, A1]
    if (any(names(estcoeff) %in% A2)) Sa2 <- mcmc[, A2]
    if (any(names(estcoeff) %in% A3)) Sa3 <- mcmc[, A3]
    if (any(names(estcoeff) %in% A4)) Sa4 <- mcmc[, A4]
    if (any(names(estcoeff) %in% Z1)) Sz1 <- mcmc[, Z1]
    if (any(names(estcoeff) %in% Z2)) Sz2 <- mcmc[, Z2]
    if (any(names(estcoeff) %in% Z3)) Sz3 <- mcmc[, Z3]
    if (any(names(estcoeff) %in% Z4)) Sz4 <- mcmc[, Z4]
    SstdZ <- sqrt(mcmc[,varZ])

    # Calculate Simulated Index MM
    if (any(names(estcoeff) %in% Z1)) IndexMM <- Sz1*Sa2*Sa3*Sa4
    if (any(names(estcoeff) %in% Z2)) IndexMM <- Sa1*Sz2*Sa3*Sa4
    if (any(names(estcoeff) %in% Z3)) IndexMM <- Sa1*Sa2*Sz3*Sa4
    if (any(names(estcoeff) %in% Z4)) IndexMM <- Sa1*Sa2*Sa3*Sz4

    #### Percentile and Bias-Corrected Confidence Intervals of Conditional Indirect Effects ####

    PCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))
    BCCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                               c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    for (X in 1:5) {
      level <- X - 3
      estX <- (a1+z1*level*stdZ)*(a2+z2*level*stdZ)*(a3+z3*level*stdZ)*(a4+z4*level*stdZ)
      abX <- (Sa1+Sz1*level*SstdZ)*(Sa2+Sz2*level*SstdZ)*(Sa3+Sz3*level*SstdZ)*(Sa4+Sz4*level*SstdZ)
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

  } # End one moderating effect for Z


  ## --- Only One Moderating Effect W --- ##
  if ((NoMod == 1) & (NoModw == 1)) {

    # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% A1)) a1 <- estcoeff[A1]
    if (any(names(estcoeff) %in% A2)) a2 <- estcoeff[A2]
    if (any(names(estcoeff) %in% A3)) a3 <- estcoeff[A3]
    if (any(names(estcoeff) %in% A4)) a4 <- estcoeff[A4]
    if (any(names(estcoeff) %in% W1)) w1 <- estcoeff[W1]
    if (any(names(estcoeff) %in% W2)) w2 <- estcoeff[W2]
    if (any(names(estcoeff) %in% W3)) w3 <- estcoeff[W3]
    if (any(names(estcoeff) %in% W4)) w4 <- estcoeff[W4]
    stdW <- sqrt(estcoeff[varW])


    # Calculate Estimated Index MM #
    if (PoMod == 1) estIMM <- w1*a2*a3*a4
    if (PoMod == 2) estIMM <- a1*w2*a3*a4
    if (PoMod == 3) estIMM <- a1*a2*w3*a4
    if (PoMod == 4) estIMM <- a1*a2*a3*w4

    # Capture simulated parameters for calculating indirect effects #
    if (any(names(estcoeff) %in% A1)) Sa1 <- mcmc[, A1]
    if (any(names(estcoeff) %in% A2)) Sa2 <- mcmc[, A2]
    if (any(names(estcoeff) %in% A3)) Sa3 <- mcmc[, A3]
    if (any(names(estcoeff) %in% A4)) Sa4 <- mcmc[, A4]
    if (any(names(estcoeff) %in% W1)) Sw1 <- mcmc[, W1]
    if (any(names(estcoeff) %in% W2)) Sw2 <- mcmc[, W2]
    if (any(names(estcoeff) %in% W3)) Sw3 <- mcmc[, W3]
    if (any(names(estcoeff) %in% W4)) Sw4 <- mcmc[, W4]
    SstdW <- sqrt(mcmc[,varW])

    # Calculate Simulated Index MM
    if (any(names(estcoeff) %in% W1)) IndexMM <- Sw1*Sa2*Sa3*Sa4
    if (any(names(estcoeff) %in% W2)) IndexMM <- Sa1*Sw2*Sa3*Sa4
    if (any(names(estcoeff) %in% W3)) IndexMM <- Sa1*Sa2*Sw3*Sa4
    if (any(names(estcoeff) %in% W4)) IndexMM <- Sa1*Sa2*Sa3*Sw4


    #### Percentile and Bias-Corrected Confidence Intervals of Conditional Indirect Effects ####

    PCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))
    BCCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                               c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    for (X in 1:5) {
      level <- X - 3
      estX <- (a1+w1*level*stdW)*(a2+w2*level*stdW)*(a3+w3*level*stdW)*(a4+w4*level*stdW)
      abX <- (Sa1+Sw1*level*SstdW)*(Sa2+Sw2*level*SstdW)*(Sa3+Sw3*level*SstdW)*(Sa4+Sw4*level*SstdW)
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

  } # End one moderating effect for W



  ## --- More Than One Moderating Effect for Z only --- ##
  if ((NoMod > 1) & (NoModw == 0)) {

    # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% A1)) a1 <- estcoeff[A1]
    if (any(names(estcoeff) %in% A2)) a2 <- estcoeff[A2]
    if (any(names(estcoeff) %in% A3)) a3 <- estcoeff[A3]
    if (any(names(estcoeff) %in% A4)) a4 <- estcoeff[A4]
    if (any(names(estcoeff) %in% Z1)) z1 <- estcoeff[Z1]
    if (any(names(estcoeff) %in% Z2)) z2 <- estcoeff[Z2]
    if (any(names(estcoeff) %in% Z3)) z3 <- estcoeff[Z3]
    if (any(names(estcoeff) %in% Z4)) z4 <- estcoeff[Z4]
    stdZ <- sqrt(estcoeff[varZ])

    # Capture simulated parameters for calculating indirect effects #
    if (any(names(estcoeff) %in% A1)) Sa1 <- mcmc[, A1]
    if (any(names(estcoeff) %in% A2)) Sa2 <- mcmc[, A2]
    if (any(names(estcoeff) %in% A3)) Sa3 <- mcmc[, A3]
    if (any(names(estcoeff) %in% A4)) Sa4 <- mcmc[, A4]
    if (any(names(estcoeff) %in% Z1)) Sz1 <- mcmc[, Z1]
    if (any(names(estcoeff) %in% Z2)) Sz2 <- mcmc[, Z2]
    if (any(names(estcoeff) %in% Z3)) Sz3 <- mcmc[, Z3]
    if (any(names(estcoeff) %in% Z4)) Sz4 <- mcmc[, Z4]
    SstdZ <- sqrt(mcmc[,varZ])

    #### Percentile and Bias-Corrected Confidence Intervals of Conditional Indirect Effects ####

    PCI <- matrix(1:40, nrow = 5, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))
    BCCI <- matrix(1:40, nrow = 5, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"),
                                               c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    for (X in 1:5) {
      level <- X - 3
      estX <- (a1+z1*level*stdZ)*(a2+z2*level*stdZ)*(a3+z3*level*stdZ)*(a4+z4*level*stdZ)
      abX <- (Sa1+Sz1*level*SstdZ)*(Sa2+Sz2*level*SstdZ)*(Sa3+Sz3*level*SstdZ)*(Sa4+Sz4*level*SstdZ)
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

  } # End more than one moderating effect for Z


  ## --- More Than One Moderating Effect for W only --- ##
  if ((NoMod > 1) & (NoModz == 0)) {

    # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% A1)) a1 <- estcoeff[A1]
    if (any(names(estcoeff) %in% A2)) a2 <- estcoeff[A2]
    if (any(names(estcoeff) %in% A3)) a3 <- estcoeff[A3]
    if (any(names(estcoeff) %in% A4)) a4 <- estcoeff[A4]
    if (any(names(estcoeff) %in% W1)) w1 <- estcoeff[W1]
    if (any(names(estcoeff) %in% W2)) w2 <- estcoeff[W2]
    if (any(names(estcoeff) %in% W3)) w3 <- estcoeff[W3]
    if (any(names(estcoeff) %in% W4)) w4 <- estcoeff[W4]
    stdW <- sqrt(estcoeff[varW])

    # Capture simulated parameters for calculating indirect effects #
    if (any(names(estcoeff) %in% A1)) Sa1 <- mcmc[, A1]
    if (any(names(estcoeff) %in% A2)) Sa2 <- mcmc[, A2]
    if (any(names(estcoeff) %in% A3)) Sa3 <- mcmc[, A3]
    if (any(names(estcoeff) %in% A4)) Sa4 <- mcmc[, A4]
    if (any(names(estcoeff) %in% W1)) Sw1 <- mcmc[, W1]
    if (any(names(estcoeff) %in% W2)) Sw2 <- mcmc[, W2]
    if (any(names(estcoeff) %in% W3)) Sw3 <- mcmc[, W3]
    if (any(names(estcoeff) %in% W4)) Sw4 <- mcmc[, W4]
    SstdW <- sqrt(mcmc[,varW])

    #### Percentile and Bias-Corrected Confidence Intervals of Conditional Indirect Effects ####

    PCI <- matrix(1:40, nrow = 5, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))
    BCCI <- matrix(1:40, nrow = 5, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd"),
                                              c("     0.5%","     2.5%","       5%"," Estimate","      95%","    97.5%","    99.5%", "  p-value")))

    for (X in 1:5) {
      level <- X - 3
      estX <- (a1+w1*level*stdW)*(a2+w2*level*stdW)*(a3+w3*level*stdW)*(a4+w4*level*stdW)
      abX <- (Sa1+Sw1*level*SstdW)*(Sa2+Sw2*level*SstdW)*(Sa3+Sw3*level*SstdW)*(Sa4+Sw4*level*SstdW)
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


  } # End more than one moderating effect for W



  ## --- More Than One Moderating Effect for Z & W --- ##
  if ((NoMod > 1) & (NoModz != 0) & (NoModw != 0)) {

  # Define estimated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% A1)) a1 <- estcoeff[A1]
    if (any(names(estcoeff) %in% A2)) a2 <- estcoeff[A2]
    if (any(names(estcoeff) %in% A3)) a3 <- estcoeff[A3]
    if (any(names(estcoeff) %in% A4)) a4 <- estcoeff[A4]
    if (any(names(estcoeff) %in% Z1)) z1 <- estcoeff[Z1]
    if (any(names(estcoeff) %in% Z2)) z2 <- estcoeff[Z2]
    if (any(names(estcoeff) %in% Z3)) z3 <- estcoeff[Z3]
    if (any(names(estcoeff) %in% Z4)) z4 <- estcoeff[Z4]
    if (any(names(estcoeff) %in% W1)) w1 <- estcoeff[W1]
    if (any(names(estcoeff) %in% W2)) w2 <- estcoeff[W2]
    if (any(names(estcoeff) %in% W3)) w3 <- estcoeff[W3]
    if (any(names(estcoeff) %in% W4)) w4 <- estcoeff[W4]
    if (any(names(estcoeff) %in% ZW1)) zw1 <- estcoeff[ZW1]
    if (any(names(estcoeff) %in% ZW2)) zw2 <- estcoeff[ZW2]
    if (any(names(estcoeff) %in% ZW3)) zw3 <- estcoeff[ZW3]
    if (any(names(estcoeff) %in% ZW4)) zw4 <- estcoeff[ZW4]
    stdZ <- sqrt(estcoeff[varZ])
    stdW <- sqrt(estcoeff[varW])


  # Capture simulated parameters for calculating indirect effects
    if (any(names(estcoeff) %in% A1)) Sa1 <- mcmc[, A1]
    if (any(names(estcoeff) %in% A2)) Sa2 <- mcmc[, A2]
    if (any(names(estcoeff) %in% A3)) Sa3 <- mcmc[, A3]
    if (any(names(estcoeff) %in% A4)) Sa4 <- mcmc[, A4]
    if (any(names(estcoeff) %in% Z1)) Sz1 <- mcmc[, Z1]
    if (any(names(estcoeff) %in% Z2)) Sz2 <- mcmc[, Z2]
    if (any(names(estcoeff) %in% Z3)) Sz3 <- mcmc[, Z3]
    if (any(names(estcoeff) %in% Z4)) Sz4 <- mcmc[, Z4]
    if (any(names(estcoeff) %in% W1)) Sw1 <- mcmc[, W1]
    if (any(names(estcoeff) %in% W2)) Sw2 <- mcmc[, W2]
    if (any(names(estcoeff) %in% W3)) Sw3 <- mcmc[, W3]
    if (any(names(estcoeff) %in% W4)) Sw4 <- mcmc[, W4]
    if (any(names(estcoeff) %in% ZW1)) Szw1 <- mcmc[, ZW1]
    if (any(names(estcoeff) %in% ZW2)) Szw2 <- mcmc[, ZW2]
    if (any(names(estcoeff) %in% ZW3)) Szw3 <- mcmc[, ZW3]
    if (any(names(estcoeff) %in% ZW4)) Szw4 <- mcmc[, ZW4]
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
        estX <- (a1+z1*levelZ*stdZ+w1*levelW*stdW+zw1*levelZ*stdZ*levelW*stdW)*
                (a2+z2*levelZ*stdZ+w2*levelW*stdW+zw2*levelZ*stdZ*levelW*stdW)*
                (a3+z3*levelZ*stdZ+w3*levelW*stdW+zw3*levelZ*stdZ*levelW*stdW)*
                (a4+z4*levelZ*stdZ+w4*levelW*stdW+zw4*levelZ*stdZ*levelW*stdW)
        abX <- (Sa1+Sz1*levelZ*SstdZ+Sw1*levelW*SstdW+Szw1*levelZ*SstdZ*levelW*SstdW)*
               (Sa2+Sz2*levelZ*SstdZ+Sw2*levelW*SstdW+Szw2*levelZ*SstdZ*levelW*SstdW)*
               (Sa3+Sz3*levelZ*SstdZ+Sw3*levelW*SstdW+Szw3*levelZ*SstdZ*levelW*SstdW)*
               (Sa4+Sz4*levelZ*SstdZ+Sw4*levelW*SstdW+Szw4*levelZ*SstdZ*levelW*SstdW)
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
    estX <- (a1+z1*stdZ+w1*stdW+zw1*stdZ*stdW)*(a2+z2*stdZ+w2*stdW+zw2*stdZ*stdW)*(a3+z3*stdZ+w3*stdW+zw3*stdZ*stdW)*(a4+z4*stdZ+w4*stdW+zw4*stdZ*stdW) -
          (a1+z1*stdZ-w1*stdW-zw1*stdZ*stdW)*(a2+z2*stdZ-w2*stdW-zw2*stdZ*stdW)*(a3+z3*stdZ-w3*stdW-zw3*stdZ*stdW)*(a4+z4*stdZ-w4*stdW-zw4*stdZ*stdW)
    abX <- (Sa1+Sz1*SstdZ+Sw1*SstdW+Szw1*SstdZ*SstdW)*(Sa2+Sz2*SstdZ+Sw2*SstdW+Szw2*SstdZ*SstdW)*(Sa3+Sz3*SstdZ+Sw3*SstdW+Szw3*SstdZ*SstdW)*
           (Sa4+Sz4*SstdZ+Sw4*SstdW+Szw4*SstdZ*SstdW) -
           (Sa1+Sz1*SstdZ-Sw1*SstdW-Szw1*SstdZ*SstdW)*(Sa2+Sz2*SstdZ-Sw2*SstdW-Szw2*SstdZ*SstdW)*(Sa3+Sz3*SstdZ-Sw3*SstdW-Szw3*SstdZ*SstdW)*
           (Sa4+Sz4*SstdZ-Sw4*SstdW-Szw4*SstdZ*SstdW)
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
    estX <- (a1+z1*stdZ+w1*stdW+zw1*stdZ*stdW)*(a2+z2*stdZ+w2*stdW+zw2*stdZ*stdW)*(a3+z3*stdZ+w3*stdW+zw3*stdZ*stdW)*(a4+z4*stdZ+w4*stdW+zw4*stdZ*stdW) -
            (a1-z1*stdZ+w1*stdW-zw1*stdZ*stdW)*(a2-z2*stdZ+w2*stdW-zw2*stdZ*stdW)*(a3-z3*stdZ+w3*stdW-zw3*stdZ*stdW)*(a4-z4*stdZ+w4*stdW-zw4*stdZ*stdW)
    abX <- (Sa1+Sz1*SstdZ+Sw1*SstdW+Szw1*SstdZ*SstdW)*(Sa2+Sz2*SstdZ+Sw2*SstdW+Szw2*SstdZ*SstdW)*(Sa3+Sz3*SstdZ+Sw3*SstdW+Szw3*SstdZ*SstdW)*
           (Sa4+Sz4*SstdZ+Sw4*SstdW+Szw4*SstdZ*SstdW) -
           (Sa1-Sz1*SstdZ+Sw1*SstdW-Szw1*SstdZ*SstdW)*(Sa2-Sz2*SstdZ+Sw2*SstdW-Szw2*SstdZ*SstdW)*(Sa3-Sz3*SstdZ+Sw3*SstdW-Szw3*SstdZ*SstdW)*
           (Sa4-Sz4*SstdZ+Sw4*SstdW-Szw4*SstdZ*SstdW)
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
    estX <- (a1+z1*stdZ-w1*stdW-zw1*stdZ*stdW)*(a2+z2*stdZ-w2*stdW-zw2*stdZ*stdW)*(a3+z3*stdZ-w3*stdW-zw3*stdZ*stdW)*(a4+z4*stdZ-w4*stdW-zw4*stdZ*stdW) -
            (a1-z1*stdZ-w1*stdW+zw1*stdZ*stdW)*(a2-z2*stdZ-w2*stdW+zw2*stdZ*stdW)*(a3-z3*stdZ-w3*stdW+zw3*stdZ*stdW)*(a4-z4*stdZ-w4*stdW+zw4*stdZ*stdW)
    abX <- (Sa1+Sz1*SstdZ-Sw1*SstdW-Szw1*SstdZ*SstdW)*(Sa2+Sz2*SstdZ-Sw2*SstdW-Szw2*SstdZ*SstdW)*(Sa3+Sz3*SstdZ-Sw3*SstdW-Szw3*SstdZ*SstdW)*
           (Sa4+Sz4*SstdZ-Sw4*SstdW-Szw4*SstdZ*SstdW) -
           (Sa1-Sz1*SstdZ-Sw1*SstdW+Szw1*SstdZ*SstdW)*(Sa2-Sz2*SstdZ-Sw2*SstdW+Szw2*SstdZ*SstdW)*(Sa3-Sz3*SstdZ-Sw3*SstdW+Szw3*SstdZ*SstdW)*
           (Sa4-Sz4*SstdZ-Sw4*SstdW+Szw4*SstdZ*SstdW)
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
    estX <- (a1-z1*stdZ+w1*stdW-zw1*stdZ*stdW)*(a2-z2*stdZ+w2*stdW-zw2*stdZ*stdW)*(a3-z3*stdZ+w3*stdW-zw3*stdZ*stdW)*(a4-z4*stdZ+w4*stdW-zw4*stdZ*stdW) -
            (a1-z1*stdZ-w1*stdW+zw1*stdZ*stdW)*(a2-z2*stdZ-w2*stdW+zw2*stdZ*stdW)*(a3-z3*stdZ-w3*stdW+zw3*stdZ*stdW)*(a4-z4*stdZ-w4*stdW+zw4*stdZ*stdW)
    abX <- (Sa1-Sz1*SstdZ+Sw1*SstdW-Szw1*SstdZ*SstdW)*(Sa2-Sz2*SstdZ+Sw2*SstdW-Szw2*SstdZ*SstdW)*(Sa3-Sz3*SstdZ+Sw3*SstdW-Szw3*SstdZ*SstdW)*
           (Sa4-Sz4*SstdZ+Sw4*SstdW-Szw4*SstdZ*SstdW) -
           (Sa1-Sz1*SstdZ-Sw1*SstdW+Szw1*SstdZ*SstdW)*(Sa2-Sz2*SstdZ-Sw2*SstdW+Szw2*SstdZ*SstdW)*(Sa3-Sz3*SstdZ-Sw3*SstdW+Szw3*SstdZ*SstdW)*
           (Sa4-Sz4*SstdZ-Sw4*SstdW+Szw4*SstdZ*SstdW)
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
    estX <- (a1+z1*stdZ+w1*stdW+zw1*stdZ*stdW)*(a2+z2*stdZ+w2*stdW+zw2*stdZ*stdW)*(a3+z3*stdZ+w3*stdW+zw3*stdZ*stdW)*(a4+z4*stdZ+w4*stdW+zw4*stdZ*stdW) -
            (a1-z1*stdZ-w1*stdW+zw1*stdZ*stdW)*(a2-z2*stdZ-w2*stdW+zw2*stdZ*stdW)*(a3-z3*stdZ-w3*stdW+zw3*stdZ*stdW)*(a4-z4*stdZ-w4*stdW+zw4*stdZ*stdW)
    abX <- (Sa1+Sz1*SstdZ+Sw1*SstdW+Szw1*SstdZ*SstdW)*(Sa2+Sz2*SstdZ+Sw2*SstdW+Szw2*SstdZ*SstdW)*(Sa3+Sz3*SstdZ+Sw3*SstdW+Szw3*SstdZ*SstdW)*
           (Sa4+Sz4*SstdZ+Sw4*SstdW+Szw4*SstdZ*SstdW) -
           (Sa1-Sz1*SstdZ-Sw1*SstdW+Szw1*SstdZ*SstdW)*(Sa2-Sz2*SstdZ-Sw2*SstdW+Szw2*SstdZ*SstdW)*(Sa3-Sz3*SstdZ-Sw3*SstdW+Szw3*SstdZ*SstdW)*
           (Sa4-Sz4*SstdZ-Sw4*SstdW+Szw4*SstdZ*SstdW)
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

    ## Hi/Hi ##
    ## estX <- (a1+z1*stdZ+w1*stdW+zw1*stdZ*stdW)*(a2+z2*stdZ+w2*stdW+zw2*stdZ*stdW)*(a3+z3*stdZ+w3*stdW+zw3*stdZ*stdW)*(a4+z4*stdZ+w4*stdW+zw4*stdZ*stdW)
    ##    abX <- (Sa1+Sz1*SstdZ+Sw1*SstdW+Szw1*SstdZ*SstdW)*(Sa2+Sz2*SstdZ+Sw2*SstdW+Szw2*SstdZ*SstdW)*(Sa3+Sz3*SstdZ+Sw3*SstdW+Szw3*SstdZ*SstdW)*
    ## (Sa4+Sz4*SstdZ+Sw4*SstdW+Szw4*SstdZ*SstdW)

    ## Hi/Lo ##
    ## estX <- (a1+z1*stdZ-w1*stdW-zw1*stdZ*stdW)*(a2+z2*stdZ-w2*stdW-zw2*stdZ*stdW)*(a3+z3*stdZ-w3*stdW-zw3*stdZ*stdW)*(a4+z4*stdZ-w4*stdW-zw4*stdZ*stdW)
    ##    abX <- (Sa1+Sz1*SstdZ-Sw1*SstdW-Szw1*SstdZ*SstdW)*(Sa2+Sz2*SstdZ-Sw2*SstdW-Szw2*SstdZ*SstdW)*(Sa3+Sz3*SstdZ-Sw3*SstdW-Szw3*SstdZ*SstdW)*
    ## (Sa4+Sz4*SstdZ-Sw4*SstdW-Szw4*SstdZ*SstdW)

    ## Lo/Hi ##
    ## estX <- (a1-z1*stdZ+w1*stdW-zw1*stdZ*stdW)*(a2-z2*stdZ+w2*stdW-zw2*stdZ*stdW)*(a3-z3*stdZ+w3*stdW-zw3*stdZ*stdW)*(a4-z4*stdZ+w4*stdW-zw4*stdZ*stdW)
    ## abX <- (Sa1-Sz1*SstdZ+Sw1*SstdW-Szw1*SstdZ*SstdW)*(Sa2-Sz2*SstdZ+Sw2*SstdW-Szw2*SstdZ*SstdW)*(Sa3-Sz3*SstdZ+Sw3*SstdW-Szw3*SstdZ*SstdW)*
    ## (Sa4-Sz4*SstdZ+Sw4*SstdW-Szw4*SstdZ*SstdW)

    ## Lo/Lo ##
    ## estX <- (a1-z1*stdZ-w1*stdW+zw1*stdZ*stdW)*(a2-z2*stdZ-w2*stdW+zw2*stdZ*stdW)*(a3-z3*stdZ-w3*stdW+zw3*stdZ*stdW)*(a4-z4*stdZ-w4*stdW+zw4*stdZ*stdW)
    ## abX <- (Sa1-Sz1*SstdZ-Sw1*SstdW+Szw1*SstdZ*SstdW)*(Sa2-Sz2*SstdZ-Sw2*SstdW+Szw2*SstdZ*SstdW)*(Sa3-Sz3*SstdZ-Sw3*SstdW+Szw3*SstdZ*SstdW)*
    ## (Sa4-Sz4*SstdZ-Sw4*SstdW+Szw4*SstdZ*SstdW)


    ## Slope Difference Test 6 "HiZ/LoW - LoZ/HiW" ##
    estX <- (a1+z1*stdZ-w1*stdW-zw1*stdZ*stdW)*(a2+z2*stdZ-w2*stdW-zw2*stdZ*stdW)*(a3+z3*stdZ-w3*stdW-zw3*stdZ*stdW)*(a4+z4*stdZ-w4*stdW-zw4*stdZ*stdW) -
            (a1-z1*stdZ+w1*stdW-zw1*stdZ*stdW)*(a2-z2*stdZ+w2*stdW-zw2*stdZ*stdW)*(a3-z3*stdZ+w3*stdW-zw3*stdZ*stdW)*(a4-z4*stdZ+w4*stdW-zw4*stdZ*stdW)
     abX <- (Sa1+Sz1*SstdZ-Sw1*SstdW-Szw1*SstdZ*SstdW)*(Sa2+Sz2*SstdZ-Sw2*SstdW-Szw2*SstdZ*SstdW)*(Sa3+Sz3*SstdZ-Sw3*SstdW-Szw3*SstdZ*SstdW)*
            (Sa4+Sz4*SstdZ-Sw4*SstdW-Szw4*SstdZ*SstdW) -
            (Sa1-Sz1*SstdZ+Sw1*SstdW-Szw1*SstdZ*SstdW)*(Sa2-Sz2*SstdZ+Sw2*SstdW-Szw2*SstdZ*SstdW)*(Sa3-Sz3*SstdZ+Sw3*SstdW-Szw3*SstdZ*SstdW)*
            (Sa4-Sz4*SstdZ+Sw4*SstdW-Szw4*SstdZ*SstdW)
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


    ## == Create 3-way Interaction Figure == ##

    ## Slopes ##
    TWx1 <- PCI[4,4,4] # Hi-Z, Hi-W
    TWx2 <- PCI[4,4,2] # Hi-Z, Lo-W
    TWx3 <- PCI[2,4,4] # Lo-Z, Hi-W
    TWx4 <- PCI[2,4,2] # Lo-Z, Lo-W

#    TWx1 <- .0837
#    TWx2 <- -.2347
#    TWx3 <- -.0099
#    TWx4 <- .2497

# library(tidyr)

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

    ggplot(data = df_long, aes(x = x_var, y = y_value, color = line_id, linetype = line_id)) +
      xlim(-2.5, 2.5) +
      ylim(y.lower, y.upper) +
      geom_line(linewidth=1) +
      scale_linetype_manual(name = "My Legend Title",
                            values = c("line1" = "solid", "line2" = "twodash", "line3" = "solid", "line4" = "twodash"),
                            labels = c("Hi-Z, Hi-W", "Hi-Z, Lo-W", "Lo-Z, Hi-W", "Lo-Z, Lo-W")) +
      scale_color_manual(values = c("line1" = "black", "line2" = "black", "line3" = "grey", "line4"="grey")) +
      guides(color = "none") +
      labs(title = "Four Lines on One Graph",
           x = "X-axis Label",
           y = "Y-axis Label") +
      theme_minimal() # Optional: Apply a theme
      ggplot2::ggsave("3-Way Interaction Figure.png", width = 22.86, height = 16.51, units = "cm")

    cat("\n")
    cat("Figure is saved in '3-Way Interaction Figure.png'", rep("\n", 2))

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

  return(list(a1=a1, a2=a2, a3=a3, a4=a4, z1=z1, z2=z2, z3=z3, z4=z4, 
              w1=w1, w2=w2, w3=w3, w4=w4, sd_z=sd_z,
              Sa1=Sa1, Sa2=Sa2, Sa3=Sa3, Sa4=Sa4, Sz1=Sz1, Sz2=Sz2, Sz3=Sz3, Sz4=Sz4, 
              Sw1=Sw1, Sw2=Sw2, Sw3=Sw3, Sw4=Sw4, 
              b.no = b.no, NoModz = NoModz, NoModw = NoModw, Z=Z)) 
}


## Simulation of Designed Function ##
mccimm_modsem_fun <- function(object = est_lms, Z="NA", W="NA",
                   A1="NA", Z1="NA", W1="NA", ZW1="NA",
                   A2="NA", Z2="NA", W2="NA", ZW2="NA",
                   A3="NA", Z3="NA", W3="NA", ZW3="NA",
                   A4="NA", Z4="NA", W4="NA", ZW4="NA",
                   Sfunction="NULL", R=5) {


  ## -- Initialize variables -- ##
  a1 <- 0
  a2 <- 0
  a3 <- 0
  a4 <- 0
  z1 <- 0
  z2 <- 0
  z3 <- 0
  z4 <- 0
  w1 <- 0
  w2 <- 0
  w3 <- 0
  w4 <- 0
  zw1 <- 0
  zw2 <- 0
  zw3 <- 0
  zw4 <- 0
  Sa1 <- matrix(1, R*1e6)
  Sa2 <- matrix(1, R*1e6)
  Sa3 <- matrix(1, R*1e6)
  Sa4 <- matrix(1, R*1e6)
  Sz1 <- matrix(0, R*1e6)
  Sz2 <- matrix(0, R*1e6)
  Sz3 <- matrix(0, R*1e6)
  Sz4 <- matrix(0, R*1e6)
  Sw1 <- matrix(0, R*1e6)
  Sw2 <- matrix(0, R*1e6)
  Sw3 <- matrix(0, R*1e6)
  Sw4 <- matrix(0, R*1e6)
  Szw1 <- matrix(0, R*1e6)
  Szw2 <- matrix(0, R*1e6)
  Szw3 <- matrix(0, R*1e6)
  Szw4 <- matrix(0, R*1e6)
  ## ------------------------------ ##


  ## Extract defined parameters and vcov ##
  varZ <- "NA"
  varW <- "NA"
  if (Z != "NA") varZ <- paste0(Z, "~~", Z)
  if (W != "NA") varW <- paste0(W, "~~", W)
  dp <- c(A1, A2, A3, A4, Z1, Z2, Z3, Z4, W1, W2, W3, W4, ZW1, ZW2, ZW3, ZW4, varZ, varW)
  dp <- dp[dp != "NA"]

  temp <- modsem_coef(object)
  estcoeff <- temp[dp]
  Temp3 <- modsem_vcov(object)
  Tech3 <- Temp3[dp, dp]


  ## -- Monte Carlo Simulation of R*1e6 samples, default: R = 5 -- ##
  mcmc <- MASS::mvrnorm(n=R*1e6, mu=estcoeff, Sigma=Tech3, tol = 1e-6)

  b.no <- nrow(mcmc)
  R.no <- format(R*1e6, scientific = FALSE)

  # ===== Print number of bootstrap samples
  cat("\n", "   Number of requested simulated sample = ", R.no)
  cat("\n", "   Number of completed simulated sample = ", b.no, rep("\n",2))


  cat("Simulated Defined Function Values", rep("\n", 2))

  # ==================================================================== #

  # Define estimated parameters for calculating indirect effects
  if (any(names(estcoeff) %in% A1)) a1 <- estcoeff[A1]
  if (any(names(estcoeff) %in% A2)) a2 <- estcoeff[A2]
  if (any(names(estcoeff) %in% A3)) a3 <- estcoeff[A3]
  if (any(names(estcoeff) %in% A4)) a4 <- estcoeff[A4]
  if (any(names(estcoeff) %in% Z1)) z1 <- estcoeff[Z1]
  if (any(names(estcoeff) %in% Z2)) z2 <- estcoeff[Z2]
  if (any(names(estcoeff) %in% Z3)) z3 <- estcoeff[Z3]
  if (any(names(estcoeff) %in% Z4)) z4 <- estcoeff[Z4]
  if (any(names(estcoeff) %in% W1)) w1 <- estcoeff[W1]
  if (any(names(estcoeff) %in% W2)) w2 <- estcoeff[W2]
  if (any(names(estcoeff) %in% W3)) w3 <- estcoeff[W3]
  if (any(names(estcoeff) %in% W4)) w4 <- estcoeff[W4]
  if (any(names(estcoeff) %in% ZW1)) zw1 <- estcoeff[ZW1]
  if (any(names(estcoeff) %in% ZW2)) zw2 <- estcoeff[ZW2]
  if (any(names(estcoeff) %in% ZW3)) zw3 <- estcoeff[ZW3]
  if (any(names(estcoeff) %in% ZW4)) zw4 <- estcoeff[ZW4]

  # Calculate estimated parameter from Dfunction
  estM  <- eval(parse(text=Sfunction))

  # Capture simulated parameters for calculating function values
  if (any(names(estcoeff) %in% A1)) Sa1 <- mcmc[, A1]
  if (any(names(estcoeff) %in% A2)) Sa2 <- mcmc[, A2]
  if (any(names(estcoeff) %in% A3)) Sa3 <- mcmc[, A3]
  if (any(names(estcoeff) %in% A4)) Sa4 <- mcmc[, A4]
  if (any(names(estcoeff) %in% Z1)) Sz1 <- mcmc[, Z1]
  if (any(names(estcoeff) %in% Z2)) Sz2 <- mcmc[, Z2]
  if (any(names(estcoeff) %in% Z3)) Sz3 <- mcmc[, Z3]
  if (any(names(estcoeff) %in% Z4)) Sz4 <- mcmc[, Z4]
  if (any(names(estcoeff) %in% W1)) Sw1 <- mcmc[, W1]
  if (any(names(estcoeff) %in% W2)) Sw2 <- mcmc[, W2]
  if (any(names(estcoeff) %in% W3)) Sw3 <- mcmc[, W3]
  if (any(names(estcoeff) %in% W4)) Sw4 <- mcmc[, W4]
  if (any(names(estcoeff) %in% ZW1)) Szw1 <- mcmc[, ZW1]
  if (any(names(estcoeff) %in% ZW2)) Szw2 <- mcmc[, ZW2]
  if (any(names(estcoeff) %in% ZW3)) Szw3 <- mcmc[, ZW3]
  if (any(names(estcoeff) %in% ZW4)) Szw4 <- mcmc[, ZW4]

  # Calculate Simulated parameter from Dfunction
  abM <- eval(parse(text=Sfunction))


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
}  ## End function mccimm_modsem_fun ##




## -- FUNCTION JN_plot to plot Johnson-Neyman Figure -- ##
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
  jnp <- ggplot2::ggplot(df_plot, ggplot2::aes(x = z, y = slope)) +
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
       ggplot2::theme_minimal()
       if (length(approx_jn) > 0) {
         top_y <- suppressWarnings(max(df_plot$slope[is.finite(df_plot$slope)], na.rm = TRUE))
         if (!is.finite(top_y)) top_y <- y_range[2]
         for (zstar in approx_jn) {
           if (is.finite(zstar) && zstar >= min_z_abs && zstar <= max_z_abs) {
             jnp <- jnp + ggplot2::geom_vline(xintercept = zstar, linetype = "dashed", color = "red") +
                    ggplot2::annotate("text", x = zstar, y = top_y, label = paste("JN point (~):",
                              round(zstar, 2)), hjust = -0.1, vjust = 1, color = "black")
           }
         }
       }

  ## -- Print the Figure -- ##      
  print(jnp)
  return(jnp)

} # End plot Johnson-Neyman Figure 
