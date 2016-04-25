#' Calculates creatinine clearance or outputs Fortran code for in-model calculations. Patient covariates are supplied as vectors of equal length. 
#' If no patient values are supplied, the function outputs the Fortran code to be used in a \emp{Pmetrics} model file.
#'
#' \code{CrCl} takes vectors with patient parameters and calculates creatinine clearance using specified formula.
#'
#' @title Calculate creatinine clearance
#' @param formula Required. Identifies the formula to be used for calculation. Available options are \code{Cocroft-Gault}, \code{MDRD}, \code{MDRD-BUN}, \code{CKD-EPI}, \code{Schwartz}, \code{Mayo}, \code{Jelliffe73}.
#' @param creat Required. A vector containing serum creatinine levels either in umol/l or mg/dL as specified by the \code{SI} parameter.
#' @param age Required. A vector containing age values in years.
#' @param male A vector specifying patients sex, 1 for male, 0 for female.
#' @param BW A vector containing body weight values in kilograms.
#' @param black A vector of 1's and 0's specifying if patient is black, used for the MDRD and CKD-EPI formulas.
#' @param BUN A vector of BUN values for the MDRD-BUN formula.
#' @param albumin A vector of albumin levels for the MDRD-BUN formula.
#' @param preterm A vector of 1's and 0's specifying if a child was born preterm. Used for the Schwartz formula.
#' @param height A vector of heights in \emph{cm}. Used for the Schwartz formula.
#' @param SI A logical operator. If \code{FALSE}, creatinine and BUN is in mg/dL, albumin in g/dL. Default is \code{TRUE}, i.e. umol/l, mmol/l, and g/l, respectively.
#' @return The output is a vector of creatinine clearance values. Please be aware that unlike other formulas, Cockroft-Gault returns ml/min!
#' @examples 
#' # this will return a vector of creatinine clearance values
#' eGFR <- CrCl(formula="Cockroft-Gault", creat=c(78,40,50), BW=c(70,50,30), age=c(50,78,30), male=c(1,0,0))
#' 
#' #this will print the Fortran code to be used in a model file
#' CrCl("Cockroft-Gault")
#' 
#' @author Jan Strojil
#' 
CrCl <- function(formula, creat, BW, age, male, black, BUN, albumin, preterm, height, SI=T){
  # define subfunctions for individual calculations, all in US units
  CG <- function(creat, age, BW, male){
    if (min(length(creat),length(age), length(male),length(BW))!=max(length(creat),length(age), length(male),length(BW))){
      stop("Vectors are not of equal length.", call. = F)
    }
    coef <- rep(1, length(creat))
    coef[male == 0] <- 0.85
    CrCl_CG <-  ((140-age)*BW*coef)/(72*creat)
  }
  CDK_EPI <- function(creat,age,male,black){
    if (min(length(creat),length(age), length(male),length(black))!=max(length(creat),length(age), length(male),length(black))){
      stop("Vectors are not of equal length.", call. = F)
    }
    if (SI) creat <- creat/88.4
    nsub <- length(creat)
    alpha <- rep(-0.329, nsub)
    kappa <- rep(0.7, nsub)
    coefF <- rep(1.018, nsub)
    coefB <- rep(1.159, nsub)
    alpha[male == 1]  <- -0.411
    kappa[male == 1]  <-  0.9
    coefF[male == 1]  <-  1
    coefB[black == 0]  <-  1
    maxi <-  creat/kappa
    mini <-  creat/kappa
    maxi[maxi < 1] <- 1
    mini[mini > 1] <- 1
    CrCl_CKD <-  141 * mini**alpha * maxi**(-1.209) * 0.993**age * coefF * coefB
  }
  MDRD <- function(creat, age, male, black){
    if (min(length(creat),length(age), length(male),length(black))!=max(length(creat),length(age), length(male),length(black))){
      stop("Vectors are not of equal length.", call. = F)
    }
    nsub <- length(creat)
    coefF <- rep(1, nsub)
    coefB <- rep(1, nsub)
    coefF[male == 0] <- 0.742
    coefB[black == 1] <- 1.210
    CrCl_MDRD <- 186 * creat**(-1.154) * age**(-0.203) * coefF * coefB
  }
  MDRD_BUN <- function(creat, age, male, black, BUN, albumin){
    if (min(length(creat),length(age), length(male),length(black),length(BUN),length(albumin))!=max(length(creat),length(age), length(male),length(black),length(BUN),length(albumin))){
      stop("Vectors are not of equal length.", call. = F)
    }
    nsub <- length(creat)
    coefF <- rep(1, nsub)
    coefB <- rep(1, nsub)
    coefF[male == 0] <-  0.742
    coefB[black == 1] <-  1.210
    CrCl_MDRD_BUN <-  170 * creat**(-1.154) * age**(-0.203) * coefF * coefB * BUN**(-0.17) * albumin**(0.318)
  }
  Mayo <- function(creat, age, male){
    if (min(length(creat),length(age), length(male))!=max(length(creat),length(age), length(male))){
      stop("Vectors are not of equal length.", call. = F)
    }
    coefF <- rep(1, length(creat))
    creat[creat < 0.8]  <- 0.8
    coefF[male == 0]  <-  0.205
    CrCl_MCQE <- 2.718**(1.911 + 5.249/creat - 2.114/creat**2 - 0.00686 * age - coefF)
  }
  Schwartz <- function(creat, age, height, preterm){
    if (min(length(creat),length(age), length(height), length(preterm))!=max(length(creat),length(age), length(height), length(preterm))){
      stop("Vectors are not of equal length.", call. = F)
    }
    if (max(age)>12)
    k <- rep(0.55, length(creat))
    k[age < 1 & preterm==1] <- 0.33
    k[age < 1 & preterm==0] <- 0.45
    CrCl_Schwartz = k * height / creat
  }
  Jelliffe <- function(age, creat, male){
    if (min(length(creat),length(age), length(male))!=max(length(creat),length(age), length(male))){
      stop("Vectors are not of equal length.", call. = F)
    }
    CrCl_Jelliffe = (98-16*((age-20)/20))/creat
    CrCl_Jelliffe[male == 0] <- 0.9 * CrCl_Jelliffe
  }
  if (missing(formula)){
    cat("No formula specified. Available options are:\n1) Cockroft-Gault\n2) MDRD\n3) MDRD-BUN\n4) CKD-EPI\n5) Schwartz\n6) Mayo\n7) Jelliffe73")
    stop("No formula specified, don't know what to do.", call. = F)
  }
  if (missing(creat) & !missing(formula)){ # print out Fortran code
    tCKD <- paste(
      "\nalpha=-0.329",
      "\nkappa=0.7",
      "\ncoefF=1.018",
      "\ncoefB=1.159",
      "\n&IF(MALE == 1) alpha = -0.411",
      "\n&IF(MALE == 1) kappa = 0.9", 
      "\n&IF(MALE == 1) coefF = 1",
      "\n&IF(BLACK == 0) coefB = 1",
      "\nmaxi = CREAT/kappa",
      "\nmini = CREAT/kappa",
      "\n&IF(maxi .LT. 1) maxi=1",
      "\n&IF(mini .GT. 1) mini=1",
      "\nCrCl_CKD =141 * mini**alpha * maxi**(-1.209) * 0.993**AGE * coefF * coefB", sep = "")
    
    formulas <- list("Cockroft-Gault"=c(paste("\nc ### Cocroft-Gault - US units (mg/dL, kg) ### returns: ml/min",
                                              "\ncoef = 1",
                                              "\n&IF(MALE == 0) coef = 0.85",
                                              "\nCrCl_CG = ((140-AGE) * BW * coef)/(72*CREAT)", sep = ""),
                                        paste("\nc ### Cocroft-Gault - SI units (umol/l, kg) ### returns: ml/min",
                                              "\ncoef = 1.23",
                                              "\n&IF(MALE == 0) coef = 1.04",
                                              "\nCrCl_CG = ((140-AGE) * BW * coef)/CREAT", sep = "")),
                     "CKD-EPI"=c(paste("\nc ### CKD-EPI - US units (mg/dL) ### returns: mL/min/1.73 m²",
                                       tCKD, sep = ""),
                                 paste("\nc ### CKD-EPI - SI units (umol/l) ### returns: mL/min/1.73 m²",
                                       "\nCREAT = CREAT / 88.4", tCKD, sep = "")),
                     "MDRD"=c(paste("\nc ### MDRD - US units (mg/dL) ### returns: mL/min/1.73 m²",
                                    "\ncoefF=1",
                                    "\ncoefB=1",
                                    "\n&IF(MALE == 0) coefF = 0.742",
                                    "\n&IF(BLACK == 1) coefB = 1.210",
                                    "\nCrCl_MDRD = 186 * CREAT**(-1.154) * AGE**(-0.203) * coefF * coefB", sep = ""),
                              paste("\nc ### MDRD - SI units (umol/L) ### returns: mL/min/1.73 m²",
                                    "\ncoefF=1",
                                    "\ncoefB=1",
                                    "\n&IF(MALE == 0) coefF = 0.742",
                                    "\n&IF(BLACK == 1) coefB = 1.210",
                                    "\nCrCl_MDRD = 32788 * CREAT**(-1.154) * AGE**(-0.203) * coefF * coefB", sep = "")),
                     "MDRD-BUN"=c(paste("\nc ### MDRD-BUN - US units (CREAT and BUN in mg/dL, ALB in g/dL) ### returns: mL/min/1.73 m²",
                                        "\ncoefF=1",
                                        "\ncoefB=1",
                                        "\n&IF(MALE == 0) coefF = 0.762",
                                        "\n&IF(BLACK == 1) coefB = 1.180", 
                                        "\nCrCl_MDRD-BUN = 170 * CREAT**(-1.154) * AGE**(-0.203) * coefF * coefB * BUN**(-0.17) * ALB**(0.318)", sep = ""),
                                  paste("\nc ### MDRD-BUN - SI units (CREAT in umol/L, BUN in mmol/L, ALB in g/L) ### returns: mL/min/1.73 m²",
                                        "\ncoefF=1",
                                        "\ncoefB=1",
                                        "\n&IF(MALE == 0) coefF = 0.742",
                                        "\n&IF(BLACK == 1) coefB = 1.210",
                                        "\nCrCl_MDRD-BUN = 170 * CREAT**(-1.154) * AGE**(-0.203) * coefF * coefB * (BUN/0.3571)**(-0.17) * (ALB/10)**(0.318)", sep = "")),
                     "Mayo"=c(paste("\nc ### Mayo Quadratic - US units (CREAT in mg/dL) ### returns: mL/min/1.73 m²",
                                    "\ncoefF=1",
                                    "\n&IF(CREAT < 0.8) CREAT=0.8",
                                    "\n&IF(MALE == 0) coefF = 0.205",
                                    "\nCrCl_MCQ = 2.718**(1.911 + 5.249/CREAT - 2.114/CREAT**2 - 0.00686 * AGE - coefF)", sep = ""),
                              paste("\nc ### Mayo Quadratic - SI units (CREAT in umol/L) ### returns: mL/min/1.73 m²",
                                    "\ncoefF=1",
                                    "\nCREAT = CREAT / 88.4",
                                    "\n&IF(CREAT < 0.8) CREAT=0.8",
                                    "\n&IF(MALE == 0) coefF = 0.205",
                                    "\nCrCl_MCQ = 2.718**(1.911 + 5.249/CREAT - 2.114/CREAT**2 - 0.00686 * AGE - coefF)", sep = "")),
                     "Schwartz"=c(paste("\nc ### Schwartz - US units (CREAT in mg/dL, height in cm) ### returns: mL/min/1.73 m²",
                                        "\nk=0.55",
                                        "\n&IF(AGE.LT.1 .AND. PRETERM==1) k=0.33",
                                        "\n&IF(AGE.LT.1 .AND. PRETERM==0) k=0.45",
                                        "\nCrCl_Schwartz = k * HEIGHT / CREAT", sep = ""),
                                  paste("\nc ### Schwartz - SI units (CREAT in umol/L, height in cm) ### returns: mL/min/1.73 m²",
                                        "\nCREAT = CREAT/88.4",
                                        "\nk=0.55",
                                        "\n&IF(AGE.LT.1 .AND. PRETERM==1) k=0.33",
                                        "\n&IF(AGE.LT.1 .AND. PRETERM==0) k=0.45",
                                        "\nCrCl_Schwartz = k * HEIGHT / CREAT", sep = "")),
                     "Jelliffe"=c(paste("\nc #### Jellife - US units (creat in mg/dL) ### returns: mL/min/1.73 m²",
                                        "\nCrCl_Jelliffe = (98-16*((AGE-20)/20))/CREAT",
                                        "\n&IF(MALE == 0) CrCl_Jelliffe = 0.9 * CrCl_Jelliffe", sep = ""),
                                  paste("\nc #### Jellife - SI units (creat in umol/l) ### returns: mL/min/1.73 m²",
                                        "\nCrCl_Jelliffe = (98-16*((AGE-20)/20))/(CREAT/88.4)",
                                        "\n&IF(MALE == 0) CrCl_Jelliffe = 0.9 * CrCl_Jelliffe", sep = ""))
    )
    
    cat("\nCopy the code below into the #sec block of your model file, renaming the covariates as needed:\n")
    cat(unlist(formulas[formula])[as.integer(SI)+1])
    
  } # end of PRINT Fortran code
  else {
    
    # checks
    if (missing(creat)){
      stop("Missing serum creatinine values.", call. = F)
    }
    if (missing(age)){
      stop("Missing age values.", call. = F)
    }
    if (missing(male) & formula!="Schwartz"){
      stop("Missing sex.", call. = F)
    }
    if (SI & mean(creat)<10){
      warning("Creatinine levels are very low, are you sure they are in SI units?", call. = F)
    }
    if (!SI & mean(creat)>10){
      warning("Creatinine levels are very high, are you sure they are in US units?", call. = F)
    }
    
    # if SI, convert to US (shameful, but most formulas are simpler in US units)
    if (SI){
      creat <- creat / 88.4
      if (!missing(albumin)) albumin <- albumin*10
      if (!missing(BUN)) BUN <- BUN*0.3571
    }
    
    # call function with corresponding formula
    switch(formula,
           "Cockroft-Gault"={rval <- CG(creat, age, BW, male)},
           "CKD-EPI"= {rval <- CDK_EPI(creat,age,male,black)},
           "MDRD"={rval <- MDRD(creat, age, male, black)},
           "MDRD-BUN"={rval <- MDRD_BUN(creat, age, male, black, BUN, albumin)},
           "Mayo"={rval <- Mayo(creat, age, male)},
           "Schwartz"={rval <- Schwartz(creat, age, height, preterm)},
           "Jelliffe73"={rval <- Jelliffe(age, creat, male)}
    )
    return(rval)
  }
}