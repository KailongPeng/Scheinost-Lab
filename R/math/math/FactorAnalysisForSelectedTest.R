
noHighCor_SelectedTest <- read.table("raw_data/noHighCor_SelectedTest_test_no_norm_no_nan_header.txt", header=TRUE, sep=",")

noHighCor_SelectedTest <- data.frame(sapply(noHighCor_SelectedTest, as.numeric))

write.csv(noHighCor_SelectedTest, file = "output/noHighCor_SelectedTest.csv", row.names = FALSE)

library(lavaan)
noHighCor_SelectedTest <- read.table("output/noHighCor_SelectedTest.csv", header=TRUE, sep=",")

# check for number of factors
psych::scree(noHighCor_SelectedTest)
psych::fa.parallel(noHighCor_SelectedTest)

# examine exploratory factor analysis (EFA) for noHighCor_SelectedTest
sequence <- seq(20)
for (i in sequence)
{
  i
  fac1 <- factanal(noHighCor_SelectedTest, i, rotation = "promax")
  print(fac1, cutoff = .30)
  fl <- round(unclass(fac1$loadings), 3)
  fl
  file <- paste("output/noHighCor_SelectedTest-",toString(i),"factor-efa-factor-loadings.csv", collapse = "")
  write.csv(fl, file)
}

{
# Confirmatory factor analysis for math data
models <- list()
fits <- list()
{# 1 factor model 
  models$m1 <- 
    '    factor1 =~ AWMA0x2DS_VerbalSTM_StS + AWMA0x2DS_VerbalWM_StS + AWMA0x2DS_VisuoSpatialSTM_StS + AWMA0x2DS_VisuoSpatialWM_StS + CMAT_BasicCalc_Comp_Quotient + CTOPP_EL_StS + CTOPP_BW_StS + KeyMath_Numeration_ScS + KeyMath_ProblemSolving_ScS + TOWRE_Total_StS + WASI_Vocab_T0x2DScore + WASI_BD_T0x2DScore + WASI_Sim_T0x2DScore + WASI_MR_T0x2DScore + WJ0x2DIII_WordID_StS + WJ0x2DIII_WA_StS + WJ0x2DIII_PassComp_StS + WJ0x2DIII_MathFluency_StS + WJ0x2DIII_SpatialRelations_StS'
  
  
  
  # 2 factor model 
  models$m2 <- 
    'factor1 =~ AWMA0x2DS_VisuoSpatialSTM_StS + AWMA0x2DS_VisuoSpatialWM_StS + CMAT_BasicCalc_Comp_Quotient + KeyMath_Numeration_ScS + KeyMath_ProblemSolving_ScS + TOMA0x2D2_Attitudes_StS + WASI_Vocab_T0x2DScore + WASI_BD_T0x2DScore + WASI_Sim_T0x2DScore + WASI_MR_T0x2DScore + WJ0x2DIII_SpatialRelations_StS
  
  factor2 =~ AWMA0x2DS_VerbalSTM_StS + CTOPP_EL_StS + CTOPP_BW_StS + CTOPP_RapidNaming_Comp + TOWRE_Total_StS + WJ0x2DIII_WordID_StS + WJ0x2DIII_WA_StS + WJ0x2DIII_PassComp_StS + WJ0x2DIII_MathFluency_StS'
  
  
  
  # 3 factor model 
  models$m3 <- 
    'factor1 =~ AWMA0x2DS_VisuoSpatialSTM_StS + AWMA0x2DS_VisuoSpatialWM_StS + CMAT_BasicCalc_Comp_Quotient + KeyMath_Numeration_ScS + KeyMath_ProblemSolving_ScS + TOMA0x2D2_Attitudes_StS + WASI_Vocab_T0x2DScore + WASI_BD_T0x2DScore + WASI_Sim_T0x2DScore + WASI_MR_T0x2DScore + WJ0x2DIII_MathFluency_StS + WJ0x2DIII_SpatialRelations_StS
  
  factor2 =~ AWMA0x2DS_VerbalSTM_StS + AWMA0x2DS_VerbalWM_StS + CTOPP_EL_StS + CTOPP_BW_StS + KeyMath_ProblemSolving_ScS + TOWRE_Total_StS + WASI_Vocab_T0x2DScore + WASI_Sim_T0x2DScore + WJ0x2DIII_WordID_StS + WJ0x2DIII_WA_StS + WJ0x2DIII_PassComp_StS
  
  factor3 =~ CTOPP_RapidNaming_Comp + TOWRE_Total_StS + WJ0x2DIII_MathFluency_StS'
  
  
  
  # 4 factor model 
  models$m4 <- 
    'factor1 =~ AWMA0x2DS_VerbalSTM_StS + AWMA0x2DS_VerbalWM_StS + CTOPP_EL_StS + CTOPP_BW_StS + KeyMath_Numeration_ScS + KeyMath_ProblemSolving_ScS + TOWRE_Total_StS + WASI_Vocab_T0x2DScore + WASI_Sim_T0x2DScore + WJ0x2DIII_WordID_StS + WJ0x2DIII_WA_StS + WJ0x2DIII_PassComp_StS
  
  factor2 =~ AWMA0x2DS_VisuoSpatialSTM_StS + AWMA0x2DS_VisuoSpatialWM_StS + CMAT_BasicCalc_Comp_Quotient + KeyMath_Numeration_ScS + KeyMath_ProblemSolving_ScS + TOMA0x2D2_Attitudes_StS + WJ0x2DIII_MathFluency_StS
  
  factor3 =~ WASI_BD_T0x2DScore + WASI_MR_T0x2DScore + WJ0x2DIII_SpatialRelations_StS
  
  factor4 =~ CTOPP_RapidNaming_Comp + TOWRE_Total_StS + WJ0x2DIII_WordID_StS + WJ0x2DIII_MathFluency_StS'
  
  
  
  # 5 factor model 
  models$m5 <- 
    'factor1 =~ AWMA0x2DS_VerbalSTM_StS + TOWRE_Total_StS + WJ0x2DIII_WordID_StS + WJ0x2DIII_WA_StS + WJ0x2DIII_PassComp_StS
  
  factor2 =~ AWMA0x2DS_VisuoSpatialWM_StS + CTOPP_BW_StS + WASI_BD_T0x2DScore + WASI_MR_T0x2DScore + WJ0x2DIII_SpatialRelations_StS
  
  factor3 =~ AWMA0x2DS_VisuoSpatialSTM_StS + AWMA0x2DS_VisuoSpatialWM_StS + CMAT_BasicCalc_Comp_Quotient + KeyMath_Numeration_ScS + KeyMath_ProblemSolving_ScS + TOMA0x2D2_Attitudes_StS + WJ0x2DIII_MathFluency_StS
  
  factor4 =~ WASI_Vocab_T0x2DScore + WASI_Sim_T0x2DScore
  
  factor5 =~ CTOPP_RapidNaming_Comp + TOWRE_Total_StS + WJ0x2DIII_MathFluency_StS'
  
  
  
  # 6 factor model 
  models$m6 <- 
    'factor1 =~ AWMA0x2DS_VerbalSTM_StS + TOWRE_Total_StS + WJ0x2DIII_WordID_StS + WJ0x2DIII_WA_StS + WJ0x2DIII_PassComp_StS
  
  factor2 =~ WASI_Vocab_T0x2DScore + WASI_Sim_T0x2DScore
  
  factor3 =~ WASI_BD_T0x2DScore + WASI_MR_T0x2DScore + WJ0x2DIII_SpatialRelations_StS
  
  factor4 =~ AWMA0x2DS_VerbalWM_StS + AWMA0x2DS_VisuoSpatialSTM_StS + AWMA0x2DS_VisuoSpatialWM_StS + KeyMath_Numeration_ScS + KeyMath_ProblemSolving_ScS
  
  factor5 =~ CMAT_BasicCalc_Comp_Quotient + KeyMath_Numeration_ScS + TOMA0x2D2_Attitudes_StS + WJ0x2DIII_MathFluency_StS
  
  factor6 =~ CTOPP_RapidNaming_Comp + TOWRE_Total_StS + WJ0x2DIII_MathFluency_StS'
  
  
  
  # 7 factor model 
  models$m7 <- 
    'factor1 =~ TOWRE_Total_StS + WJ0x2DIII_WordID_StS + WJ0x2DIII_WA_StS + WJ0x2DIII_PassComp_StS
  
  factor2 =~ AWMA0x2DS_VerbalWM_StS + AWMA0x2DS_VisuoSpatialSTM_StS + AWMA0x2DS_VisuoSpatialWM_StS + KeyMath_Numeration_ScS + KeyMath_ProblemSolving_ScS
  
  factor3 =~ WASI_BD_T0x2DScore + WASI_MR_T0x2DScore + WJ0x2DIII_SpatialRelations_StS
  
  factor4 =~ WASI_Vocab_T0x2DScore + WASI_Sim_T0x2DScore
  
  factor5 =~ CMAT_BasicCalc_Comp_Quotient + KeyMath_Numeration_ScS + TOMA0x2D2_Attitudes_StS + WJ0x2DIII_MathFluency_StS
  
  factor6 =~ CTOPP_RapidNaming_Comp + TOWRE_Total_StS + WJ0x2DIII_MathFluency_StS
  
  factor7 =~ AWMA0x2DS_VerbalSTM_StS'
  
  
  
  fits$m1 <- lavaan::cfa(models$m1, data = noHighCor_SelectedTest)
  fits$m2 <- lavaan::cfa(models$m2, data = noHighCor_SelectedTest)
  fits$m3 <- lavaan::cfa(models$m3, data = noHighCor_SelectedTest)
  fits$m4 <- lavaan::cfa(models$m4, data = noHighCor_SelectedTest)
  fits$m5 <- lavaan::cfa(models$m5, data = noHighCor_SelectedTest)
  fits$m6 <- lavaan::cfa(models$m6, data = noHighCor_SelectedTest)
  fits$m7 <- lavaan::cfa(models$m7, data = noHighCor_SelectedTest)
}
summary(fits$m1, fit.measures = TRUE) 
standardizedSolution(fits$m1)
v$fitindicies <- c("npar",  "chisq", "df", "pvalue", "cfi", "rmsea", 
                   "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
round(sapply(fits, function(X) fitmeasures(X)[v$fitindicies]), 3)


summary(fits$m2, fit.measures = TRUE) 
standardizedSolution(fits$m2)
v$fitindicies <- c("npar",  "chisq", "df", "pvalue", "cfi", "rmsea", 
                   "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
round(sapply(fits, function(X) fitmeasures(X)[v$fitindicies]), 3)


summary(fits$m3, fit.measures = TRUE) 
standardizedSolution(fits$m3)
v$fitindicies <- c("npar",  "chisq", "df", "pvalue", "cfi", "rmsea", 
                   "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
round(sapply(fits, function(X) fitmeasures(X)[v$fitindicies]), 3)
a <- coef(fits$m3)

str_eval=function(x) {return(eval(parse(text=x)))}
sequence <- seq(7)
for (i in sequence)
{
  i
  #FactorScores <- predict(fits$m3)
  #FactorScores <- predict(eval(paste("fits$m",sep = "")))
  str <- paste("FactorScores <- predict(fits$m",i,")",sep = "")
  FactorScores = str_eval(str)
  file <- paste("output/noHighCor_SelectedTest_",toString(i),"FactorScores.csv",sep = "")
  write.csv(FactorScores, file)
}

}
