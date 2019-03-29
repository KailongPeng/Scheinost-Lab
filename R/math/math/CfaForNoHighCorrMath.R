# Confirmatory factor analysis for all data
```{r}
models <- list()
fits <- list()

# 1 factor model 
models$m1 <- 
  '    factor1 =~ CMAT_Addition_Raw + CMAT_Subtraction_Raw + CMAT_Multiplication_Raw + KeyMath_Numeration_Raw + KeyMath_ProblemSolving_Raw + TOMA0x2D2_Attitudes_Raw'



# 2 factor model 
models$m2 <- 
  'factor1 =~ CMAT_Addition_Raw + CMAT_Subtraction_Raw + CMAT_Multiplication_Raw + TOMA0x2D2_Attitudes_Raw

factor2 =~ KeyMath_Numeration_Raw + KeyMath_ProblemSolving_Raw'



# 3 factor model 
models$m3 <- 
  'factor1 =~ KeyMath_Numeration_Raw + KeyMath_ProblemSolving_Raw

factor2 =~ CMAT_Addition_Raw + CMAT_Subtraction_Raw + CMAT_Multiplication_Raw

factor3 =~ CMAT_Subtraction_Raw + KeyMath_Numeration_Raw + TOMA0x2D2_Attitudes_Raw'


#lavaan WARNING: some observed variances are (at least) a factor 1000 times larger than others; use varTable(fit) to investiga
fits$m1 <- lavaan::cfa(models$m1, data = noHighCor_math)
fits$m2 <- lavaan::cfa(models$m2, data = noHighCor_math)
fits$m3 <- lavaan::cfa(models$m3, data = noHighCor_math)

summary(fits$m2, fit.measures = TRUE)
standardizedSolution(fits$m2)

v$fitindicies <- c("npar",  "chisq", "df", "pvalue", "cfi", "rmsea", 
                   "rmsea.ci.lower", "rmsea.ci.upper", "srmr")

round(sapply(fits, function(X) fitmeasures(X)[v$fitindicies]), 3)
