# examine 5 factor EFA for noHighCor_all_data
sequence <- seq(20)
for (i in sequence)
{
  i
  fac1 <- factanal(noHighCor_all_data, i, rotation = "promax")
  print(fac1, cutoff = .30)
  fl <- round(unclass(fac1$loadings), 3)
  fl
  file <- paste("output/noHighCor_all_data-",toString(i),"factor-efa-factor-loadings.csv", collapse = "")
  write.csv(fl, file)
}

#write.csv(fl, "output/noHighCor_all_data-2factor-efa-factor-loadings.csv")
#  fac1 <- factanal(noHighCor_all_data, 1, rotation = "promax")
#cbind("output/noHighCor_all_data-",toString(i),"factor-efa-factor-loadings.csv")



# examine 5 factor EFA for noHighCor_math
sequence <- seq(20)
for (i in sequence)
{
  i
  fac1 <- factanal(noHighCor_math, i, rotation = "promax")
  print(fac1, cutoff = .30)
  fl <- round(unclass(fac1$loadings), 3)
  fl
  file <- paste("output/noHighCor_math-",toString(i),"factor-efa-factor-loadings.csv", collapse = "")
  write.csv(fl, file)
}