library(psych)
library(GPArotation)
## Not run:
test.data <- Harman74.cor$cov
# if(!require(GPArotation)) {message("Omega requires GPA rotation" )} else {
my.omega <- omega(test.data)
print(my.omega,digits=2)
#}

####################################################################

#simulate 9 variables with a hierarchical structure
v9 <- sim.hierarchical()
#with correlations of
round(v9,2)
#find omega
v9.omega <- omega(v9,digits=2)
v9.omega

####################################################################

#create 8 items with a two factor solution, showing the use of the flip option
sim2 <- item.sim(8)
omega(sim2) #an example of misidentification-- remember to look at the loadings matrices. 
omega(sim2,2) #this shows that in fact there is no general factor
omega(sim2,2,option="first") #but, if we define one of the two group factors
                             #as a general factor, we get a falsely high omega
#apply omega to analyze 6 mental ability tests
data(ability.cov)   #has a covariance matrix
omega(ability.cov$cov)
#om <- omega(Thurstone)
#round(om$omega.group,2)
#round(om$omega.group[2]/om$omega.group[1],2) #fraction of reliable that is general variance 
# round(om$omega.group[3]/om$omega.group[1],2) #fraction of reliable that is group variance
#To find factor score estimates for the hierarchical model it is necessary to
#do two extra steps.
#Consider the case of the raw data in an object data.  (An example from simulation)
# set.seed(42)
# gload <- matrix(c(.9,.8,.7),nrow=3)
# fload <- matrix(c(.8,.7,.6,rep(0,9),.7,.6,.5,rep(0,9),.7,.6,.4),   ncol=3)
# data <- sim.hierarchical(gload=gload,fload=fload, n=100000, raw=TRUE)
#
# f3 <- fa(data$observed,3,scores="tenBerge", oblique.scores=TRUE)
# f1 <- fa(f3$scores)
# om <- omega(data$observed,sl=FALSE) #draw the hierarchical figure
# The scores from om are based upon the Schmid-Leiman factors and although the g factor
# is identical, the group factors are not.
# This is seen in the following correlation matrix
# hier.scores <- cbind(om$scores,f1$scores,f3$scores)
# lowerCor(hier.scores)
#
#jensen <- sim.hierarchical() #create a hierarchical structure
#om.jen <- omegaSem(jensen,lavaan=TRUE) #do the exploratory omega with confirmatory as well


#lav.mod <- om.jen$omegaSem$model$lavaan #get the lavaan code or create it yourself
# lav.mod <- 'g =~ +V1+V2+V3+V4+V5+V6+V7+V8+V9
# F1=~ +V1+V2+V3
# F2=~ +V4+V5+V6
# F3=~ +V7+V8+V9'
#lav.jen <- cfa(lav.mod,sample.cov=jensen,sample.nobs=500,orthogonal=TRUE,std.lv=TRUE) 
# omegaFromSem(lav.jen,jensen)
#try a one factor solution -- this is not recommended, but sometimes done
#it will just give omega_total
# lav.mod.1 <- 'g =~ +V1+V2+V3+V4+V5+V6+V7+V8+V9 '
#lav.jen.1<- cfa(lav.mod.1,sample.cov=jensen,sample.nobs=500,orthogonal=TRUE,std.lv=TRUE) 
# omegaFromSem(lav.jen.1,jensen)
## End(Not run)