##### 

### A Taste of Fishery Science
## Code to support Lecture 2
#### H. Bowlby; Feb 2023

#setwd('C:/mydocs/_courses/DALmodule')


###--------------------------------------------------------

## population growth with and without fishing

### combining survival, reproductive output and Fishing

###--------------------------------------------------------

lotka.r<-function(age.maturity,litter.size,sex.ratio,gestation,max.age,u,sel) 
{
  
  temp<-c(0,c(rep(0,age.maturity+gestation-1),litter.size*sex.ratio,rep(c(litter.size*sex.ratio),max.age))) 
  mx<-temp[1:(max.age+1)]

  # partition survival into natural and fishing mortality  
  M<-2.56*max.age^-0.873  ### just trust me - you can approximate M from longevity data for sharks.
  f<--log(1-u)
  f.vec<-c(rep(0,sel),rep(f,max.age-sel))   ## this accounts for age at selectivity to the fishery   
  si<-exp(-(M+f.vec)) ## combine the two sources of mortality
  
  ## calculate survival at age
  x<-0:(max.age)
  lx<-0
  lx[1]<-1
  for(i in 2:(length(si)+1))
  {lx[i]<-lx[i-1]*si[i-1]}

  lxmx<-lx*mx   ### calculate reproductive output by age

  ### use minimization to estimate r from the parameters; minimization done based on the sum of squares.
  minimize<-function(start.r)
  {
    rx<-exp(-start.r*x)
    lotka<-sum(lxmx*rx)
    sumsq <- sum((lotka-1)^2) 
    return(sumsq)
  }
  junk<-nlminb(start = 0.1,objective = minimize, lower=-1,upper=1)
  return(junk)
}

### Running the function - all you need to do is specify parameters

lotka.r(10,8,0.5,2,40,0,1) # e.g.age.maturity = 10,litter.size= 8,sex.ratio=0.5,gestation=2,max.age=40,u=0,sel=1 NOTE: age 1 animals are selected to the fishery (i.e. the fishery can catch age 1 to max age)
## The value you are interested in is $par
## this represents the expected rate of change in population size each year.
## remember Nt+1 = r*Nt

# New aging data has come out for the Dusky Scallop Shark that was validated with bomb radiocarbon.
# this was long-awaited so that population productivity could be estimated for this data-deficient species prior to any fishery
# reproductive output (litter size and gestation) were already known from samples of pregnant females

### Here are the parameters: 
age.maturity= 30  ## females mature at 30 years old
litter.size= 8  ## females have litters of 8 pups
sex.ratio=0.5  ## half of the pups are female
gestation=2   ## it takes 2 years from conception to birth
max.age=40  ## Dusky scallop shark have a maximum age of 40 
## there is no fishery right now
u=0  ## the exploitation rate is zero
sel=1 ## NOTE: selectivity only matters if u > 0. 
      ##A value of 1 means age 1 animals are selected to the fishery (i.e. the fishery can catch age 1 to max age)


## QUESTION 1: Would it be meaningful to do a stock assessment on this species given our understanding of life history?
###             Why or why not? [3-4 sentences]



##----------------------------------------

## exploring S-R functions

##----------------------------------------
## Fit a beverton-holt function - 2 parameters
## parameterized in terms of alpha and alphaK
srfbh2 <- function(a,aK,S) {
  # Stock recruitment function: Beverton-Holt
  #  return(a*S/(S/K + 1))  
  return((1/S)*((1/a) + (S/aK)))
}

srfbhpv2 <- function(pv,S) {
  # Stock recruitment function (parameter vector version): Beverton-Holt
  return(srfbh2(pv[1],pv[2],S))
}
## initialize values
bhinit <- function(s,r) {
  # Calculate initial parameters for fitting Beverton-Holt model - uses the median values of the S-R data to calculate an initial alpha and R0.
  middle <- length(r)/2
  orderRec <- order(r)
  sortedRec <- r[orderRec]
  recsortedStock <- s[orderRec]
  K <- recsortedStock[middle]
  medRec <- sortedRec[middle]
  a <- 2*medRec/K
  param <- c(a,K)
  return(param)
}
## use maximum likelihood to estimate parameters
lnlnegloglik <- function(fitted,observed) {
  # Lognormal [lnl] negative [neg] log likelihood [lik]
  # Arguments:
  #   fitted   : fitted median recruitment (assumed to be free of missing values)
  #   observed : observed recruitment (assumed to be free of missing values)
  len <- length(fitted)
  sumlogr <- sum(log(observed))
  sumsqrlogratio <- sum((log(observed/fitted))^2)  
  result <- 0.5*len*log( (2*pi/len)*sumsqrlogratio) + sumlogr + len/2
  return(result)  
}
## fit the function
ml.bhlnl2 <- function(s,r,ip, returnlog=F) {
  #
  # alternate formulation for bh model (alpha and alphaK instead
  # of alpha and K)
  
  choice <- !is.na(r) & !is.na(s)
  r <- r[choice]
  s <- s[choice]
  # frame argument is unneccessary in R
  assign(".Stock",s)
  assign(".Recruit",r)
  assign(".Evaluations",0)
  
  if (missing(ip)) {
    ip <- bhinit(s,r)
  }
  logalpha <- log(ip[1])
  logbeta  <- log(ip[1]*ip[2])
  
  logip <- c(logalpha,logbeta)
  
  bh2.nll <- function(x) {
    assign(".Evaluations",.Evaluations+1) # frame argument is unnecesary in R
    
    invfitted <- srfbhpv2(exp(x),.Stock)
    fitted <- 1/(invfitted)
    return(lnlnegloglik(fitted,.Recruit))  
  }
  
  # modified to use nlminb function in R instead of nlmin in S (starting values supplied first, then function to minimize)
  nlmin.out <- nlminb(logip,bh2.nll) 
  # see R documentation for description of the control parameters.  Have left them as defaults until can chat with Jamie.
  #  x <- exp(nlmin.out$x) #x was the value at which the optimizer converged in the S nlmin function
  x<-exp(nlmin.out$par) # par is the best set of parameters
  invfitted <- srfbhpv2(x,.Stock)
  fitted <- 1/(invfitted)
  nll <- lnlnegloglik(fitted,.Recruit)
  sigma.squared <- sum( (log(.Recruit/fitted))^2 )/length(.Recruit)
  
  if(returnlog==T){
    return(log(x))
  } else {
    #careful with the transformation (not right?) 
    ## THIS WAS IN THE ORIGINAL NOTES - it has to deal with transformation bias - how the exponentiated square root of a log transformed squared value is not the same as the original value 
    # sigma.squared is on the log scale.
    x[1] <- x[1]*exp(sigma.squared/2)
    
    # The output from the function nlminb is different from nlmin. Here 0 indicates successful convergence, and have asked for the number of iterations as well.
    return(list(pv=x,sigma=sqrt(sigma.squared),nll=nll,
                evaluations=.Evaluations, converged=nlmin.out$convergence,
                iterations=nlmin.out$iterations))
  }
}

#### hypothetical spawner-recruit data
spawners<-c(56, 68,100,runif(20,90,150))
recruits<-c(76,111,122, runif(20,100,200))
plot(spawners,recruits,xlim=c(0,250),ylim=c(0,250),xlab='Spawners',ylab='Recruits')

newdata<-c(seq(0:250))
BHfit<-ml.bhlnl2(spawners,recruits)
## add the fit
junk<-BHfit$pv[1]*newdata/(newdata/(BHfit$pv[2]/BHfit$pv[1])+1)
lines(newdata,junk,col="red",lty=4)

### approximate carrying capacity implied by the plot. 
abline(h=180)

### historical data was found and the highest number of spawners in one year ever recorded for dusky scallop shark was 800.
## there was no recruit data to go with this number. 
plot(spawners,recruits,xlim=c(0,800),ylim=c(0,800),xlab='Spawners',ylab='Recruits')
junk<-BHfit$pv[1]*newdata/(newdata/(BHfit$pv[2]/BHfit$pv[1])+1)
lines(newdata,junk,col="red",lty=4)

# Question 2a: Does the S-R curve meaningfully describe population size and growth potential?
## This can be argued either way. Give me one answer for 'yes' and one for 'no' [1-2 sentences each]

# Question 2b: If I told you that dusky scallop shark had been fished since recorded history started, but our S-R data was from 2010-2020;
# Does the S-R curve meaningfully describe population size and growth potential? [This one has a correct answer - either yes or no]