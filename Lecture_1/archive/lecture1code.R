##### 

### A Taste of Fishery Science
## Code to support Lecture 1
#### H. Bowlby; Jan 2022

setwd('C:/mydocs/_courses/DALmodule')



###---------------------------------------------------------

## the relationship between M and annual survival

##----------------------------------------------------------

### annual survival scenarios - assuming constant survival rates
### survival  = (# animals alive at the end of the time period) / (# animals alive at the beginning); S = Nt/Nt-1

high.surv<-0.95
low.surv<-0.72
med.surv<-0.88

lx<-c()
lx.low<-c()
lx.med<-c()

x<-0:100
lx[1]<-1
lx.low[1]<-1
lx.med[1]<-1
for(i in 2:(length(x)+1))
{lx[i]<-lx[i-1]*high.surv
lx.low[i]<-lx.low[i-1]*low.surv
lx.med[i]<-lx.med[i-1]*med.surv
}

survival<-data.frame(age=c(x,101),lx,lx.med,lx.low)
### compare the three scenarios - proportion of the original cohort surviving at age

plot(survival$age,survival$lx,type = 'l',xlab='Age',ylab='Survival') # black - high survival
points(survival$age,survival$lx.med,type='l',col='blue',lty=2) # blue - medium survival
points(survival$age,survival$lx.low,type='l',col='red',lty=3) # red - low survival

### comparison of lifespans - approximate age when ~ 10% of the original population remains.

abline(v=min(survival$age[survival$lx<0.01]))
no.F=min(survival$age[survival$lx<0.01])
    
abline(v=min(survival$age[survival$lx.med<0.01]),col='blue',lty=2)
abline(v=min(survival$age[survival$lx.low<0.01]),col='red',lty=3)

#### why are we talking about survival anyways? I thought we cared about mortality.
### in the high survival scenario: 5% mortality
### cumulative mortality is NOT simply 5% multiplied by the number of years

WRONG<-0.05*100  ### implies that more animals die than are available in the first place
CORRECT<-0.95^100  ### check: survival$lx[survival$age==100]


#### what is the probability of death over a very short period of time?
# Instantaneous natural mortality (M)

M<-log(1-0.05)  ## NOTE: same as M=log(survival) ;  annual mortality = 1-exp(M)   ### often written: 1-exp(-M);  1-exp(-0.05129329)=0.05
M.med<-log(1-0.12)  # medium survival scenario
M.low<-log(1-0.28)  ## low survival scenario

### NOTE that M estimates are ALWAYS negative


# QUESTION 1a: What would a positive value for M imply for survival (i.e. how is abundance at changing as a cohort ages)? [one sentence answer]
# QUESTION 1b: Does a long-lived species like Greenland shark experience higher or lower natural mortality than a short-lived species like alewife? [one sentence answer]

######



##--------------------------------------------------------------

# fishing mortality

###-------------------------------------------------------------

## M is instantaneous natural mortality; F is instantaneous fishing mortality; u is annual exploitation (% of population harvested each year)

## comparison of different sizes of fishery
u.low=0.01 # 1% exploitation rate
u.med<-0.1 # 10% exploitation rate
u.high=0.2 # 20% exploitation rate

F.low<--log(1-u.low) 
F.med<--log(1-u.med)
F.high<--log(1-u.high)

M<--log(1-0.05)

### NOTE: take the negative values for M and F (i.e. converting them to positive numbers) to make it easier to add in the following
high.surv<-exp(-(M+F.low))
med.surv<-exp(-(M+F.med))
low.surv<-exp(-(M+F.high))

### redo your plots - now comparing different levels of fishing. 
### the value for M is the same, the value for F is changing. 

lx<-c()

x<-0:100
lx[1]<-1
lx.low[1]<-1
lx.med[1]<-1
for(i in 2:(length(x)+1))
{lx[i]<-lx[i-1]*high.surv
lx.low[i]<-lx.low[i-1]*low.surv
lx.med[i]<-lx.med[i-1]*med.surv
}

survival<-data.frame(age=c(x,101),lx,lx.med,lx.low)

plot(survival$age,survival$lx,type = 'l',xlab='Age',ylab='Survival')
points(survival$age,survival$lx.med,type='l',col='blue',lty=2)
points(survival$age,survival$lx.low,type='l',col='red',lty=3)

abline(v=min(survival$age[survival$lx<0.01]))
abline(v=min(survival$age[survival$lx.med<0.01]),col='blue',lty=2)
abline(v=min(survival$age[survival$lx.low<0.01]),col='red',lty=3)

### compare with original high survival scenario (low M) and no F
abline(v=no.F, lwd=2)


# QUESTION 2: 
# What would you expect to happen to the age distribution of the catches (the number of animals at each age) if fishing mortality is very high?
# [1-2 sentences]

### notice that Z (total mortality) = M + F


###-----------------------------------------------------------

## Reproductive output at age

###-----------------------------------------------------------

### hypothetical life history parameters

age.maturity = 8
litter.size = 20 
sex.ratio = 0.5  ### the following code requires reproductive output of female pups
gestation = 2 ## in years
max.age = 40

# what is reproductive output at age?

temp<-c(0,c(rep(0,age.maturity+gestation-1),litter.size*sex.ratio,rep(c(litter.size*sex.ratio),max.age))) 
mx<-temp[1:(max.age+1)]

plot(mx,type='l',xlab='Age',ylab='Reproductive output (pups)')
sum(mx)


## BONUS: explore how female reproductive output might change for a species with a shorter/longer lifespan or gestation period [no answer required]
### aka - change the values above. 

###--------------------------------------------------------

## population growth with and without fishing

### combining survival, reproductive output and Fishing

###--------------------------------------------------------

lotka.r<-function(age.maturity,litter.size,sex.ratio,gestation,max.age,u,sel) 
{
  
  temp<-c(0,c(rep(0,age.maturity+gestation-1),litter.size*sex.ratio,rep(c(litter.size*sex.ratio),max.age))) 
  mx<-temp[1:(max.age+1)]

  # partition survival into natural and fishing mortality  
  M<-rep(2.56*max.age^-0.873,max.age)  ### just trust me - you can approximate M from longevity data for sharks.
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


## QUESTION 3: Would it be meaningful to do a stock assessment on this species given our understanding of life history?
###             Why or why not? [3-4 sentences]












