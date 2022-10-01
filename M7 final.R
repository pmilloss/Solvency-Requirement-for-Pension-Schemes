#rm(list=ls())

#library(polynom)
#library(StMoMo)

#fortime<-function()
#{
#number of insureds

#install.packages("StMoMo")
#install.packages("polynom")

#library(polynom)
#library(StMoMo)
set.seed(31)
totalruns=1000
myresults<-list()
for (trun in 1:totalruns) {
  
  
  
  n=1000
  
  #Interest rate
  InterestRate=0.03
  
  #number of simulations
  runs=10000
  
  #Annual payment ammount
  R=100
  
  #Entry age
  EntryAge=65
  
  #deferral period
  deferral=5
  
  #Annuity age
  AnnuityAge=EntryAge+deferral
  
  agess <- EntryAge:120
  
  #Here we first fit M7 model to EW date. The code is from StMoMo vignette

  
  EWMaleIniData <- central2initial(EWMaleData)
  ages.fit <- EntryAge:100
  wxt <- genWeightMat(ages = ages.fit, years = EWMaleIniData$years, clip = 3)
  
  M7 <- m7()
  
  M7fit <- fit(M7, data = EWMaleIniData, ages.fit = ages.fit, wxt = wxt)
  
  M7res <- residuals(M7fit)
  
  
  
  M7for <- forecast(M7fit, h = 37, gc.order = c(2, 0, 0))
  
  
  #Here I create the column of mortality rates (q_x) from the central forecasts for the cohort that is 65 in 2012. This is used to caM7ulate premium and reserves
  q_xM7 <- extractCohort(A = M7for$rates, age = 65, period = 2012)
  
  ######## Here I linearly interpolated the  mortality rates from age 100 to age 120 (when q_x=1)
  
  #coefficient<-(1-q_xM71[36,1])/20
  
  # quadratic interpolation
  
  
  quadrinterp <- poly.calc(x = c(99, 100, 120), y = c(q_xM7[35], q_xM7[36], 1))
  q_xM7[36 + (1 : 20)] <- predict(quadrinterp, 101 : 120)
  #points(q_xM7, col = "red")
  
  
  #q_xM73=c()
  #for (i in 2:20) {
  # q_xM73[1]<-q_xM71[36,1]+coefficient
  # q_xM73[i]<-q_xM73[(i-1)]+coefficient
  
  #}
  
  #q_xM72<-as.data.frame(q_xM73)
  #colnames(q_xM72)<-"CentralM7"
  #q_xM7<-rbind(q_xM71,q_xM72)
  
  ######## 
  
  #Survival rates
  p_xM7<-1-q_xM7
  
  
  # Names is differen here. "Central M7"
  names(q_xM7) <- names(p_xM7) <- seq(65, 120)
  
  
  
  tt <- seq(0, length(p_xM7))
  
  # cumulative survival rates (USE q_xM7 AS A VECTOR RATHER THAN A DATA FRAME)
  tp_xM7 <- cumprod(c(1, p_xM7))
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows <- c(rep(0, deferral + 1), rep(R, length.out = length(p_xM7) - deferral))
  # annuity_cashflows <- c(rep(0, deferral), rep(1, length.out = length(p_xM7) - deferral + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  # discount factors
  B0t <- (1 + InterestRate) ^ - tt
  
  names(tp_xM7) <- names(annuity_cashflows) <- names(B0t) <- tt
  
  premium <- n*sum(annuity_cashflows * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V <- vector(mode = "numeric", length = length(tt))
  names(V) <- tt
  V[length(tt)] <- 0 # final reserve
  
  Pt <- c(premium, rep(0, length = length(tt) - 1))
  names(Pt) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V[ti0] <- ((V[ti1] + (annuity_cashflows[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt[ti0]
  }
  
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  M7sim <- simulate(M7fit, nsim = runs, h = 37, gc.order = c(1, 1, 0))
  
  #Here I extract the simulated mortality rates for the cohort that is 65 in 2012 and make quatratic interpolation
  # In some models there are mortality rates which are higher than 1
  M7sim2012<- extractCohort(M7sim$rates, age = EntryAge, period = 2012)
  forinterpolation<-as.matrix(rbind(M7sim2012[35,], M7sim2012[36,], rep(1,runs)))
  quadrinterpsim <- as.polylist(poly.calc(x = c(99, 100, 120), forinterpolation))
  
  forsimulations=matrix(data=NA, nrow=20, ncol=runs)
  for (i in 1:runs) {
    forsimulations[,i]<- predict(quadrinterpsim[[i]], 101 : 120)
  }
  
  rownames(forsimulations)<-101:120
  
  
  
  
  M7cohort2012<-rbind(M7sim2012, forsimulations)
  M7cohort2012[M7cohort2012 > 1] <-1 
  rownames(M7cohort2012)<-65:120
  
  ########
  
  #Stochastic numbers of dead people each year. I randomly generate number of deaths(binomial) 
  survivals=matrix(nrow=57,ncol=runs)
  survivals[1,]<-n
  for (i in 1:56) {
    survivals[(i+1),]<-rbinom(runs,survivals[i,],1-M7cohort2012[i,])
    
  }
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic<-R*V*survivals
  lr1st<-rep(0,runs)
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros<-matrix(0,(deferral+1),runs)
  forPaymentStochastic<-survivals[(deferral+2):(121-EntryAge),]
  colnames(forPaymentStochastic)<-1:runs
  colnames(zeros)<-1:runs
  PaymentsStochastic<-R*rbind(zeros,forPaymentStochastic,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund= matrix(data=NA,nrow=57,ncol=runs)
  Fund[1,]<- rep(premium*(1+InterestRate)^(0),runs)-PaymentsStochastic[1,]
  for (i in 2:57) {
    
    Fund[i,]<-Fund[(i-1),]*(1+InterestRate)-PaymentsStochastic[i,]
    
  } 
  colnames(Fund)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic<-Fund-ActualReservesStochastic[(1):(122-EntryAge),]
  
  #Epsilonss
  epsilon3=0.05
  epsilon2=0.025
  epsilon1=0.01
  
  #Time horizons
  t1=7
  t2=12
  t3=56
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1<-(quantile(SolvencyMarginsStochastic[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2<-(quantile(SolvencyMarginsStochastic[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3<-(quantile(SolvencyMarginsStochastic[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  MyMarginStochasticPerc1<-((quantile(SolvencyMarginsStochastic[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium
  
  
  MyMarginStochasticPerc2<-((quantile(SolvencyMarginsStochastic[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium
  
  
  
  
  MyMarginStochasticPerc3<-((quantile(SolvencyMarginsStochastic[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium
  
  
  
  #deferral2 period
  deferral2=0
  
  #Annuity age
  AnnuityAge2=EntryAge+deferral2
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows2 <- c(rep(0, deferral2 + 1), rep(R, length.out = length(p_xM7) - deferral2))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows2) <- names(B0t) <- tt
  
  
  
  
  premium2 <- n*sum(annuity_cashflows2 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V2 <- vector(mode = "numeric", length = length(tt))
  names(V2) <- tt
  V2[length(tt)] <- 0 # final reserve
  
  Pt2 <- c(premium2, rep(0, length = length(tt) - 1))
  names(Pt2) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V2[ti0] <- ((V2[ti1] + (annuity_cashflows2[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt2[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic2<-R*V2*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros2<-matrix(0,(deferral2+1),runs)
  forPaymentStochastic2<-survivals[(deferral2+2):(121-EntryAge),]
  colnames(forPaymentStochastic2)<-1:runs
  colnames(zeros2)<-1:runs
  PaymentsStochastic2<-R*rbind(zeros2,forPaymentStochastic2,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund2= matrix(data=NA,nrow=57,ncol=runs)
  Fund2[1,]<- rep(premium2*(1+InterestRate)^(0),runs)-PaymentsStochastic2[1,]
  for (i in 2:57) {
    
    Fund2[i,]<-Fund2[(i-1),]*(1+InterestRate)-PaymentsStochastic2[i,]
    
  } 
  colnames(Fund2)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic2<-Fund2-ActualReservesStochastic2[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_2<-(quantile(SolvencyMarginsStochastic2[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_2<-(quantile(SolvencyMarginsStochastic2[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_2<-(quantile(SolvencyMarginsStochastic2[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_2<-((quantile(SolvencyMarginsStochastic2[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium2
  
  
  MyMarginStochasticPerc2_2<-((quantile(SolvencyMarginsStochastic2[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium2
  
  
  
  
  MyMarginStochasticPerc3_2<-((quantile(SolvencyMarginsStochastic2[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium2
  
  #deferral3 period
  deferral3=10
  
  #Annuity age
  AnnuityAge3=EntryAge+deferral3
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows3 <- c(rep(0, deferral3 + 1), rep(R, length.out = length(p_xM7) - deferral3))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows3) <- names(B0t) <- tt
  
  
  
  
  premium3 <- n*sum(annuity_cashflows3 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V3 <- vector(mode = "numeric", length = length(tt))
  names(V3) <- tt
  V3[length(tt)] <- 0 # final reserve
  
  Pt3 <- c(premium2, rep(0, length = length(tt) - 1))
  names(Pt3) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V3[ti0] <- ((V3[ti1] + (annuity_cashflows3[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt3[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic3<-R*V3*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros3<-matrix(0,(deferral3+1),runs)
  forPaymentStochastic3<-survivals[(deferral3+2):(121-EntryAge),]
  colnames(forPaymentStochastic3)<-1:runs
  colnames(zeros3)<-1:runs
  PaymentsStochastic3<-R*rbind(zeros3,forPaymentStochastic3,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund3= matrix(data=NA,nrow=57,ncol=runs)
  Fund3[1,]<- rep(premium3*(1+InterestRate)^(0),runs)-PaymentsStochastic3[1,]
  for (i in 2:57) {
    
    Fund3[i,]<-Fund3[(i-1),]*(1+InterestRate)-PaymentsStochastic3[i,]
    
  } 
  colnames(Fund3)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic3<-Fund3-ActualReservesStochastic3[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_3<-(quantile(SolvencyMarginsStochastic3[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_3<-(quantile(SolvencyMarginsStochastic3[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_3<-(quantile(SolvencyMarginsStochastic3[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_3<-((quantile(SolvencyMarginsStochastic3[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium3
  
  
  MyMarginStochasticPerc2_3<-((quantile(SolvencyMarginsStochastic3[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium3
  
  
  
  
  MyMarginStochasticPerc3_3<-((quantile(SolvencyMarginsStochastic3[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium3
  
  
  #deferral4 period
  deferral4=33
  
  #Annuity age
  AnnuityAge4=EntryAge+deferral4
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows4 <- c(rep(0, deferral4 + 1), rep(R, length.out = length(p_xM7) - deferral4))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows4) <- names(B0t) <- tt
  
  
  
  
  premium4 <- n*sum(annuity_cashflows4 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V4 <- vector(mode = "numeric", length = length(tt))
  names(V4) <- tt
  V4[length(tt)] <- 0 # final reserve
  
  Pt4 <- c(premium4, rep(0, length = length(tt) - 1))
  names(Pt4) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V4[ti0] <- ((V4[ti1] + (annuity_cashflows4[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt4[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic4<-R*V4*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros4<-matrix(0,(deferral4+1),runs)
  forPaymentStochastic4<-survivals[(deferral4+2):(121-EntryAge),]
  colnames(forPaymentStochastic4)<-1:runs
  colnames(zeros4)<-1:runs
  PaymentsStochastic4<-R*rbind(zeros4,forPaymentStochastic4,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund4= matrix(data=NA,nrow=57,ncol=runs)
  Fund4[1,]<- rep(premium4*(1+InterestRate)^(0),runs)-PaymentsStochastic4[1,]
  for (i in 2:57) {
    
    Fund4[i,]<-Fund4[(i-1),]*(1+InterestRate)-PaymentsStochastic4[i,]
    
  } 
  colnames(Fund4)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic4<-Fund4-ActualReservesStochastic4[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_4<-(quantile(SolvencyMarginsStochastic4[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_4<-(quantile(SolvencyMarginsStochastic4[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_4<-(quantile(SolvencyMarginsStochastic4[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_4<-((quantile(SolvencyMarginsStochastic4[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium4
  
  
  MyMarginStochasticPerc2_4<-((quantile(SolvencyMarginsStochastic4[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium4
  
  MyMarginStochasticPerc3_4<-((quantile(SolvencyMarginsStochastic4[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium4
  
  #deferral5 period
  deferral5=15
  
  #Annuity age
  AnnuityAge5=EntryAge+deferral5
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows5 <- c(rep(0, deferral5 + 1), rep(R, length.out = length(p_xM7) - deferral5))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows5) <- names(B0t) <- tt
  
  
  
  
  premium5 <- n*sum(annuity_cashflows5 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V5 <- vector(mode = "numeric", length = length(tt))
  names(V5) <- tt
  V5[length(tt)] <- 0 # final reserve
  
  Pt5 <- c(premium5, rep(0, length = length(tt) - 1))
  names(Pt5) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V5[ti0] <- ((V5[ti1] + (annuity_cashflows5[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt5[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic5<-R*V5*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros5<-matrix(0,(deferral5+1),runs)
  forPaymentStochastic5<-survivals[(deferral5+2):(121-EntryAge),]
  colnames(forPaymentStochastic5)<-1:runs
  colnames(zeros5)<-1:runs
  PaymentsStochastic5<-R*rbind(zeros5,forPaymentStochastic5,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund5= matrix(data=NA,nrow=57,ncol=runs)
  Fund5[1,]<- rep(premium5*(1+InterestRate)^(0),runs)-PaymentsStochastic5[1,]
  for (i in 2:57) {
    
    Fund5[i,]<-Fund5[(i-1),]*(1+InterestRate)-PaymentsStochastic5[i,]
    
  } 
  colnames(Fund5)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic5<-Fund5-ActualReservesStochastic5[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_5<-(quantile(SolvencyMarginsStochastic5[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_5<-(quantile(SolvencyMarginsStochastic5[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_5<-(quantile(SolvencyMarginsStochastic5[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_5<-((quantile(SolvencyMarginsStochastic5[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium5
  
  
  MyMarginStochasticPerc2_5<-((quantile(SolvencyMarginsStochastic5[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium5
  
  
  
  
  MyMarginStochasticPerc3_5<-((quantile(SolvencyMarginsStochastic5[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium5
  
  
  #deferral6 period
  deferral6=20
  
  #Annuity age
  AnnuityAge6=EntryAge+deferral6
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows6 <- c(rep(0, deferral6 + 1), rep(R, length.out = length(p_xM7) - deferral6))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows6) <- names(B0t) <- tt
  
  
  
  
  premium6 <- n*sum(annuity_cashflows6 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V6 <- vector(mode = "numeric", length = length(tt))
  names(V6) <- tt
  V6[length(tt)] <- 0 # final reserve
  
  Pt6 <- c(premium6, rep(0, length = length(tt) - 1))
  names(Pt6) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V6[ti0] <- ((V6[ti1] + (annuity_cashflows6[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt6[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic6<-R*V6*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros6<-matrix(0,(deferral6+1),runs)
  forPaymentStochastic6<-survivals[(deferral6+2):(121-EntryAge),]
  colnames(forPaymentStochastic6)<-1:runs
  colnames(zeros6)<-1:runs
  PaymentsStochastic6<-R*rbind(zeros6,forPaymentStochastic6,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund6= matrix(data=NA,nrow=57,ncol=runs)
  Fund6[1,]<- rep(premium6*(1+InterestRate)^(0),runs)-PaymentsStochastic6[1,]
  for (i in 2:57) {
    
    Fund6[i,]<-Fund6[(i-1),]*(1+InterestRate)-PaymentsStochastic6[i,]
    
  } 
  colnames(Fund6)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic6<-Fund6-ActualReservesStochastic6[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_6<-(quantile(SolvencyMarginsStochastic6[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_6<-(quantile(SolvencyMarginsStochastic6[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_6<-(quantile(SolvencyMarginsStochastic6[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_6<-((quantile(SolvencyMarginsStochastic6[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium6
  
  
  MyMarginStochasticPerc2_6<-((quantile(SolvencyMarginsStochastic6[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium6
  
  
  
  
  MyMarginStochasticPerc3_6<-((quantile(SolvencyMarginsStochastic6[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium6
  
  
  #deferral7 period
  deferral7=25
  #Annuity age
  AnnuityAge7=EntryAge+deferral7
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows7 <- c(rep(0, deferral7 + 1), rep(R, length.out = length(p_xM7) - deferral7))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows7) <- names(B0t) <- tt
  
  
  
  
  premium7 <- n*sum(annuity_cashflows7 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V7 <- vector(mode = "numeric", length = length(tt))
  names(V7) <- tt
  V7[length(tt)] <- 0 # final reserve
  
  Pt7 <- c(premium7, rep(0, length = length(tt) - 1))
  names(Pt7) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V7[ti0] <- ((V7[ti1] + (annuity_cashflows7[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt7[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic7<-R*V7*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros7<-matrix(0,(deferral7+1),runs)
  forPaymentStochastic7<-survivals[(deferral7+2):(121-EntryAge),]
  colnames(forPaymentStochastic7)<-1:runs
  colnames(zeros7)<-1:runs
  PaymentsStochastic7<-R*rbind(zeros7,forPaymentStochastic7,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund7= matrix(data=NA,nrow=57,ncol=runs)
  Fund7[1,]<- rep(premium7*(1+InterestRate)^(0),runs)-PaymentsStochastic7[1,]
  for (i in 2:57) {
    
    Fund7[i,]<-Fund7[(i-1),]*(1+InterestRate)-PaymentsStochastic7[i,]
    
  } 
  colnames(Fund7)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic7<-Fund7-ActualReservesStochastic7[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_7<-(quantile(SolvencyMarginsStochastic7[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_7<-(quantile(SolvencyMarginsStochastic7[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_7<-(quantile(SolvencyMarginsStochastic7[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_7<-((quantile(SolvencyMarginsStochastic7[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium7
  
  
  MyMarginStochasticPerc2_7<-((quantile(SolvencyMarginsStochastic7[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium7
  
  
  
  
  MyMarginStochasticPerc3_7<-((quantile(SolvencyMarginsStochastic7[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium7
  
  
  #deferral8 period
  deferral8=30
  
  #Annuity age
  AnnuityAge8=EntryAge+deferral8
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows8 <- c(rep(0, deferral8 + 1), rep(R, length.out = length(p_xM7) - deferral8))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows8) <- names(B0t) <- tt
  
  
  
  
  premium8 <- n*sum(annuity_cashflows8 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V8 <- vector(mode = "numeric", length = length(tt))
  names(V8) <- tt
  V8[length(tt)] <- 0 # final reserve
  
  Pt8 <- c(premium8, rep(0, length = length(tt) - 1))
  names(Pt8) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V8[ti0] <- ((V8[ti1] + (annuity_cashflows8[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt8[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic8<-R*V8*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros8<-matrix(0,(deferral8+1),runs)
  forPaymentStochastic8<-survivals[(deferral8+2):(121-EntryAge),]
  colnames(forPaymentStochastic8)<-1:runs
  colnames(zeros8)<-1:runs
  PaymentsStochastic8<-R*rbind(zeros8,forPaymentStochastic8,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund8= matrix(data=NA,nrow=57,ncol=runs)
  Fund8[1,]<- rep(premium8*(1+InterestRate)^(0),runs)-PaymentsStochastic8[1,]
  for (i in 2:57) {
    
    Fund8[i,]<-Fund8[(i-1),]*(1+InterestRate)-PaymentsStochastic8[i,]
    
  } 
  colnames(Fund8)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic8<-Fund8-ActualReservesStochastic8[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_8<-(quantile(SolvencyMarginsStochastic8[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_8<-(quantile(SolvencyMarginsStochastic8[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_8<-(quantile(SolvencyMarginsStochastic8[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_8<-((quantile(SolvencyMarginsStochastic8[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium8
  
  
  MyMarginStochasticPerc2_8<-((quantile(SolvencyMarginsStochastic8[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium8
  
  
  
  
  MyMarginStochasticPerc3_8<-((quantile(SolvencyMarginsStochastic8[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium8
  
  #deferral9 period
  deferral9=35
  
  #Annuity age
  AnnuityAge9=EntryAge+deferral9
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows9 <- c(rep(0, deferral9 + 1), rep(R, length.out = length(p_xM7) - deferral9))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows9) <- names(B0t) <- tt
  
  
  
  
  premium9 <- n*sum(annuity_cashflows9 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V9 <- vector(mode = "numeric", length = length(tt))
  names(V9) <- tt
  V9[length(tt)] <- 0 # final reserve
  
  Pt9 <- c(premium9, rep(0, length = length(tt) - 1))
  names(Pt9) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V9[ti0] <- ((V9[ti1] + (annuity_cashflows9[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt9[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic9<-R*V9*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros9<-matrix(0,(deferral9+1),runs)
  forPaymentStochastic9<-survivals[(deferral9+2):(121-EntryAge),]
  colnames(forPaymentStochastic9)<-1:runs
  colnames(zeros9)<-1:runs
  PaymentsStochastic9<-R*rbind(zeros9,forPaymentStochastic9,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund9= matrix(data=NA,nrow=57,ncol=runs)
  Fund9[1,]<- rep(premium9*(1+InterestRate)^(0),runs)-PaymentsStochastic9[1,]
  for (i in 2:57) {
    
    Fund9[i,]<-Fund9[(i-1),]*(1+InterestRate)-PaymentsStochastic9[i,]
    
  } 
  colnames(Fund9)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic9<-Fund9-ActualReservesStochastic9[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_9<-(quantile(SolvencyMarginsStochastic9[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_9<-(quantile(SolvencyMarginsStochastic9[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_9<-(quantile(SolvencyMarginsStochastic9[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_9<-((quantile(SolvencyMarginsStochastic9[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium9
  
  
  MyMarginStochasticPerc2_9<-((quantile(SolvencyMarginsStochastic9[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium9
  
  
  
  
  MyMarginStochasticPerc3_9<-((quantile(SolvencyMarginsStochastic9[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium9
  
  #deferral10 period
  deferral10=40
  
  #Annuity age
  AnnuityAge10=EntryAge+deferral10
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows10 <- c(rep(0, deferral10 + 1), rep(R, length.out = length(p_xM7) - deferral10))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows10) <- names(B0t) <- tt
  
  
  
  
  premium10 <- n*sum(annuity_cashflows10 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V10 <- vector(mode = "numeric", length = length(tt))
  names(V10) <- tt
  V10[length(tt)] <- 0 # final reserve
  
  Pt10 <- c(premium10, rep(0, length = length(tt) - 1))
  names(Pt10) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V10[ti0] <- ((V10[ti1] + (annuity_cashflows10[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt10[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic10<-R*V10*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros10<-matrix(0,(deferral10+1),runs)
  forPaymentStochastic10<-survivals[(deferral10+2):(121-EntryAge),]
  colnames(forPaymentStochastic10)<-1:runs
  colnames(zeros10)<-1:runs
  PaymentsStochastic10<-R*rbind(zeros10,forPaymentStochastic10,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund10= matrix(data=NA,nrow=57,ncol=runs)
  Fund10[1,]<- rep(premium10*(1+InterestRate)^(0),runs)-PaymentsStochastic10[1,]
  for (i in 2:57) {
    
    Fund10[i,]<-Fund10[(i-1),]*(1+InterestRate)-PaymentsStochastic10[i,]
    
  } 
  colnames(Fund10)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic10<-Fund10-ActualReservesStochastic10[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_10<-(quantile(SolvencyMarginsStochastic10[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_10<-(quantile(SolvencyMarginsStochastic10[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_10<-(quantile(SolvencyMarginsStochastic10[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_10<-((quantile(SolvencyMarginsStochastic10[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium10
  
  
  MyMarginStochasticPerc2_10<-((quantile(SolvencyMarginsStochastic10[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium10
  
  
  
  
  MyMarginStochasticPerc3_10<-((quantile(SolvencyMarginsStochastic10[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium10
  
  
  #deferral11 period
  deferral11=45
  
  #Annuity age
  AnnuityAge11=EntryAge+deferral11
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows11 <- c(rep(0, deferral11 + 1), rep(R, length.out = length(p_xM7) - deferral11))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows11) <- names(B0t) <- tt
  
  
  
  
  premium11 <- n*sum(annuity_cashflows11 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V11 <- vector(mode = "numeric", length = length(tt))
  names(V11) <- tt
  V11[length(tt)] <- 0 # final reserve
  
  Pt11 <- c(premium11, rep(0, length = length(tt) - 1))
  names(Pt11) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V11[ti0] <- ((V11[ti1] + (annuity_cashflows11[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt11[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic11<-R*V11*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros11<-matrix(0,(deferral11+1),runs)
  forPaymentStochastic11<-survivals[(deferral11+2):(121-EntryAge),]
  colnames(forPaymentStochastic11)<-1:runs
  colnames(zeros11)<-1:runs
  PaymentsStochastic11<-R*rbind(zeros11,forPaymentStochastic11,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund11= matrix(data=NA,nrow=57,ncol=runs)
  Fund11[1,]<- rep(premium11*(1+InterestRate)^(0),runs)-PaymentsStochastic11[1,]
  for (i in 2:57) {
    
    Fund11[i,]<-Fund11[(i-1),]*(1+InterestRate)-PaymentsStochastic11[i,]
    
  } 
  colnames(Fund11)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic11<-Fund11-ActualReservesStochastic11[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_11<-(quantile(SolvencyMarginsStochastic11[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_11<-(quantile(SolvencyMarginsStochastic11[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_11<-(quantile(SolvencyMarginsStochastic11[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_11<-((quantile(SolvencyMarginsStochastic11[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium11
  
  
  MyMarginStochasticPerc2_11<-((quantile(SolvencyMarginsStochastic11[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium11
  
  
  
  
  MyMarginStochasticPerc3_11<-((quantile(SolvencyMarginsStochastic11[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium11
  
  #deferral14 period
  deferral14=48
  
  #Annuity age
  AnnuityAge14=EntryAge+deferral14
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows14 <- c(rep(0, deferral14 + 1), rep(R, length.out = length(p_xM7) - deferral14))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows14) <- names(B0t) <- tt
  
  
  
  
  premium14 <- n*sum(annuity_cashflows14 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V14 <- vector(mode = "numeric", length = length(tt))
  names(V14) <- tt
  V14[length(tt)] <- 0 # final reserve
  
  Pt14 <- c(premium14, rep(0, length = length(tt) - 1))
  names(Pt14) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V14[ti0] <- ((V14[ti1] + (annuity_cashflows14[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt14[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic14<-R*V14*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros14<-matrix(0,(deferral14+1),runs)
  forPaymentStochastic14<-survivals[(deferral14+2):(121-EntryAge),]
  colnames(forPaymentStochastic14)<-1:runs
  colnames(zeros14)<-1:runs
  PaymentsStochastic14<-R*rbind(zeros14,forPaymentStochastic14,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund14= matrix(data=NA,nrow=57,ncol=runs)
  Fund14[1,]<- rep(premium14*(1+InterestRate)^(0),runs)-PaymentsStochastic14[1,]
  for (i in 2:57) {
    
    Fund14[i,]<-Fund14[(i-1),]*(1+InterestRate)-PaymentsStochastic14[i,]
    
  } 
  colnames(Fund14)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic14<-Fund14-ActualReservesStochastic14[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_14<-(quantile(SolvencyMarginsStochastic14[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_14<-(quantile(SolvencyMarginsStochastic14[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_14<-(quantile(SolvencyMarginsStochastic14[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_14<-((quantile(SolvencyMarginsStochastic14[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium14
  
  
  MyMarginStochasticPerc2_14<-((quantile(SolvencyMarginsStochastic14[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium14
  
  
  
  
  MyMarginStochasticPerc3_14<-((quantile(SolvencyMarginsStochastic14[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium14
  
  
  
  #deferral12 period
  deferral12=50
  
  #Annuity age
  AnnuityAge12=EntryAge+deferral12
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows12 <- c(rep(0, deferral12 + 1), rep(R, length.out = length(p_xM7) - deferral12))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows12) <- names(B0t) <- tt
  
  
  
  
  premium12 <- n*sum(annuity_cashflows12 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V12 <- vector(mode = "numeric", length = length(tt))
  names(V12) <- tt
  V12[length(tt)] <- 0 # final reserve
  
  Pt12 <- c(premium12, rep(0, length = length(tt) - 1))
  names(Pt12) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V12[ti0] <- ((V12[ti1] + (annuity_cashflows12[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt12[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic12<-R*V12*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros12<-matrix(0,(deferral12+1),runs)
  forPaymentStochastic12<-survivals[(deferral12+2):(121-EntryAge),]
  colnames(forPaymentStochastic12)<-1:runs
  colnames(zeros12)<-1:runs
  PaymentsStochastic12<-R*rbind(zeros12,forPaymentStochastic12,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund12= matrix(data=NA,nrow=57,ncol=runs)
  Fund12[1,]<- rep(premium12*(1+InterestRate)^(0),runs)-PaymentsStochastic12[1,]
  for (i in 2:57) {
    
    Fund12[i,]<-Fund12[(i-1),]*(1+InterestRate)-PaymentsStochastic12[i,]
    
  } 
  colnames(Fund12)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic12<-Fund12-ActualReservesStochastic12[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_12<-(quantile(SolvencyMarginsStochastic12[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_12<-(quantile(SolvencyMarginsStochastic12[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_12<-(quantile(SolvencyMarginsStochastic12[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_12<-((quantile(SolvencyMarginsStochastic12[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium12
  
  
  MyMarginStochasticPerc2_12<-((quantile(SolvencyMarginsStochastic12[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium12
  
  
  
  
  MyMarginStochasticPerc3_12<-((quantile(SolvencyMarginsStochastic12[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium12
  #deferral15 period
  deferral15=52
  
  #Annuity age
  AnnuityAge15=EntryAge+deferral15
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows15 <- c(rep(0, deferral15 + 1), rep(R, length.out = length(p_xM7) - deferral15))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows15) <- names(B0t) <- tt
  
  
  
  
  premium15 <- n*sum(annuity_cashflows15 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V15 <- vector(mode = "numeric", length = length(tt))
  names(V15) <- tt
  V15[length(tt)] <- 0 # final reserve
  
  Pt15 <- c(premium15, rep(0, length = length(tt) - 1))
  names(Pt15) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V15[ti0] <- ((V15[ti1] + (annuity_cashflows15[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt15[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic15<-R*V15*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros15<-matrix(0,(deferral15+1),runs)
  forPaymentStochastic15<-survivals[(deferral15+2):(121-EntryAge),]
  colnames(forPaymentStochastic15)<-1:runs
  colnames(zeros15)<-1:runs
  PaymentsStochastic15<-R*rbind(zeros15,forPaymentStochastic15,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund15= matrix(data=NA,nrow=57,ncol=runs)
  Fund15[1,]<- rep(premium15*(1+InterestRate)^(0),runs)-PaymentsStochastic15[1,]
  for (i in 2:57) {
    
    Fund15[i,]<-Fund15[(i-1),]*(1+InterestRate)-PaymentsStochastic15[i,]
    
  } 
  colnames(Fund15)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic15<-Fund15-ActualReservesStochastic15[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_15<-(quantile(SolvencyMarginsStochastic15[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_15<-(quantile(SolvencyMarginsStochastic15[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_15<-(quantile(SolvencyMarginsStochastic15[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_15<-((quantile(SolvencyMarginsStochastic15[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium15
  
  
  MyMarginStochasticPerc2_15<-((quantile(SolvencyMarginsStochastic15[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium15
  
  
  
  
  MyMarginStochasticPerc3_15<-((quantile(SolvencyMarginsStochastic15[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium15
  
  
  
  
  
  
  #deferral13 period
  deferral13=55
  
  #Annuity age
  AnnuityAge13=EntryAge+deferral13
  
  
  
  
  # cashflows at times 0, 1, 2, ... (I DON'T REMEMBER IF IT WAS IMMEDIATE OR DUE :) THIS ONE BELOW IS IMMEDIATE)
  annuity_cashflows13 <- c(rep(0, deferral13 + 1), rep(R, length.out = length(p_xM7) - deferral13))
  # annuity_cashflows <- c(rep(0, deferral2), rep(1, length.out = length(p_xM7) - deferral2 + 1)) # THIS IF YOU WANT ANNUITY DUE
  
  
  names(tp_xM7) <- names(annuity_cashflows13) <- names(B0t) <- tt
  
  
  
  
  premium13 <- n*sum(annuity_cashflows13 * tp_xM7 * B0t)
  
  
  
  ############################
  # CaM7ulation backward #####
  ############################
  V13 <- vector(mode = "numeric", length = length(tt))
  names(V13) <- tt
  V13[length(tt)] <- 0 # final reserve
  
  Pt13 <- c(premium13, rep(0, length = length(tt) - 1))
  names(Pt13) <- tt
  
  for (i in tail(tt, 1) : 1)
  {
    ti1 <- as.character(i)
    ti0 <- as.character(i - 1)
    xt <- as.character(EntryAge + i - 1)
    
    # caM7ulate V(t-1) from V(t)
    V13[ti0] <- ((V13[ti1] + (annuity_cashflows13[ti1])/R) * p_xM7[xt]) / (1 + InterestRate) - Pt13[ti0]
  }
  
  
  #Here I simulate (1000 simulations) future mortality rates
  
  #Stochastic reserves based on the number of alives
  ActualReservesStochastic13<-R*V13*survivals
  
  #CaM7ulate random payments to policyholders considering random number of alives
  zeros13<-matrix(0,(deferral13+1),runs)
  forPaymentStochastic13<-survivals[(deferral13+2):(121-EntryAge),]
  colnames(forPaymentStochastic13)<-1:runs
  colnames(zeros13)<-1:runs
  PaymentsStochastic13<-R*rbind(zeros13,forPaymentStochastic13,lr1st)
  
  #CaM7ulate random portfolio fund at times t, i.e. Z_t
  
  Fund13= matrix(data=NA,nrow=57,ncol=runs)
  Fund13[1,]<- rep(premium13*(1+InterestRate)^(0),runs)-PaymentsStochastic13[1,]
  for (i in 2:57) {
    
    Fund13[i,]<-Fund13[(i-1),]*(1+InterestRate)-PaymentsStochastic13[i,]
    
  } 
  colnames(Fund13)<-1:runs
  
  #Matrix of deficits in case no solvency margin is added at the beginning
  SolvencyMarginsStochastic13<-Fund13-ActualReservesStochastic13[(1):(122-EntryAge),]
  
  
  #CaM7ulate solvency margin that needs to be added in the beginning in order to ensure solvence at time t. This is basically a quantile of dicounted deficit
  MyMarginStochastic1_13<-(quantile(SolvencyMarginsStochastic13[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1)
  
  
  MyMarginStochastic2_13<-(quantile(SolvencyMarginsStochastic13[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1)
  
  
  MyMarginStochastic3_13<-(quantile(SolvencyMarginsStochastic13[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1)
  
  
  
  MyMarginStochasticPerc1_13<-((quantile(SolvencyMarginsStochastic13[t1,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t1-1))/premium13
  
  
  MyMarginStochasticPerc2_13<-((quantile(SolvencyMarginsStochastic13[t2,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t2-1))/premium13
  
  
  
  
  MyMarginStochasticPerc3_13<-((quantile(SolvencyMarginsStochastic13[t3,],c(epsilon1,epsilon2,epsilon3)))/(1+InterestRate)^(t3-1))/premium13
  
  
  
  
  ES1<-c((mean(SolvencyMarginsStochastic[t1,][rank(SolvencyMarginsStochastic[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium,(mean(SolvencyMarginsStochastic[t1,][rank(SolvencyMarginsStochastic[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium,(mean(SolvencyMarginsStochastic[t1,][rank(SolvencyMarginsStochastic[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium)
  ES2<-c((mean(SolvencyMarginsStochastic[t2,][rank(SolvencyMarginsStochastic[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium,(mean(SolvencyMarginsStochastic[t2,][rank(SolvencyMarginsStochastic[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium,(mean(SolvencyMarginsStochastic[t2,][rank(SolvencyMarginsStochastic[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium)
  ES3<-c((mean(SolvencyMarginsStochastic[t3,][rank(SolvencyMarginsStochastic[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium,(mean(SolvencyMarginsStochastic[t3,][rank(SolvencyMarginsStochastic[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium,(mean(SolvencyMarginsStochastic[t3,][rank(SolvencyMarginsStochastic[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium)
  ES2_1<-c((mean(SolvencyMarginsStochastic2[t1,][rank(SolvencyMarginsStochastic2[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium2,(mean(SolvencyMarginsStochastic2[t1,][rank(SolvencyMarginsStochastic2[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium2,(mean(SolvencyMarginsStochastic2[t1,][rank(SolvencyMarginsStochastic2[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium2)
  ES2_2<-c((mean(SolvencyMarginsStochastic2[t2,][rank(SolvencyMarginsStochastic2[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium2,(mean(SolvencyMarginsStochastic2[t2,][rank(SolvencyMarginsStochastic2[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium2,(mean(SolvencyMarginsStochastic2[t2,][rank(SolvencyMarginsStochastic2[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium2)
  ES2_3<-c((mean(SolvencyMarginsStochastic2[t3,][rank(SolvencyMarginsStochastic2[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium2,(mean(SolvencyMarginsStochastic2[t3,][rank(SolvencyMarginsStochastic2[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium2,(mean(SolvencyMarginsStochastic2[t3,][rank(SolvencyMarginsStochastic2[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium2)
  ES3_1<-c((mean(SolvencyMarginsStochastic3[t1,][rank(SolvencyMarginsStochastic3[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium3,(mean(SolvencyMarginsStochastic3[t1,][rank(SolvencyMarginsStochastic3[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium3,(mean(SolvencyMarginsStochastic3[t1,][rank(SolvencyMarginsStochastic3[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium3)
  ES3_2<-c((mean(SolvencyMarginsStochastic3[t2,][rank(SolvencyMarginsStochastic3[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium3,(mean(SolvencyMarginsStochastic3[t2,][rank(SolvencyMarginsStochastic3[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium3,(mean(SolvencyMarginsStochastic3[t2,][rank(SolvencyMarginsStochastic3[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium3)
  ES3_3<-c((mean(SolvencyMarginsStochastic3[t3,][rank(SolvencyMarginsStochastic3[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium3,(mean(SolvencyMarginsStochastic3[t3,][rank(SolvencyMarginsStochastic3[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium3,(mean(SolvencyMarginsStochastic3[t3,][rank(SolvencyMarginsStochastic3[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium3)
  ES4_1<-c((mean(SolvencyMarginsStochastic4[t1,][rank(SolvencyMarginsStochastic4[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium4,(mean(SolvencyMarginsStochastic4[t1,][rank(SolvencyMarginsStochastic4[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium4,(mean(SolvencyMarginsStochastic4[t1,][rank(SolvencyMarginsStochastic4[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium4)
  ES4_2<-c((mean(SolvencyMarginsStochastic4[t2,][rank(SolvencyMarginsStochastic4[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium4,(mean(SolvencyMarginsStochastic4[t2,][rank(SolvencyMarginsStochastic4[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium4,(mean(SolvencyMarginsStochastic4[t2,][rank(SolvencyMarginsStochastic4[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium4)
  ES4_3<-c((mean(SolvencyMarginsStochastic4[t3,][rank(SolvencyMarginsStochastic4[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium4,(mean(SolvencyMarginsStochastic4[t3,][rank(SolvencyMarginsStochastic4[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium4,(mean(SolvencyMarginsStochastic4[t3,][rank(SolvencyMarginsStochastic4[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium4)
  ES5_1<-c((mean(SolvencyMarginsStochastic5[t1,][rank(SolvencyMarginsStochastic5[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium5,(mean(SolvencyMarginsStochastic5[t1,][rank(SolvencyMarginsStochastic5[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium5,(mean(SolvencyMarginsStochastic5[t1,][rank(SolvencyMarginsStochastic5[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium5)
  ES5_2<-c((mean(SolvencyMarginsStochastic5[t2,][rank(SolvencyMarginsStochastic5[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium5,(mean(SolvencyMarginsStochastic5[t2,][rank(SolvencyMarginsStochastic5[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium5,(mean(SolvencyMarginsStochastic5[t2,][rank(SolvencyMarginsStochastic5[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium5)
  ES5_3<-c((mean(SolvencyMarginsStochastic5[t3,][rank(SolvencyMarginsStochastic5[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium5,(mean(SolvencyMarginsStochastic5[t3,][rank(SolvencyMarginsStochastic5[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium5,(mean(SolvencyMarginsStochastic5[t3,][rank(SolvencyMarginsStochastic5[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium5)
  ES6_1<-c((mean(SolvencyMarginsStochastic6[t1,][rank(SolvencyMarginsStochastic6[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium6,(mean(SolvencyMarginsStochastic6[t1,][rank(SolvencyMarginsStochastic6[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium6,(mean(SolvencyMarginsStochastic6[t1,][rank(SolvencyMarginsStochastic6[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium6)
  ES6_2<-c((mean(SolvencyMarginsStochastic6[t2,][rank(SolvencyMarginsStochastic6[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium6,(mean(SolvencyMarginsStochastic6[t2,][rank(SolvencyMarginsStochastic6[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium6,(mean(SolvencyMarginsStochastic6[t2,][rank(SolvencyMarginsStochastic6[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium6)
  ES6_3<-c((mean(SolvencyMarginsStochastic6[t3,][rank(SolvencyMarginsStochastic6[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium6,(mean(SolvencyMarginsStochastic6[t3,][rank(SolvencyMarginsStochastic6[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium6,(mean(SolvencyMarginsStochastic6[t3,][rank(SolvencyMarginsStochastic6[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium6)
  ES7_1<-c((mean(SolvencyMarginsStochastic7[t1,][rank(SolvencyMarginsStochastic7[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium7,(mean(SolvencyMarginsStochastic7[t1,][rank(SolvencyMarginsStochastic7[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium7,(mean(SolvencyMarginsStochastic7[t1,][rank(SolvencyMarginsStochastic7[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium7)
  ES7_2<-c((mean(SolvencyMarginsStochastic7[t2,][rank(SolvencyMarginsStochastic7[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium7,(mean(SolvencyMarginsStochastic7[t2,][rank(SolvencyMarginsStochastic7[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium7,(mean(SolvencyMarginsStochastic7[t2,][rank(SolvencyMarginsStochastic7[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium7)
  ES7_3<-c((mean(SolvencyMarginsStochastic7[t3,][rank(SolvencyMarginsStochastic7[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium7,(mean(SolvencyMarginsStochastic7[t3,][rank(SolvencyMarginsStochastic7[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium7,(mean(SolvencyMarginsStochastic7[t3,][rank(SolvencyMarginsStochastic7[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium7)
  ES8_1<-c((mean(SolvencyMarginsStochastic8[t1,][rank(SolvencyMarginsStochastic8[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium8,(mean(SolvencyMarginsStochastic8[t1,][rank(SolvencyMarginsStochastic8[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium8,(mean(SolvencyMarginsStochastic8[t1,][rank(SolvencyMarginsStochastic8[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium8)
  ES8_2<-c((mean(SolvencyMarginsStochastic8[t2,][rank(SolvencyMarginsStochastic8[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium8,(mean(SolvencyMarginsStochastic8[t2,][rank(SolvencyMarginsStochastic8[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium8,(mean(SolvencyMarginsStochastic8[t2,][rank(SolvencyMarginsStochastic8[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium8)
  ES8_3<-c((mean(SolvencyMarginsStochastic8[t3,][rank(SolvencyMarginsStochastic8[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium8,(mean(SolvencyMarginsStochastic8[t3,][rank(SolvencyMarginsStochastic8[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium8,(mean(SolvencyMarginsStochastic8[t3,][rank(SolvencyMarginsStochastic8[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium8)
  ES9_1<-c((mean(SolvencyMarginsStochastic9[t1,][rank(SolvencyMarginsStochastic9[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium9,(mean(SolvencyMarginsStochastic9[t1,][rank(SolvencyMarginsStochastic9[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium9,(mean(SolvencyMarginsStochastic9[t1,][rank(SolvencyMarginsStochastic9[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium9)
  ES9_2<-c((mean(SolvencyMarginsStochastic9[t2,][rank(SolvencyMarginsStochastic9[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium9,(mean(SolvencyMarginsStochastic9[t2,][rank(SolvencyMarginsStochastic9[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium9,(mean(SolvencyMarginsStochastic9[t2,][rank(SolvencyMarginsStochastic9[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium9)
  ES9_3<-c((mean(SolvencyMarginsStochastic9[t3,][rank(SolvencyMarginsStochastic9[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium9,(mean(SolvencyMarginsStochastic9[t3,][rank(SolvencyMarginsStochastic9[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium9,(mean(SolvencyMarginsStochastic9[t3,][rank(SolvencyMarginsStochastic9[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium9)
  ES10_1<-c((mean(SolvencyMarginsStochastic10[t1,][rank(SolvencyMarginsStochastic10[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium10,(mean(SolvencyMarginsStochastic10[t1,][rank(SolvencyMarginsStochastic10[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium10,(mean(SolvencyMarginsStochastic10[t1,][rank(SolvencyMarginsStochastic10[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium10)
  ES10_2<-c((mean(SolvencyMarginsStochastic10[t2,][rank(SolvencyMarginsStochastic10[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium10,(mean(SolvencyMarginsStochastic10[t2,][rank(SolvencyMarginsStochastic10[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium10,(mean(SolvencyMarginsStochastic10[t2,][rank(SolvencyMarginsStochastic10[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium10)
  ES10_3<-c((mean(SolvencyMarginsStochastic10[t3,][rank(SolvencyMarginsStochastic10[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium10,(mean(SolvencyMarginsStochastic10[t3,][rank(SolvencyMarginsStochastic10[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium10,(mean(SolvencyMarginsStochastic10[t3,][rank(SolvencyMarginsStochastic10[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium10)
  ES11_1<-c((mean(SolvencyMarginsStochastic11[t1,][rank(SolvencyMarginsStochastic11[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium11,(mean(SolvencyMarginsStochastic11[t1,][rank(SolvencyMarginsStochastic11[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium11,(mean(SolvencyMarginsStochastic11[t1,][rank(SolvencyMarginsStochastic11[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium11)
  ES11_2<-c((mean(SolvencyMarginsStochastic11[t2,][rank(SolvencyMarginsStochastic11[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium11,(mean(SolvencyMarginsStochastic11[t2,][rank(SolvencyMarginsStochastic11[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium11,(mean(SolvencyMarginsStochastic11[t2,][rank(SolvencyMarginsStochastic11[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium11)
  ES11_3<-c((mean(SolvencyMarginsStochastic11[t3,][rank(SolvencyMarginsStochastic11[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium11,(mean(SolvencyMarginsStochastic11[t3,][rank(SolvencyMarginsStochastic11[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium11,(mean(SolvencyMarginsStochastic11[t3,][rank(SolvencyMarginsStochastic11[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium11)
  ES14_1<-c((mean(SolvencyMarginsStochastic14[t1,][rank(SolvencyMarginsStochastic14[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium14,(mean(SolvencyMarginsStochastic14[t1,][rank(SolvencyMarginsStochastic14[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium14,(mean(SolvencyMarginsStochastic14[t1,][rank(SolvencyMarginsStochastic14[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium14)
  ES14_2<-c((mean(SolvencyMarginsStochastic14[t2,][rank(SolvencyMarginsStochastic14[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium14,(mean(SolvencyMarginsStochastic14[t2,][rank(SolvencyMarginsStochastic14[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium14,(mean(SolvencyMarginsStochastic14[t2,][rank(SolvencyMarginsStochastic14[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium14)
  ES14_3<-c((mean(SolvencyMarginsStochastic14[t3,][rank(SolvencyMarginsStochastic14[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium14,(mean(SolvencyMarginsStochastic14[t3,][rank(SolvencyMarginsStochastic14[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium14,(mean(SolvencyMarginsStochastic14[t3,][rank(SolvencyMarginsStochastic14[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium14)
  ES12_1<-c((mean(SolvencyMarginsStochastic12[t1,][rank(SolvencyMarginsStochastic12[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium12,(mean(SolvencyMarginsStochastic12[t1,][rank(SolvencyMarginsStochastic12[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium12,(mean(SolvencyMarginsStochastic12[t1,][rank(SolvencyMarginsStochastic12[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium12)
  ES12_2<-c((mean(SolvencyMarginsStochastic12[t2,][rank(SolvencyMarginsStochastic12[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium12,(mean(SolvencyMarginsStochastic12[t2,][rank(SolvencyMarginsStochastic12[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium12,(mean(SolvencyMarginsStochastic12[t2,][rank(SolvencyMarginsStochastic12[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium12)
  ES12_3<-c((mean(SolvencyMarginsStochastic12[t3,][rank(SolvencyMarginsStochastic12[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium12,(mean(SolvencyMarginsStochastic12[t3,][rank(SolvencyMarginsStochastic12[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium12,(mean(SolvencyMarginsStochastic12[t3,][rank(SolvencyMarginsStochastic12[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium12)
  ES15_1<-c((mean(SolvencyMarginsStochastic15[t1,][rank(SolvencyMarginsStochastic15[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium15,(mean(SolvencyMarginsStochastic15[t1,][rank(SolvencyMarginsStochastic15[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium15,(mean(SolvencyMarginsStochastic15[t1,][rank(SolvencyMarginsStochastic15[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium15)
  ES15_2<-c((mean(SolvencyMarginsStochastic15[t2,][rank(SolvencyMarginsStochastic15[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium15,(mean(SolvencyMarginsStochastic15[t2,][rank(SolvencyMarginsStochastic15[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium15,(mean(SolvencyMarginsStochastic15[t2,][rank(SolvencyMarginsStochastic15[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium15)
  ES15_3<-c((mean(SolvencyMarginsStochastic15[t3,][rank(SolvencyMarginsStochastic15[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium15,(mean(SolvencyMarginsStochastic15[t3,][rank(SolvencyMarginsStochastic15[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium15,(mean(SolvencyMarginsStochastic15[t3,][rank(SolvencyMarginsStochastic15[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium15)
  ES13_1<-c((mean(SolvencyMarginsStochastic13[t1,][rank(SolvencyMarginsStochastic13[t1,]) <= 100])/(1+InterestRate)^(t1-1))/premium13,(mean(SolvencyMarginsStochastic13[t1,][rank(SolvencyMarginsStochastic13[t1,]) <= 250])/(1+InterestRate)^(t1-1))/premium13,(mean(SolvencyMarginsStochastic13[t1,][rank(SolvencyMarginsStochastic13[t1,]) <= 500])/(1+InterestRate)^(t1-1))/premium13)
  ES13_2<-c((mean(SolvencyMarginsStochastic13[t2,][rank(SolvencyMarginsStochastic13[t2,]) <= 100])/(1+InterestRate)^(t2-1))/premium13,(mean(SolvencyMarginsStochastic13[t2,][rank(SolvencyMarginsStochastic13[t2,]) <= 250])/(1+InterestRate)^(t2-1))/premium13,(mean(SolvencyMarginsStochastic13[t2,][rank(SolvencyMarginsStochastic13[t2,]) <= 500])/(1+InterestRate)^(t2-1))/premium13)
  ES13_3<-c((mean(SolvencyMarginsStochastic13[t3,][rank(SolvencyMarginsStochastic13[t3,]) <= 100])/(1+InterestRate)^(t3-1))/premium13,(mean(SolvencyMarginsStochastic13[t3,][rank(SolvencyMarginsStochastic13[t3,]) <= 250])/(1+InterestRate)^(t3-1))/premium13,(mean(SolvencyMarginsStochastic13[t3,][rank(SolvencyMarginsStochastic13[t3,]) <= 500])/(1+InterestRate)^(t3-1))/premium13)
  
  
  
  
  MyfinalresultsM7<-c(MyMarginStochastic1,MyMarginStochastic2,MyMarginStochastic3,MyMarginStochasticPerc1,MyMarginStochasticPerc2,MyMarginStochasticPerc3,MyMarginStochastic1_2,MyMarginStochastic2_2,MyMarginStochastic3_2,MyMarginStochasticPerc1_2,MyMarginStochasticPerc2_2,MyMarginStochasticPerc3_2,MyMarginStochastic1_3,MyMarginStochastic2_3,MyMarginStochastic3_3,MyMarginStochasticPerc1_3,MyMarginStochasticPerc2_3,MyMarginStochasticPerc3_3,
                        MyMarginStochastic1_4,MyMarginStochastic2_4,MyMarginStochastic3_4,MyMarginStochasticPerc1_4,MyMarginStochasticPerc2_4,MyMarginStochasticPerc3_4,MyMarginStochastic1_5,MyMarginStochastic2_5,MyMarginStochastic3_5,MyMarginStochasticPerc1_5,MyMarginStochasticPerc2_5,MyMarginStochasticPerc3_5,
                        MyMarginStochastic1_6,MyMarginStochastic2_6,MyMarginStochastic3_6,MyMarginStochasticPerc1_6,MyMarginStochasticPerc2_6,MyMarginStochasticPerc3_6,
                        MyMarginStochastic1_7,MyMarginStochastic2_7,MyMarginStochastic3_7,MyMarginStochasticPerc1_7,MyMarginStochasticPerc2_7,MyMarginStochasticPerc3_7,
                        MyMarginStochastic1_8,MyMarginStochastic2_8,MyMarginStochastic3_8,MyMarginStochasticPerc1_8,MyMarginStochasticPerc2_8,MyMarginStochasticPerc3_8,
                        MyMarginStochastic1_9,MyMarginStochastic2_9,MyMarginStochastic3_9,MyMarginStochasticPerc1_9,MyMarginStochasticPerc2_9,MyMarginStochasticPerc3_9,
                        MyMarginStochastic1_10,MyMarginStochastic2_10,MyMarginStochastic3_10,MyMarginStochasticPerc1_10,MyMarginStochasticPerc2_10,MyMarginStochasticPerc3_10,
                        MyMarginStochastic1_11,MyMarginStochastic2_11,MyMarginStochastic3_11,MyMarginStochasticPerc1_11,MyMarginStochasticPerc2_11,MyMarginStochasticPerc3_11,
                        MyMarginStochastic1_14,MyMarginStochastic2_14,MyMarginStochastic3_14,MyMarginStochasticPerc1_14,MyMarginStochasticPerc2_14,MyMarginStochasticPerc3_14,
                        MyMarginStochastic1_12,MyMarginStochastic2_12,MyMarginStochastic3_12,MyMarginStochasticPerc1_12,MyMarginStochasticPerc2_12,MyMarginStochasticPerc3_12,
                        MyMarginStochastic1_15,MyMarginStochastic2_15,MyMarginStochastic3_15,MyMarginStochasticPerc1_15,MyMarginStochasticPerc2_15,MyMarginStochasticPerc3_15,
                        MyMarginStochastic1_13,MyMarginStochastic2_13,MyMarginStochastic3_13,MyMarginStochasticPerc1_13,MyMarginStochasticPerc2_13,MyMarginStochasticPerc3_13,
                        ES1,ES2,ES3,ES2_1,ES2_2,ES2_3,ES3_1,ES3_2,ES3_3,ES4_1,ES4_2,ES4_3,ES5_1,ES5_2,ES5_3,ES6_1,ES6_2,ES6_3,ES7_1,ES7_2,ES7_3,ES8_1,ES8_2,ES8_3,
                        ES9_1,ES9_2,ES9_3,ES10_1,ES10_2,ES10_3,ES11_1,ES11_2,ES11_3,ES14_1,ES14_2,ES14_3,ES12_1,ES12_2,ES12_3,ES15_1,ES15_2,ES15_3,ES13_1,ES13_2,ES13_3)
  
  
  
  myresults[[trun]]<-MyfinalresultsM7
  
}
outputM7<-matrix(unlist(myresults),ncol=405,byrow = TRUE)
write.table(rbind(c(1:405),outputM7),file="OutputM710002018.csv",append=TRUE,sep=",",col.names = FALSE,row.names = FALSE)
#dim(output)
#}
#system.time(fortime())

