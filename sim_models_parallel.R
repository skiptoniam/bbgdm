library(bbgdm)
form <- y~1+x
nsp <- 250
Nsite <- 25
ndissim <- (Nsite*(Nsite-1))/2
dat <- data.frame(y=rep(1,ndissim),x=seq(0,4,length.out = ndissim))#imaginary difference in a covariate
set.seed(42)
sim1 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),.1),1,2),nsp)
sim2 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),.5),1,2),nsp)
sim3 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),1),1,2),nsp)
sim4 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),2),1,2),nsp)
sim5 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),5),1,2),nsp)

sim6 <- simulate_dism_data(form,dat, matrix(c(qlogis(.1),1),1,2),nsp)
sim7 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),1),1,2),nsp)
sim8 <- simulate_dism_data(form,dat, matrix(c(qlogis(.5),1),1,2),nsp)
sim9 <- simulate_dism_data(form,dat, matrix(c(qlogis(.65),1),1,2),nsp)
sim10 <- simulate_dism_data(form,dat, matrix(c(qlogis(.8),1),1,2),nsp)

form <- y~1+x
nsp <- 250
Nsite <- 50
ndissim <- (Nsite*(Nsite-1))/2
dat <- data.frame(y=rep(1,ndissim),x=seq(0,4,length.out = ndissim))#imaginary difference in a covariate
set.seed(42)
sim11 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),.1),1,2),nsp)
sim12 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),.5),1,2),nsp)
sim13 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),1),1,2),nsp)
sim14 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),2),1,2),nsp)
sim15 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),5),1,2),nsp)

sim16 <- simulate_dism_data(form,dat, matrix(c(qlogis(.1),1),1,2),nsp)
sim17 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),1),1,2),nsp)
sim18 <- simulate_dism_data(form,dat, matrix(c(qlogis(.5),1),1,2),nsp)
sim19 <- simulate_dism_data(form,dat, matrix(c(qlogis(.65),1),1,2),nsp)
sim20 <- simulate_dism_data(form,dat, matrix(c(qlogis(.8),1),1,2),nsp)

form <- y~1+x
nsp <- 250
Nsite <- 100
ndissim <- (Nsite*(Nsite-1))/2
dat <- data.frame(y=rep(1,ndissim),x=seq(0,4,length.out = ndissim))#imaginary difference in a covariate
set.seed(42)
sim21 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),.1),1,2),nsp)
sim22 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),.5),1,2),nsp)
sim23 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),1),1,2),nsp)
sim24 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),2),1,2),nsp)
sim25 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),5),1,2),nsp)

sim26 <- simulate_dism_data(form,dat, matrix(c(qlogis(.1),1),1,2),nsp)
sim27 <- simulate_dism_data(form,dat, matrix(c(qlogis(.25),1),1,2),nsp)
sim28 <- simulate_dism_data(form,dat, matrix(c(qlogis(.5),1),1,2),nsp)
sim29 <- simulate_dism_data(form,dat, matrix(c(qlogis(.65),1),1,2),nsp)
sim30 <- simulate_dism_data(form,dat, matrix(c(qlogis(.8),1),1,2),nsp)

library(doParallel)#library(doMC)
library(foreach)
library(iterators)

# A two core cluster 
registerDoParallel(22)

models <- 
  list(
    a=list(data = sim1,nsites=25),
    b=list(data = sim2,nsites=25),
    c=list(data = sim3,nsites=25),
    d=list(data = sim4,nsites=25),
    e=list(data = sim5,nsites=25),
    f=list(data = sim6,nsites=25),
    g=list(data = sim7,nsites=25),
    h=list(data = sim8,nsites=25),
    i=list(data = sim9,nsites=25),
    j=list(data = sim10,nsites=25),
    aa=list(data = sim11,nsites=50),
    bb=list(data = sim12,nsites=50),
    cc=list(data = sim13,nsites=50),
    dd=list(data = sim14,nsites=50),
    ee=list(data = sim15,nsites=50),
    ff=list(data = sim16,nsites=50),
    gg=list(data = sim17,nsites=50),
    hh=list(data = sim18,nsites=50),
    ii=list(data = sim19,nsites=50),
    jj=list(data = sim20,nsites=50),
    aaa=list(data = sim21,nsites=100),
    bbb=list(data = sim22,nsites=100),
    ccc=list(data = sim23,nsites=100),
    ddd=list(data = sim24,nsites=100),
    eee=list(data = sim25,nsites=100),
    fff=list(data = sim26,nsites=100),
    ggg=list(data = sim27,nsites=100),
    hhh=list(data = sim28,nsites=100),
    iii=list(data = sim29,nsites=100),
    jjj=list(data = sim30,nsites=100)
    )

results <-
  foreach(i = iter(models),.packages='bbgdm') %dopar% {
    sim_bbgdm(X=i$data$X,y=i$data$y,Nsite = i$nsites,nboots = 100)
  }

save(results,file = 'sim_results_all.Rdata')
stopImplicitCluster()

for(i in 1:length(results)){
  plotResponseSim(results[[1]],plotdim = c(1,1))
  plotResponseSim(results[[3]],plotdim = c(1,1))



