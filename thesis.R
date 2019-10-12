m.files <- list.files(path = "lifetables", full.names = FALSE)

#prepare the data by making all numbers numeric 

#install.packages("HMDHFDplus")
library(HMDHFDplus)
forceType <- function(z){
  z <- z[-1,]
  colnames(z) <- as.character(unlist(z[1,]))
  z = z[-1, 1:10]
  z$Age <- age2int(z$Age)
  z$Year <- as.integer(as.character(z$Year))
  for(i in 3: 10){
    z[,i] <- as.numeric(as.character(z[,i]))
  }
  z
}
fdata <- forceType(read.delim(m.files[1], header = FALSE, sep = "", dec = "."))
mdata <- forceType(read.delim(m.files[2], header = FALSE, sep = "", dec = "."))

#merge male and female to obtain unisex data
library(plyr)
udata <- join(fdata, mdata, by = c("Year", "Age"))
udata[,6] <- udata[,6] + udata[,14]
udata <- udata[,c(1:2, 6)]

#this function creates a cohort for a given dataset 'z', with people bornin 'year', 
#of initial age 'yearMin' and maximum age 'yearMax'
cohort <- function(z, year, yearMin, yearMax){
  z <- z[ which(z$Year - z$Age == year), ]
  z <- z[ which(z$Age > yearMin & z$Age < yearMax), ]
  z
}
fcohort <- cohort(fdata,1915, 65, 110)
mcohort <- cohort(mdata,1915, 65, 110)
ucohort <- cohort(udata,1915, 65, 110)

#1900, 65, 110
#1885, 65, 110
#1915, 65, 110

#below the necessary functions are defined to compute the survival functions for both simple OU and jump OU processes
alpha.OU <- function(p,t){
  mu <- p[1]
  sigma <- abs(p[2])
  (sigma^2/(2*mu^2))*t - (sigma^2/mu^3)*exp(mu*t) + (sigma^2/(4*mu^3))*exp(2*mu*t) + (3*sigma^2)/(4*mu^3) 
}


alpha.OU.jump <- function(p,t){
  mu <- p[1]
  sigma <- p[2]
  eta <- abs(p[3])
  nu <- abs(p[4])
  
  (sigma^2/(2*mu^2))*t - (sigma^2/mu^3)*exp(mu*t) + (sigma^2/(4*mu^3))*exp(2*mu*t) + (3*sigma^2)/(4*mu^3) + ((eta*nu)/(mu - nu))*t - (eta/(mu - nu)) * log(1 - (nu/mu) + (nu/mu)*exp(mu*t))
}


beta <- function(p, t){
  a <- p[1]
  (1/a) * (1-exp(a*t))
}

survival.OU <- function(p,t,z){
  lambda <- -log(1- z$qx[1])
  exp(alpha.OU(p,t) + beta(p, t)*lambda)
}

survival.OU.jump <- function(p,t,z){
  lambda <- -log(1- z$qx[1])
  exp(alpha.OU.jump(p,t) + beta(p, t)*lambda)
}


#computes the mean square error of the survival function
MSE.OU <- function(p, z){
  difference <- rep(NA, nrow(z) - 1)
  tpx <- rep(NA, nrow(z) - 1)
  for(i in 1: nrow(z) - 1){
  tpx[i] <- z$lx[i+1]/z$lx[1]
  diff <- (tpx[i] - survival.OU(p,i,z))^2
  difference[i] <- diff
  }
  sum(difference)/length(difference)
}

MSE.OU.jump <- function(p,z){
  difference <- rep(NA, nrow(z) - 1)
  tpx <- rep(NA, nrow(z) - 1)
  for(i in 1: nrow(z) - 1){
    tpx[i] <- z$lx[i+1]/z$lx[1]
    diff <- (tpx[i] -survival.OU.jump(p,i, z))^2
    difference[i] <- diff
  }
  sum(difference)/length(difference)
}

#calibration of the parameters using the optim function and minimizing the MSE w.r.t the parameters to be estimated 
calibration.OU <- function(p.OU, z){
  p.OU <- c(0.01, 0.001)
  m.MSE.OU <- optim(p.OU, MSE.OU, z = z, hessian = TRUE)
  par.OU <- c(m.MSE.OU$par[1], abs(m.MSE.OU$par[2]))
  return(par.OU)
}

calibration.OU.j <- function(p.OU.j, z){
  p.OU.j <- c(0.01, 0.0001, 0.01, 0.003)
  m.MSE.OU.jump <- optim(p=p.OU.j, MSE.OU.jump,z= z)
  par.OU.jump <- c( m.MSE.OU.jump$par[1], abs(m.MSE.OU.jump$par[2]), abs(m.MSE.OU.jump$par[3]), abs(m.MSE.OU.jump$par[4]))
  return(par.OU.jump)
}

#observed survival probabilities
tpx <- function(z){
  tpx <- rep(NA, nrow(z) - 1)
  for(i in 1: nrow(z) - 1){
    tpx[i] <- z$lx[i+1]/z$lx[1]
  }
  return(tpx)
}


#computes estimated survival probabilities given the calibrated parameters
estimated.density.OU <- function(z, t){
 survival.OU(calibration.OU(p.OU,z), t, z)
}

estimated.density.OU.jump <- function(z, t){
  survival.OU.jump(calibration.OU.j(z = z), t, z)
}

#Graphs the observed and estimated survival probabilities for any cohort 'z'
library(ggplot2)
graphs <- function(z){
  obser <- data.frame(cbind(tpx(z), 1:(nrow(z) - 1)))
  x.grid <- seq(0,to=nrow(z) - 1,by=0.01)
  dd <- data.frame(x.grid,x2 = estimated.density.OU(z, t = x.grid), x3 = estimated.density.OU.jump(z, t = x.grid))
  ggplot(z) + 
  geom_line(data=dd,aes(x=x.grid,y=x2, colour = "OU")) + 
  geom_line(data=dd,aes(x=x.grid,y=x3, colour = "OU Jump")) +
  geom_line(data=obser, aes(x = X2, y = X1, colour = "Observed")) + 
  ylab("Survival probability") + xlab("Remaining lifetime") +
  scale_colour_manual(values=c("red","blue", "green"))
}

#probability of intensity becoming negative
prob <- function(z, t){
  mu <- calibration.OU.j(z = z)[1]
  sigma <- calibration.OU.j(z = z)[2]
  pnorm(log(1-z$qx[1])*exp(mu*t)/(sigma * sqrt((exp(2* mu*t )-1)/(2*mu))))
}

#make a table of estimates of the OU process for 3 cohorts
mestimates <- matrix(ncol = 3, nrow = 4)
colnames(mestimates) <- c("1885", "1900", "1915")
rownames(mestimates) <- c("mu", "sigma", "MSE", "lambda_0")
festimates <- matrix(ncol = 3, nrow = 4)
colnames(festimates) <- c("1885", "1900", "1915")
rownames(festimates) <- c("mu", "sigma", "MSE", "lambda_0")
year <- c(1885, 1900, 1915)
for(i in 1:3){
  mestimates[,i ] <- c(calibration.OU(z = cohort(mdata,year[i], 65, 110)), 
                    MSE.OU(calibration.OU(z = cohort(mdata,year[i], 65, 110)), z = cohort(mdata,year[i], 65, 110)), 
                    -log(1-cohort(mdata, year[i], 65, 110)$qx[1]))
  festimates[,i ] <- c(calibration.OU(z = cohort(fdata,year[i], 65, 110)), 
                       MSE.OU(calibration.OU(z = cohort(fdata,year[i], 65, 110)), z = cohort(fdata,year[i], 65, 110)), 
                     -log(1-cohort(fdata, year[i], 65, 110)$qx[1]))
}

#now for the OU process with jump
mestimates.j <- matrix(ncol = 3, nrow = 6)
colnames(mestimates.j) <- c("1885", "1900", "1915")
rownames(mestimates.j) <- c("mu", "sigma", "eta", "nu", "MSE", "lambda_0")
festimates.j <- matrix(ncol = 3, nrow = 6)
colnames(festimates.j) <- c("1885", "1900", "1915")
rownames(festimates.j) <- c("mu", "sigma", "eta", "nu", "MSE", "lambda_0")
year <- c(1885, 1900, 1915)
for(i in 1:3){
  mestimates.j[,i ] <- c(calibration.OU.j(z = cohort(mdata,year[i], 65, 110)), 
                       MSE.OU.jump(calibration.OU.j(z = cohort(mdata,year[i], 65, 110)), z = cohort(mdata,year[i], 65, 110)), 
                       -log(1-cohort(mdata, year[i], 65, 110)$qx[1]))
  festimates.j[,i ] <- c(calibration.OU.j(z = cohort(fdata,year[i], 65, 110)), 
                       MSE.OU.jump(calibration.OU.j(z = cohort(fdata,year[i], 65, 110)), z = cohort(fdata,year[i], 65, 110)), 
                       -log(1-cohort(fdata, year[i], 65, 110)$qx[1]))
}

#defines the shocked survival probability density
shockSurvival <- function(t, e, z){
  estimated.density.OU.jump(z=z, t =t)^(1+e)
}

#computes the estimated shocked mortality as a function of the shock e and using the graph, 
#we can visualize the severity of the shock in order for the shocked mortality to be 15 percent higher after 1 year
par.OU.j = calibration.OU.j(z = fcohort)
x <- seq(0, to = 0.2, by = 0.001)
est.mortality <- 1 - shockSurvival(x, t = 1, z = fcohort)
SS <- data.frame(x, est.mortality, qs = (1+0.15)*(1-survival.OU.jump(par.OU.j, t = 1, z = fcohort)))
ggplot(SS) + 
  geom_line(data = SS, aes(x = x, y = est.mortality))+
  geom_line(data = SS, aes(x = x, y = qs )) #the intercept is somewhere around e = 0.15 for female data

#function to define the shock parameter epsilon
epsilon <- function(z){
  shock <- function(epsilon, s, z){ 
    ( 1- estimated.density.OU.jump(z=z, t =1)^(1+epsilon) - (1 + s)*(1 - estimated.density.OU.jump(z=z, t=1)))^2 #minimize (1-Sx(1)^(1+e) - (1+s)Qx)^2 wrt e
  }
  return(optimize(shock, interval= c(0,1), s = 0.15, z= z)$minimum)
}

#creates a graph of the shocked and regular survival function
graphShock <- function(z){
  e <- epsilon(z)
  t <- seq(0, to = nrow(z) - 1, by = 0.01)
  est.SS <- shockSurvival(t=t, e = e, z=z)
  datashock <- data.frame(t, x1 = est.SS, Sxt =estimated.density.OU.jump(z = z, t = t))
  ggplot(datashock) + 
    geom_line(aes(x=t,y=Sxt, colour = "without shock")) + 
    geom_line(aes(x=t,y=x1, colour = "shock")) 
}

#n year term life insurance
nya <- function(n, r=0.02, z){
  units <- rep(NA, n)
  par <- calibration.OU.j(z= z)
  for(i in 1:n){
    units[i] <- 1/(1+r)^(i+1)*(survival.OU.jump(z=z, t = i-1, p  = par) - survival.OU.jump(z=z, t= i, p = par))
  }
  return(sum(units))
}

#shocked n-year term life insurance
nyas <- function(n, r=0.02, z){
  units <- rep(NA, n)
  par <- calibration.OU.j(z= z)
  e <- epsilon(z)
  for(i in 1:n){
    units[i] <- 1/(1+r)^(i+1)*(survival.OU.jump(z=z, t = i-1, p  = par)^(1+e) - survival.OU.jump(z=z, t= i, p = par)^(1+e))
  }
  return(sum(units))
}

#SCR n year term life insurance
SCRnya <- function(n, z){
  nyas(n, z=z) - nya(n, z=z)
}

#n year endowment  
nye <- function(n, r=0.02, z){
  units <- rep(NA, n)
  par <- calibration.OU.j(z= z)
  for(i in 1:n){
    units[i] <- 1/(1+r)^(i+1)*(survival.OU.jump(par, z = z, t =i-1) - survival.OU.jump(p = par, z=z, t=i)) 
  }
  return(sum(units) + 1/(1+r)^n*(survival.OU.jump(par, z=z, t= n)))
}

#shocked n year  endowment 
nyes <- function(n, r=0.02, z){
  units <- rep(NA, n)
  par <- calibration.OU.j(z= z)
  e <- epsilon(z)
  for(i in 1:n){
    units[i] <- 1/(1+r)^(i+1)*(survival.OU.jump(z=z, t = i-1, p  = par)^(1+e) - survival.OU.jump(z=z, t= i, p = par)^(1+e)) 
  }
  return(sum(units) + 1/(1+r)^n*(survival.OU.jump(z=z, t = n, p  = par)^(1+e)))
}

#SCR n year endowment
SCRnye <- function(n, z){
  nyes(n, z=z) - nye(n, z=z)
}

#to check calculations, by definition, we must have SCRnya - SCRnye = SCRend
SCRnya(2, fcohort) - SCRnye(2, fcohort)
par <- calibration.OU.j(z = fcohort)
1/(1+0.02)^2*(survival.OU.jump(z=fcohort, t = 2, p  = par)) - 1/(1+0.02)^2*(survival.OU.jump(z=fcohort, t = 2, p  = par))^(1+epsilon(fcohort)) 
#since they are indeed equal, the calculations must have gone right

#alpha for unisex data and survival function for unisex data, using estimates from separate gender processes
alpha2 <- function(p, t){
  mu1 <- p[1] #female
  s1 <- p[2]
  mu2 <- p[3] #male
  s2 <- p[4]
  g <- p[5] #weight
  r <- p[6] #correlation
  
  ((s2)^2 * g^2 / (4*mu2^3 ) )* ((exp(mu2*t) - 2)^2 + 2*mu2*t - 1) + 
  ((s1)^2 * (g-1)^2 / (4*mu1^3 ) )* ((exp(mu1*t) - 2)^2 + 2*mu1*t - 1) - 
  (r * s1 * s2 * g * (g-1)) / (mu1^2 * mu2^2* (mu1 + mu2)) * 
  (mu2^2*(1-exp(mu1*t)) + mu1^2*(1-exp(mu2 * t)) + 
  mu1 * mu2 * ((1 - exp(mu2 * t))*(1-exp(mu1*t)) + (mu1 + mu2)*t))
}

survival.OU.uni <- function(p, t, z1, z2){
  lambda.f <- - log(1- z1$qx[1])
  lambda.m <- - log(1- z2$qx[1])
  exp( alpha2(p, t) + (1-p[5]) * beta(p = p[1], t) * lambda.f + p[5]*beta(p = p[3], t)* lambda.m)
}

graphSurvivalFunctions.OU <- function(z1, z2){
  obserfemale <- data.frame(f = cbind(tpx(z1), 1:(nrow(z1) - 1)))
  obsermale <- data.frame(m = cbind(tpx(z2), 1:(nrow(z2) - 1)))
  obserunisex <- data.frame(u = cbind(tpx(ucohort), 1:(nrow(ucohort) - 1)))
  obser <- cbind( obserfemale, obsermale, obserunisex)
  x.grid <- seq(0,to=nrow(z1) - 1,by=0.01)
  dd <- data.frame(x.grid,female = estimated.density.OU(z1, t = x.grid), male = estimated.density.OU(z2, t = x.grid), 
                   unisex = survival.OU.uni(p = c(calibration.OU(round(maxlik(z1), digits = 4), z1), calibration.OU(round(maxlik(z2), digits = 4), z2), 0.5, 0.95), t = x.grid, z1 = z1, z2 = z2))
  
  library(ggplot2)
  ggplot(dd) + 
    geom_line(data=dd,aes(x=x.grid,y=female, colour = "female")) + 
    geom_line(data=dd,aes(x=x.grid,y=male, colour = "male")) +
    #geom_line(data=dd,aes(x=x.grid,y=unisex, colour = "unisex")) +
    geom_line(data=obser, aes(x = f.2, y = f.1, colour = "Observed female")) + 
    geom_line(data=obser, aes(x = m.2, y = m.1, colour = "Observed male")) + 
    #geom_line(data=obser, aes(x = u.2, y = u.1, colour = "Observed unisex")) + 
    ylab("Survival probability") + xlab("Remaining lifetime") +
    scale_colour_manual(values=c("red","blue", "black", "yellow")) + 
    theme(legend.position="bottom")
}

#alpha for unisex data with jump component 
alpha3 <- function(p, t){
  mu1 <- p[1] #female
  s1 <- p[2]
  mu2 <- p[5] #male
  s2 <- p[6]
  g <- p[9] #weight
  r <- p[10] #correlation
  eta1 <- p[3] #jump arrival
  eta2 <- p[7]
  nu1 <- p[4] #jump size
  nu2 <- p[8]

  g^2*((s2^2/(2*mu2^2))*t - (s2^2/mu2^3)*exp(mu2*t) + (s2^2/(4*mu2^3))*exp(2*mu2*t) + (3*s2^2)/(4*mu2^3)) + 
    (1-g)^2*((s1^2/(2*mu1^2))*t - (s1^2/mu1^3)*exp(mu1*t) + (s1^2/(4*mu1^3))*exp(2*mu1*t) + (3*s1^2)/(4*mu1^3)) - 
    (r * s1 * s2 * g * (g-1)) / (mu1^2 * mu2^2* (mu1 + mu2)) * 
    (mu2^2*(1-exp(mu1*t)) + mu1^2*(1-exp(mu2 * t)) + 
       mu1 * mu2 * ((1 - exp(mu2 * t))*(1-exp(mu1*t)) + (mu1 + mu2)*t)) +
    eta1/( mu1 - (1-g) *nu1)*((1-g)*nu1*t - log(1 - nu1* (1-g)/mu1 *(1- exp(mu1*t)))) +
    eta2/( mu2 - (g) *nu2)*((g)*nu2*t - log(1 - nu2* (g)/mu2 *(1- exp(mu2*t))))
}


survival.OU.uni.jump <- function(p, t, z1, z2){
  lambda.f <- - log(1- z1$qx[1])
  lambda.m <- - log(1- z2$qx[1])
  exp( alpha3(p, t) + (1-p[9]) * beta(p = p[1], t) * lambda.f + p[9]*beta(p = p[5], t)* lambda.m)
}


graphSurvivalFunctions.OUj <- function(z1, z2){
  obserfemale <- data.frame(f = cbind(tpx(z1), 1:(nrow(z1) - 1)))
  obsermale <- data.frame(m = cbind(tpx(z2), 1:(nrow(z2) - 1)))
  obserunisex <- data.frame(u = cbind(tpx(ucohort), 1:(nrow(ucohort) - 1)))
  obser <- cbind( obserfemale, obsermale, obserunisex)
  x.grid <- seq(0,to=nrow(z1) - 1,by=0.01)
  dd <- data.frame(x.grid,female = estimated.density.OU.jump(z1, t = x.grid), male = estimated.density.OU.jump(z2, t = x.grid), 
                   unisex.jump = survival.OU.uni.jump(p = c(calibration.OU.j(z = z1), calibration.OU.j(z = z2), 0.5, 0.95), t = x.grid, z1 = z1, z2 = z2)) 
  
  library(ggplot2)
  ggplot(dd) + 
    geom_line(data=dd,aes(x=x.grid,y=female, colour = "female")) + 
    geom_line(data=dd,aes(x=x.grid,y=male, colour = "male")) +
    geom_line(data=dd, aes(x = x.grid, y = unisex.jump, colour = "unisex jump"))+
    geom_line(data=obser, aes(x = f.2, y = f.1, colour = "Observed female")) + 
    geom_line(data=obser, aes(x = m.2, y = m.1, colour = "Observed male")) + 
    geom_line(data=obser, aes(x = u.2, y = u.1, colour = "Observed unisex")) + 
    ylab("Survival probability") + xlab("Remaining lifetime") +
    scale_colour_manual(values=c("red","blue", "black", "yellow", "green", "purple"))
}



differenceNYE <- function(n, z1, z2, w){ #calculates the relative difference for uniform SCR and weighted SCR given the data, duration and proportion
  
  #Calculates the n-year term endowment for unisex data 
  nye.uni <- function(n , r=0.02, z1, z2, g ){
    par <- c(calibration.OU.j(z = z1), calibration.OU.j(z = z2), g, 0.95)
    units <- rep(NA, n)
    for(i in 1:n){
      units[i] <- 1/(1+r)^(i+1)* (survival.OU.uni.jump(p = par, t = i-1, z1 = z1, z2 = z2)-survival.OU.uni.jump(p = par, t = i, z1 = z1, z2 = z2))
    }
  return(sum(units) + 1/(1+r)^n*survival.OU.uni.jump(p = par, t = n, z1 = z1, z2 = z2)) }

  #Calculates shocked n-year term endowment for unisex data 
  nyes.uni <- function(n, r=0.02, z1, z2, e, g){
    par <- c(calibration.OU.j(z = z1), calibration.OU.j(z = z2), g, 0.95)
    units <- rep(NA, n)
    for(i in 1:n){
      units[i] <- 1/(1+r)^(i+1)* (survival.OU.uni.jump(p = par, t = i-1, z1 = z1, z2 = z2)^(1+e)  -
      survival.OU.uni.jump(p = par, t = i, z1 = z1, z2 = z2)^(1+e))
    }
  return(sum(units) + 1/(1+r)^n*survival.OU.uni.jump(p = par, t = n, z1 = z1, z2 = z2)^( 1+e))
  }
  
  #calculates the SCR of the n-year term endowment 
  SCR.nye.u <- function(n, z1, z2, g, e){
    nyes.uni(n=n, z1 = z1, z2 = z2, e = e, g=g) - nye.uni(n =n, z1 = z1, z2 = z2, g=g)
  }
  
  #given the proportion of males in the dataset, calculate g such that premium(unisex) = w*male + (1-w)*female
  #w <- z2$lx[1]/(z2$lx[1] + z1$lx[1]) proportion of males in dataset
  wNYE <- function(z1, z2, n, g){ #z1 female, z2 male
    (nye.uni(n, z1 = z1, z2 = z2, g = g) - (w * nye(n, z = z2) + (1-w)*nye(n, z = z1)))^2
  }
  g <- optimize(wNYE, interval = c(0,1), z1 = z1, z2 = z2, n = n)$minimum #g such that unisex premium satisfies actuarial fairness principle
  
  
  shock <- function(epsilon, s, z1, z2, g){
    par <- c(calibration.OU.j(z = z1), calibration.OU.j( z= z2), g, 0.95)
    return(( 1- survival.OU.uni.jump(par, t =1, z1, z2)^(1+epsilon) - (1 + s)*(1 - survival.OU.uni.jump(par, t=1, z1, z2)))^2) #minimize (1-Sx(1)^(1+e) - (1+s)Qx)^2 wrt e
  }
  e.u <- optimize(shock, interval= c(0,1), s = 0.15, z1 = z1, z2 = z2, g=g)$minimum
  
  return(c(((SCR.nye.u(n= n, z1 = z1, z2 = z2, g = g, e = e.u) - (w*SCRnye(n =n, z = z2) + (1-w)*(SCRnye(n=n, z =z1))))/(w*SCRnye(n =n, z = z2) + (1-w)*(SCRnye(n=n, z =z1)))), 
         g, e.u))
}

differenceNYA <- function(n, z1, z2, w){ #calculates the relative difference for uniform SCR and weighted SCR given the data, duration and proportion
  
  #n year term life insurance
  nya.uni <- function(n, r=0.02, z1, z2, g){
    units <- rep(NA, n)
    par <- c(calibration.OU.j(z = z1), calibration.OU.j( z= z2), g, 0.95)
    for(i in 1:n){
      units[i] <- 1/(1+r)^(i+1)*(survival.OU.uni.jump(p = par, t = i-1, z1 = z1, z2 = z2) - survival.OU.uni.jump(p = par, t = i, z1 = z1, z2 = z2))
    }
    return(sum(units))
  }
  
  #shocked n-year term life insurance
  nyas.uni <- function(n, r=0.02, z1, z2, g, e){
    units <- rep(NA, n)
    par <- c(calibration.OU.j(z = z1), calibration.OU.j( z= z2), g, 0.95)
    for(i in 1:n){
      units[i] <- 1/(1+r)^(i+1)*(survival.OU.uni.jump(p = par, t = i-1, z1 = z1, z2 = z2)^(1+e) - survival.OU.uni.jump(p = par, t = i, z1 = z1, z2 = z2)^(1+e))
    }
    return(sum(units))
  }
  
  #SCR n year term life insurance
  SCR.nya.u <- function(n, z1, z2, g, e){
    nyas.uni(n=n, z1 = z1, z2 = z2, e = e, g=g) - nya.uni(n =n, z1 = z1, z2 = z2, g=g)
  }
  
  #given the proportion of males in the dataset, calculate g such that premium(unisex) = w*male + (1-w)*female
  #w <- z2$lx[1]/(z2$lx[1] + z1$lx[1]) proportion of males in dataset
  wNYE <- function(z1, z2, n, g){ #z1 female, z2 male
    (nya.uni(n, z1 = z1, z2 = z2, g = g) - (w * nya(n, z = z2) + (1-w)*nya(n, z = z1)))^2
  }
  g <- optimize(wNYE, interval = c(0,1), z1 = z1, z2 = z2, n = n)$minimum #g such that unisex premium satisfies actuarial fairness principle
  
  
  shock <- function(epsilon, s, z1, z2, g){
    par <- c(calibration.OU.j(z = z1), calibration.OU.j( z= z2), g, 0.95)
    return(( 1- survival.OU.uni.jump(par, t =1, z1, z2)^(1+epsilon) - (1 + s)*(1 - survival.OU.uni.jump(par, t=1, z1, z2)))^2) #minimize (1-Sx(1)^(1+e) - (1+s)Qx)^2 wrt e
  }
  e.u <- optimize(shock, interval= c(0,1), s = 0.15, z1 = z1, z2 = z2, g=g)$minimum
  
  return(c((SCR.nya.u(n= n, z1 = z1, z2 = z2, g = g, e = e.u) - (w*SCRnya(n =n, z = z2) + (1-w)*(SCRnya(n=n, z =z1))))
         /(w*SCRnya(n =n, z = z2) + (1-w)*(SCRnya(n=n, z =z1))), g, e.u))
}


#inputW <- c(0, 0.1, 0.2, 0.3, 0.4,  0.5, 0.6, 0.7, 0.8, 0.9, 1)
#outputNYA <- matrix(nrow = length(inputW), ncol = 3)
#colnames(outputNYA) <- c("SCR", "g", "e")
#rownames(outputNYA) <- inputW
#for(i in 1:length(inputW)){
   #outputNYA[i,] <- differenceNYA(n = 20, z1 = fcohort, z2 = mcohort, w = inputW[i])
   #print(outputNYA[i,])
#}

#outputNYE <- matrix(nrow = length(inputW), ncol = 3)
#colnames(outputNYE) <- c("SCR", "g", "e")
#rownames(outputNYE) <- inputW
#for(i in 1:length(inputW)){
  #outputNYE[i,] <- differenceNYE(n = 20, z1 = fcohort, z2 = mcohort, w = inputW[i])
  #print(outputNYE[i,])
#}


