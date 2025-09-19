#1
S1 <- function(t){
	return(exp(-t))
}
#2
S2 <- function(t){
	return(exp(-t^0.5))
}
#3
S3 <- function(t){
	return(exp(-t^2))
}
#4
S4 <- function(t){
	return(exp(-2*t))
}
#たぶん教科書のスケールパラメータlambdaの逆数がRのスケールパラメータ
#5
S5 <- function(t){
	return(exp(-0.5*t))
}
#6
h6 <- function(x){
	if (0<x && x<=0.5){
		return(1)
	}else if(0.5<x && x<=1){
		return(2)
	}else return(0)
}
S6 <- function(t){
	if (0<t && t<=0.5){
		return(exp(-t))
	}else if(0.5<t && t<=1){
		return(exp(0.5-2*t))
	}else return(0)
}
S6_1 <- function(t){
	return(exp(-t))
}
S6_2 <- function(t){
	return(exp(0.5-2*t))
}
dS6 <-function(n){
	Surv <- seq(0,1,length=n)
	Survival <- numeric(n)
	Survival[0<=Surv & Surv <= 0.5] <- exp(-Surv[0<=Survival & Surv <= 0.5])
	Survival[0.5<Surv & Surv <= 1] <- exp(0.5-2*Surv[0.5<Surv & Surv <= 1])
	return(invisible(Survival))
}
#1e-6でいいのか問題
rf6 <- function(n){
	randomunifs <- 1 - runif(n,1e-6,1)
	randoms <- numeric(n)
	randoms[randomunifs>=exp(-0.5) & randomunifs < 1] <- -log(randomunifs[randomunifs>=exp(-0.5) & randomunifs < 1])
	randoms[randomunifs>=0.25 & randomunifs < exp(-0.5)] <- 0.25-0.5*log(randomunifs[randomunifs>=0.25 & randomunifs < exp(-0.5)])
	randoms[randomunifs<0.25] <- Inf
	return(randoms)
}
#7
h7 <- function(x){
	if (0<x && x<=0.5){
		return(2)
	}else if(0.5<x && x<=1){
		return(1)
	}else return(0)
}
S7 <- function(t){
	if (0<t && t<=0.5){
		return(exp(-2*t))
	}else if(0.5<t && t<=1){
		return(exp(-t-0.5))
	}else return(0)
}
S7_1 <- function(t){
	return(exp(-2*t))
}
S7_2 <- function(t){
	return(exp(-t-0.5))
}
dS7 <-function(n){
	Surv <- seq(0,1,length=n)
	Survival <- numeric(n)
	Survival[0<=Surv & Surv <= 0.5] <- exp(-2*Surv[0<=Survival & Surv <= 0.5])
	Survival[0.5<Surv & Surv <= 1] <- exp(-0.5-Surv[0.5<Surv & Surv <= 1])
	return(invisible(Survival))
}
rf7 <- function(n){
	randomunifs <- 1 - runif(n,1e-6,1)
	randoms <- numeric(n)
	randoms[randomunifs>=exp(-1) & randomunifs < 1] <- -0.5*log(randomunifs[randomunifs>=exp(-1) & randomunifs < 1])
	randoms[randomunifs>=exp(-1.5) & randomunifs < exp(-1)] <- -0.5-log(randomunifs[randomunifs>=exp(-1.5) & randomunifs < exp(-1)])
	randoms[randomunifs<exp(-1.5)] <- Inf
	return(randoms)
}
#8
h8 <- function(x){
	if (0<x && x<=0.5){
		return(2*t+1)
	}else if(0.5<x && x<=1){
		return(3-2*t)
	}else return(0)
}
S8 <- function(t){
	if (0<t && t<=0.5){
		return(exp(-t^2-t))
	}else if(0.5<t && t<=1){
		return(exp(t^2-3*t+0.5))
	}else return(0)
}
S8_1 <- function(t){
	return(exp(-t^2-t))
}
S8_2 <- function(t){
	return(exp(t^2-3*t+0.5))
}
dS8 <-function(n){
	Surv <- seq(0,1,length=n)
	Survival <- numeric(n)
	Survival[0<=Surv & Surv <= 0.5] <- exp(-Surv[0<=Survival & Surv <= 0.5]^2-Surv[0<=Survival & Surv <= 0.5])
	Survival[0.5<Surv & Surv <= 1] <- exp(Surv[0.5<Surv & Surv <= 1]^2-3*Surv[0.5<Surv & Surv <= 1]+0.5)
	return(invisible(Survival))
}
rf8 <- function(n){
	randomunifs <- 1 - runif(n,1e-6,1)
	randoms <- numeric(n)
	randoms[randomunifs>=exp(-0.75) & randomunifs < 1] <- -0.5 + sqrt(0.25-log(randomunifs[randomunifs>=exp(-0.75) & randomunifs < 1]))
	randoms[randomunifs>=exp(-1.5) & randomunifs < exp(-0.75)] <- 1.5-sqrt(1.75+log(randomunifs[randomunifs>=exp(-1.5) & randomunifs < exp(-0.75)]))
	randoms[randomunifs<exp(-1.5)] <- Inf
	return(randoms)
}
#9
h9 <- function(x){
	if (0<x && x<=0.5){
		return(-2*t+2)
	}else if(0.5<x && x<=1){
		return(2*t)
	}else return(0)
}
S9 <- function(t){
	if (0<t && t<=0.5){
		return(exp(t^2-2*t))
	}else if(0.5<t && t<=1){
		return(exp(-t^2-0.5))
	}else return(0)
}
S9_1 <- function(t){
	return(exp(t^2-2*t))
}
S9_2 <- function(t){
	return(exp(-t^2-0.5))
}
dS9 <-function(n){
	Surv <- seq(0,1,length=n)
	Survival <- numeric(n)
	Survival[0<=Surv & Surv <= 0.5] <- exp(Surv[0<=Survival & Surv <= 0.5]^2-2*Surv[0<=Survival & Surv <= 0.5])
	Survival[0.5<Surv & Surv <= 1] <- exp(-Surv[0.5<Surv & Surv <= 1]^2-0.5)
	return(invisible(Survival))
}
rf9 <- function(n){
	randomunifs <- 1 - runif(n,1e-6,1)
	randoms <- numeric(n)
	randoms[randomunifs>=exp(-0.75) & randomunifs < 1] <- 1 - sqrt(1+log(randomunifs[randomunifs>=exp(-0.75) & randomunifs < 1]))
	randoms[randomunifs>=exp(-1.5) & randomunifs < exp(-0.75)] <- sqrt(-0.5-log(randomunifs[randomunifs>=exp(-1.5) & randomunifs < exp(-0.75)]))
	randoms[randomunifs<exp(-1.5)] <- Inf
	return(randoms)
}
#10
h10 <- function(x){
	if (0<x && x<=0.5){
		return(1)
	}else if(0.5<x && x<=1){
		return(2*t)
	}else return(0)
}
S10 <- function(t){
	if (0<t && t<=0.5){
		return(exp(-t))
	}else if(0.5<t && t<=1){
		return(exp(-t^2-0.25))
	}else return(0)
}
S10_1 <- function(t){
	return(exp(-t))
}
S10_2 <- function(t){
	return(exp(-t^2-0.25))
}
dS10 <-function(n){
	Surv <- seq(0,1,length=n)
	Survival <- numeric(n)
	Survival[0<=Surv & Surv <= 0.5] <- exp(-Surv[0<=Survival & Surv <= 0.5])
	Survival[0.5<Surv & Surv <= 1] <- exp(-Surv[0.5<Surv & Surv <= 1]^2-0.25)
	return(invisible(Survival))
}
rf10 <- function(n){
	randomunifs <- 1 - runif(n,1e-6,1)
	randoms <- numeric(n)
	randoms[randomunifs>=exp(-0.5) & randomunifs < 1] <- -log(randomunifs[randomunifs>=exp(-0.5) & randomunifs < 1])
	randoms[randomunifs>=exp(-1.25) & randomunifs < exp(-0.5)] <- sqrt(-0.25-log(randomunifs[randomunifs>=exp(-1.25) & randomunifs < exp(-0.5)]))
	randoms[randomunifs<exp(-1.25)] <- Inf
	return(randoms)
}
#11
h11 <- function(x){
	if (0<x && x<=0.5){
		return(2)
	}else if(0.5<x && x<=1){
		return(3-2*x)
	}else return(0)
}
S11 <- function(t){
	if (0<t && t<=0.5){
		return(exp(-2*t))
	}else if(0.5<t && t<=1){
		return(exp(t^2-3*t+0.25))
	}else return(0)
}
S11_1 <- function(t){
	return(exp(-2*t))
}
S11_2 <- function(t){
	return(exp(t^2-3*t+0.25))
}
dS11 <-function(n){
	Surv <- seq(0,1,length=n)
	Survival <- numeric(n)
	Survival[0<=Surv & Surv <= 0.5] <- exp(-2*Surv[0<=Survival & Surv <= 0.5])
	Survival[0.5<Surv & Surv <= 1] <- exp(Surv[0.5<Surv & Surv <= 1]^2-3*Surv[0.5<Surv & Surv <= 1]+0.25)
	return(invisible(Survival))
}
rf11 <- function(n){
	randomunifs <- 1 - runif(n,1e-6,1)
	randoms <- numeric(n)
	randoms[randomunifs>=exp(-1) & randomunifs < 1] <- -0.5*log(randomunifs[randomunifs>=exp(-1) & randomunifs < 1])
	randoms[randomunifs>=exp(-1.75) & randomunifs < exp(-1)] <- 1.5-sqrt(2+log(randomunifs[randomunifs>=exp(-1.75) & randomunifs < exp(-1)]))
	randoms[randomunifs<exp(-1.75)] <- Inf
	return(randoms)
}
png("plot1.png", width = 800, height = 600)
curve(1-pweibull(x,1,1), 0, 1, type='l', ylim=c(0,1.0), ylab="S(x)", main="1")
dev.off()

