#1
T1 <- function(t){
  return(exp(-t))
}

#2
T2 <- function(t){
  return(exp(-t^0.5))
}

#3
T3 <- function(t){
  return(exp(-t^2))
}

#4
#[0,0.5]
T4_1 <- function(t){
  return(exp(-2*t))
}
#[0.5,1.0]
T4_2 <- function(t){
  return(exp(-t-0.5))
}
rT4 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-1) & randomunifs < 1] <- -0.5*log(randomunifs[randomunifs>=exp(-1) & randomunifs < 1])
  randoms[randomunifs>=exp(-1.5) & randomunifs < exp(-1)] <- -0.5-log(randomunifs[randomunifs>=exp(-1.5) & randomunifs < exp(-1)])
  randoms[randomunifs<exp(-1.5)] <- Inf
  return(randoms)
}
#5
#[0,0.5]
T5_1 <- function(t){
  return(exp(-t))
}
#[0.5,1.0]
T5_2 <- function(t){
  return(exp(-t^2-0.25))
}
rT5 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-0.5) & randomunifs < 1] <- -log(randomunifs[randomunifs>=exp(-0.5) & randomunifs < 1])
  randoms[randomunifs>=exp(-1.25) & randomunifs < exp(-0.5)] <- sqrt(-0.25-log(randomunifs[randomunifs>=exp(-1.25) & randomunifs < exp(-0.5)]))
  randoms[randomunifs<exp(-1.25)] <- Inf
  return(randoms)
}
#6
#Control
#T1
#Treatment
T6 <- function(t){
  return(exp(-2*t))
}

#7
#Control
#T2
#Treatment
T7 <- function(t){
  return(exp(-(2*t)^0.5))
}

#8
#Control
#T3
#Treatment
T8 <- function(t){
  return(exp(-(4/3*t)^2))
}

#9
#Control
#T4_1.T4_2
#Treatment
#[0,0.5]
T9_1 <- function(t){
  return(exp(-3*t))
}
#[0.5,1.0]
T9_2 <- function(t){
  return(exp(-1.5*t-0.75))
}
rT9 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-1.5) & randomunifs < 1] <- -(1/3)*log(randomunifs[randomunifs>=exp(-1.5) & randomunifs < 1])
  randoms[randomunifs>=exp(-2.25) & randomunifs < exp(-1.5)] <- -0.5-(2/3)*log(randomunifs[randomunifs>=exp(-2.25) & randomunifs < exp(-1.5)])
  randoms[randomunifs<exp(-2.25)] <- Inf
  return(randoms)
}
#10
#Control
#T5_1.T5_2
#Treatment
#[0,0.5]
T10_1 <- function(t){
  return(exp(-1.5*t))
}
#[0.5,1.0]
T10_2 <- function(t){
  return(exp(-1.5*t^2-3/8))
}
rT10 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-0.75) & randomunifs < 1] <- -(2/3)*log(randomunifs[randomunifs>=exp(-0.75) & randomunifs < 1])
  randoms[randomunifs>=exp(-(15/8)) & randomunifs < exp(-0.75)] <- sqrt(-0.25-(2/3)*log(randomunifs[randomunifs>=exp(-(15/8)) & randomunifs < exp(-0.75)]))
  randoms[randomunifs<exp(-(15/8))] <- Inf
  return(randoms)
}
#11
#Control
#[0,0.5]
T11_1 <- function(t){
  return(exp(-1.5*t^2-0.25*t))
}
#[0.5,1.0]
T11_2 <- function(t){
  return(exp(-0.5*t^2-1.25*t+0.25))
}
rT11 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-0.5) & randomunifs < 1] <- -(1/12)+sqrt(1/144-(2/3)*log(randomunifs[randomunifs>=exp(-0.5) & randomunifs < 1]))
  randoms[randomunifs>=exp(-1.5) & randomunifs < exp(-0.5)] <- -1.25+sqrt(33/16-2*log(randomunifs[randomunifs>=exp(-1.5) & randomunifs < exp(-0.5)]))
  randoms[randomunifs<exp(-1.5)] <- Inf
  return(randoms)
}
#Treatment
#[0,0.5]
T12_1 <- function(t){
  return(exp(-1.75*t))
}
#[0.5,1.0]
T12_2 <- function(t){
  return(exp(-0.5*t^2-1.25*t-1/8))
}
rT12 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-7/8) & randomunifs < 1] <- -(4/7)*log(randomunifs[randomunifs>=exp(-7/8) & randomunifs < 1])
  randoms[randomunifs>=exp(-15/8) & randomunifs < exp(-7/8)] <- -1.25+sqrt(21/16-2*log(randomunifs[randomunifs>=exp(-15/8) & randomunifs < exp(-7/8)]))
  randoms[randomunifs<exp(-15/8)] <- Inf
  return(randoms)
}
#12
#Control
#[0,0.2]
T13_1 <- function(t){
  return(exp(-2*t))
}
#[0.2,1.0]
T13_2 <- function(t){
  return(exp(t^2-2.4*t+1/25))
}
rT13 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-0.4) & randomunifs < 1] <- -0.5*log(randomunifs[randomunifs>=exp(-0.4) & randomunifs < 1])
  randoms[randomunifs>=exp(-34/25) & randomunifs < exp(-0.4)] <- 1.2-sqrt(7/5+log(randomunifs[randomunifs>=exp(-34/25) & randomunifs < exp(-0.4)]))
  randoms[randomunifs<exp(-34/25)] <- Inf
  return(randoms)
}
#Treatment
#[0,0.2]
#T13_1
#[0.2,1.0]
T14_2 <- function(t){
  return(exp(-2*t^2-(6/5)*t-2/25))
}
rT14 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-0.4) & randomunifs < 1] <- -0.5*log(randomunifs[randomunifs>=exp(-0.4) & randomunifs < 1])
  randoms[randomunifs>=exp(-82/25) & randomunifs < exp(-0.4)] <- -0.3+sqrt(1/20-0.5*log(randomunifs[randomunifs>=exp(-82/25) & randomunifs < exp(-0.4)]))
  randoms[randomunifs<exp(-82/25)] <- Inf
  return(randoms)
}
#13
#Control
T15 <- function(t){
  return(exp(-(3/4)*t^2-0.5*t))
}
rT15 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-1.25) & randomunifs < 1] <- -(1/3)+sqrt(1/9-(4/3)*log(randomunifs[randomunifs>=exp(-1.25) & randomunifs < 1]))
  randoms[randomunifs<exp(-1.25)] <- Inf
  return(randoms)
}
#Treatment
T16 <- function(t){
  return(exp((3/4)*t^2-2*t))
}
rT16 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-1.25) & randomunifs < 1] <- 4/3-sqrt(16/9+(4/3)*log(randomunifs[randomunifs>=exp(-1.25) & randomunifs < 1]))
  randoms[randomunifs<exp(-1.25)] <- Inf
  return(randoms)
}
#14
#Control
T17 <- function(t){
  return(exp(-0.5*t^2-0.5*t))
}
rT17 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-1) & randomunifs < 1] <- -0.5+sqrt(0.25-2*log(randomunifs[randomunifs>=exp(-1) & randomunifs < 1]))
  randoms[randomunifs<exp(-1)] <- Inf
  return(randoms)
}
#Treatment
#[0,0.5]
T18_1 <- function(t){
  return(exp(-1.5*t))
}
#[0.5,1.0]
T18_2 <- function(t){
  return(exp(-0.5*t-0.5))
}
rT18 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-0.75) & randomunifs < 1] <- -(2/3)*log(randomunifs[randomunifs>=exp(-0.75) & randomunifs < 1])
  randoms[randomunifs>=exp(-1) & randomunifs < exp(-0.75)] <- -1-2*log(randomunifs[randomunifs>=exp(-1) & randomunifs < exp(-0.75)])
  randoms[randomunifs<exp(-1)] <- Inf
  return(randoms)
}
#15
#Control
#[0,0.2]
T19_1 <- function(t){
  return(exp(-5*t^2-t))
}
#[0.2,1.0]
T19_2 <- function(t){
  return(exp(-3*t+0.2))
}
rT19 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-0.4) & randomunifs < 1] <- -0.1+0.1*sqrt(1-20*log(randomunifs[randomunifs>=exp(-0.4) & randomunifs < 1]))
  randoms[randomunifs>=exp(-2.8) & randomunifs < exp(-0.4)] <- 1/15-(1/3)*log(randomunifs[randomunifs>=exp(-2.8) & randomunifs < exp(-0.4)])
  randoms[randomunifs<exp(-2.8)] <- Inf
  return(randoms)
}
#Treatment
#[0,0.2]
T20_1 <- function(t){
  return(exp(5*t^2-3*t))
}
#[0.2,1.0]
T20_2 <- function(t){
  return(exp(-t-0.2))
}
rT20 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-0.4) & randomunifs < 1] <- 0.3-0.1*sqrt(9+20*log(randomunifs[randomunifs>=exp(-0.4) & randomunifs < 1]))
  randoms[randomunifs>=exp(-1.2) & randomunifs < exp(-0.4)] <- -1/5-log(randomunifs[randomunifs>=exp(-1.2) & randomunifs < exp(-0.4)])
  randoms[randomunifs<exp(-1.2)] <- Inf
  return(randoms)
}
#16
#Control
#[0,0.25]
T21_1 <- function(t){
  return(exp(-t))
}
#[0.25,1.0]
T21_2 <- function(t){
  return(exp(-3*t+0.5))
}
rT21 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-0.25) & randomunifs < 1] <- -log(randomunifs[randomunifs>=exp(-0.25) & randomunifs < 1])
  randoms[randomunifs>=exp(-2.5) & randomunifs < exp(-0.25)] <- 1/6-(1/3)*log(randomunifs[randomunifs>=exp(-2.5) & randomunifs < exp(-0.25)])
  randoms[randomunifs<exp(-2.5)] <- Inf
  return(randoms)
}
#Treatment
T22 <- function(t){
  return(exp(-2*t))
}
rT22 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-2) & randomunifs < 1] <- -0.5*log(randomunifs[randomunifs>=exp(-2) & randomunifs < 1])
  randoms[randomunifs<exp(-2)] <- Inf
  return(randoms)
}
#17
#Control
#[0,0.8]
T23_1 <- function(t){
  return(exp(-1.25*t^2-t))
}
#[0.8,1.0]
T23_2 <- function(t){
  return(exp(-3*t+0.8))
}
rT23 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-1.6) & randomunifs < 1] <- -0.4+sqrt(4/25-0.8*log(randomunifs[randomunifs>=exp(-1.6) & randomunifs < 1]))
  randoms[randomunifs>=exp(-2.2) & randomunifs < exp(-1.6)] <- 4/15-(1/3)*log(randomunifs[randomunifs>=exp(-2.2) & randomunifs < exp(-1.6)])
  randoms[randomunifs<exp(-2.2)] <- Inf
  return(randoms)
}
#Treatment
#[0,0.8]
T24_1 <- function(t){
  return(exp(1.25*t^2-3*t))
}
#[0.8,1.0]
T24_2 <- function(t){
  return(exp(-t-0.8))
}
rT24 <- function(n){
  randomunifs <- runif(n,1e-9,1)
  randoms <- numeric(n)
  randoms[randomunifs>=exp(-1.6) & randomunifs < 1] <- 1.2-sqrt(36/25+0.8*log(randomunifs[randomunifs>=exp(-1.6) & randomunifs < 1]))
  randoms[randomunifs>=exp(-1.8) & randomunifs < exp(-1.6)] <- -0.8-log(randomunifs[randomunifs>=exp(-1.8) & randomunifs < exp(-1.6)])
  randoms[randomunifs<exp(-1.8)] <- Inf
  return(randoms)
}