library(DATRAS)
load("EBcod.RData"); d <- dAll
library(lgc)
grid <- gridConstruct(d,km=50)
## plot(grid)
## map("worldHires",add=TRUE)
d <- addSpectrum(d,cm.breaks=seq(4,120,by=2))
d$haulid <- d$haul.id
d <- subset(d, Quarter == 4, Gear != "GRT")
d <- subset(d, Year %in% 2000:2015 )
d <- subset(d, 25<HaulDur & HaulDur<35 )
d <- as.data.frame(d)
## 24 * 78 * 58
library(lgc)
library(mapdata)

## Set up time factor (careful with empty factor levels ! )
years <- as.numeric(levels(d$Year))
d$time <- factor(d$Year, levels=min(years):max(years))

## Set up spatial factor and grid
## grid <- gridConstruct(d,nearestObs=100)
d$position <- gridFactor(d,grid)
Q0 <- -attr(grid,"pattern")
diag(Q0) <- 0; diag(Q0) <- -rowSums(Q0)
I <- .symDiagonal(nrow(Q0))

## Set up design matrix
## TODO: Gear by sizeGroup
##A <- sparse.model.matrix( ~ sizeGroup:time + Gear - 1, data=d)
A0 <- sparse.model.matrix( ~ sizeGroup:time - 1, data=d)
A <- sparse.model.matrix( ~ Gear - 1, data=d)
A <- A[,-which.max(table(d$Gear)),drop=FALSE]

B <- cbind2(A,A0); B <- t(B)%*%B
if(min(eigen(B)$val)<1e-8)stop("Singular B")

data <- as.list(d)
data$A <- A
data$I <- I
data$Q0 <- Q0

data <- data[!sapply(data,is.character)]
data <- data[!sapply(data,is.logical)]

library(TMB)
openmp(8)
compile("model.cpp")
dyn.load(dynlib("model"))
obj <- MakeADFun(
    data=data,
    parameters=list(
        logdelta= 0 ,
        logkappa= 0 ,
        tphi_time= 0 ,
        tphi_size = 0 ,
        logsigma= 0 ,
        beta= rep(0, ncol(A)) ,
        eta= array(0, c(nrow(Q0), nlevels(d$time), nlevels(d$sizeGroup) ) ),
        etanug= array(0, c(nlevels(d$sizeGroup), nlevels(d$haulid) ) ),
        etamean = array(0, c(nlevels(d$sizeGroup) , nlevels(d$time)) )
        ),
    DLL="model",
    random=c("eta","etanug","etamean","beta")
    )

print(obj$par)
runSymbolicAnalysis(obj)
system.time(obj$fn())

system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))

## system.time(sdr <- sdreport(obj))
## summary(sdr,"fixed")
## summary(sdr,"report")

logindex <- obj$report()$logindex
rownames(logindex) <- levels(d$sizeGroup)
colnames(logindex) <- levels(d$time)

save(logindex,file="result_Q4.RData")
