## Default inputs
SPECIES  <- "Gadus morhua"
QUARTER  <- 4
KM       <- 6
MINSIZE  <- 4
MAXSIZE  <- 120
MINYEAR  <- 1991
MAXYEAR  <- 2019
BY       <- 1
CMGROUP <- MINSIZE+1
DATFILE  <- "EBcodProcessedData.RData"
PVECFILE <- "pvecs.RData"
OUTFILE  <- paste0("results", QUARTER,"-cm",CMGROUP,".RData")

## For scripting
input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
print(input)
eval(input)

## Load data
library(DATRAS)
d <- local({
    load(DATFILE)
    dQ14
})
stopifnot( class(d) == "DATRASraw" )

## Load fine scale prediction grid
area <- local({
    load(DATFILE)
    EBarea.s
})
stopifnot( class(area) == "data.frame" )

## Integration weights
pvecs <- local({
    load(PVECFILE)
    pvecs
})
stopifnot( class(pvecs) == "list" )
pvecs <- lapply(pvecs, function(x) x[CMGROUP, , ])

## Make grid
library(gridConstruct)
grid <- gridConstruct(d,km=KM)
## plot(grid)
## map("worldHires",add=TRUE)

## Lookup grid cell of fine scale (lon,lat) coords
area$position <- gridFactor(area, grid)

## Common transformations applied to both sets:
## (re-scaled depth)
rescale <- function(x) (x - 50) / 20
d   $DepthRS <- rescale(   d$Depth)
area$DepthRS <- rescale(area$Depth)
d   $logHaulDur <- log(   d$HaulDur)
area$logHaulDur <- log(area$HaulDur)
d$timeOfDay <- factor(d$timeOfDay)
area$timeOfDay <- 1
area$timeOfDay <- factor(area$timeOfDay, levels=levels(d$timeOfDay))
area$Quarter <- 1
area$Quarter <- factor(area$Quarter, levels=levels(d$Quarter))

## Data subset
d$N <- d$N[,CMGROUP,drop=FALSE]
##d <- addSpectrum(d,cm.breaks=seq(CMGROUP-1,CMGROUP+1,by=BY))
d$haulid <- d$haul.id
##d <- subset(d, Quarter == QUARTER, Gear != "GRT")
d <- subset(d, Year %in% MINYEAR:MAXYEAR )
##d <- subset(d, 25<HaulDur & HaulDur<35 )
d <- as.data.frame(d)

library(mapdata)

## Set up time factor (careful with empty factor levels ! )
years <- as.numeric(levels(d$Year))
## d$time <- factor(d$Year, levels=min(years):max(years))
## HACK: Time = year + quarter / 10
timeLevels <- sort(outer(min(years):max(years), c(1,4)/10, "+"))
time <- as.numeric(as.character(d$Year)) + as.numeric(as.character(d$Quarter)) / 10
timeLevels <- timeLevels[which(min(time)==timeLevels):which(max(time)==timeLevels)]
d$time <- factor(time, levels=timeLevels)

## Set up spatial factor and grid
## grid <- gridConstruct(d,nearestObs=100)
d$position <- gridFactor(d,grid)
Q0 <- -attr(grid,"pattern")
diag(Q0) <- 0; diag(Q0) <- -rowSums(Q0)
I <- .symDiagonal(nrow(Q0))

## Set up design matrix
## TODO: Gear by sizeGroup
##A <- sparse.model.matrix( ~ sizeGroup:time + Gear - 1, data=d)
A0 <- sparse.model.matrix( ~ time - 1, data=d)
form <- ~ Gear - 1 + Quarter:DepthRS + Quarter:I(DepthRS^2) + logHaulDur + timeOfDay
## If you want to test the effect of Quarter specific depth:
## form <- ~ Gear - 1 + DepthRS + I(DepthRS^2) + Quarter:DepthRS + Quarter:I(DepthRS^2) + logHaulDur + timeOfDay
A <- sparse.model.matrix(form, data=d)
##A <- A[ , colnames(A) != "GearTVL" , drop=FALSE]

beta     <- rep(0, ncol(A)); names(beta) <- colnames(A)
beta_map <- beta * 0 + seq_along(beta)
beta["GearTVL"]        <- 0  ## Reference !!!
beta["logHaulDur"]     <- 1  ## Offset !!!
beta_map["logHaulDur"] <- NA ## FIXED
beta_map["GearTVL"]    <- NA ## FIXED
map <- list(beta = factor(beta_map))

B <- cbind2(A[, names(beta_map[!is.na(beta_map)]) ],A0); B <- t(B)%*%B
if(min(eigen(B)$val)<1e-8)stop("Singular B")

## Design matrix for *spatial* prediction
area$Gear <- factor(area$Gear, levels = levels(d$Gear))
Apredict <- sparse.model.matrix( form, data=area)

data <- as.list(d)
data$A <- A
data$I <- I
data$Q0 <- Q0

data$Apredict <- Apredict
data$pos_predict <- area$pos
data$p <- matrix(1, nrow(data$Apredict), nlevels(data$time))

data <- data[!sapply(data,is.character)]
data <- data[!sapply(data,is.logical)]

## Prior std dev on fixed effects (for robustness only)
data$huge_sd <- 100

## Can't estimate size corr if only one sizegroup:
if (nlevels(d$sizeGroup) <= 1) {
    map$tphi_size <- factor(NA)
}

## Assert pvec dimensions are right
stopifnot(all( sapply(pvecs, nrow) == nrow(data$Apredict) ) )
stopifnot(all( sapply(pvecs, ncol) == nlevels(data$time) ) )

library(TMB)
compile("model.cpp")
dyn.load(dynlib("model"))
obj <- MakeADFun(
    data=data,
    parameters=list(
        logdelta= 0 ,
        logkappa= 0 ,
        tphi_time= c(0.01, 0.8) ,
        tphi_size = 0 ,
        logsigma= 0 ,
        beta= beta ,
        eta= array(0, c(nrow(Q0), nlevels(d$time), nlevels(d$sizeGroup) ) ),
        etanug= array(0, c(nlevels(d$sizeGroup), nlevels(d$haulid) ) ),
        etamean = array(0, c(nlevels(d$sizeGroup) , nlevels(d$time)) )
        ),
    DLL="model",
    random=c("eta","etanug","etamean","beta"),
    map=map
    )

print(obj$par)
runSymbolicAnalysis(obj)
system.time(obj$fn())

system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
hessian <- optimHess(opt$par, obj$fn, obj$gr)

## Quarter 1 sdreport
sdr1 <- lapply(pvecs, function(p) {
    obj$env$data$p <- p
    sdreport(obj, hessian.fixed=hessian)
})

## Quarter 4 sdreport
sdr4 <- lapply(pvecs, function(p) {
    obj$env$data$p <- p
    area$Quarter[] <- "4"
    obj$env$data$Apredict <- as(sparse.model.matrix( form, data=area), "dgTMatrix")
    sdreport(obj, hessian.fixed=hessian)
})

if(FALSE) { ## Skip bias correction for now
    system.time(sdr2 <- sdreport(obj,bias.correct=TRUE,hessian.fixed=solve(sdr$cov.fixed)))
    mat <- summary(sdr2,"report")
    rownames(mat) <- NULL
    df <- as.data.frame(mat)
    df$unbiased <- sdr2$unbiased$value
    df <- cbind(df,
                expand.grid(sizeGroup=levels(d$sizeGroup),time=levels(d$time))
                )
    xtabs(Estimate ~ sizeGroup + time,data=df)
    x1 <- xtabs(N ~ sizeGroup + time,data=d) / xtabs( ~ sizeGroup + time,data=d)
    x2 <- xtabs(Estimate ~ sizeGroup + time,data=df)
    x3 <- xtabs(unbiased ~ sizeGroup + time,data=df)
}

obj$env$L.created.by.newton <- NULL
obj$env$spHess <- NULL

save(sdr1, sdr4, obj, grid, file=OUTFILE)
## Set choleski and sp_hess to NULL before saving obj.

## pdf("plotQ4.pdf")
## matplot(x1,type="l",main="Raw average")
## matplot(x2,type="l",main="Posterior mode")
## matplot(x3,type="l",main="Posterior mean")
## dev.off()
## rep <- obj$report(obj$env$last.par.best)
