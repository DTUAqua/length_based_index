RESFILE <- "results/cod-Q14-cm0008.RData"
OUTFILE <- "anim.pdf"

## For scripting
input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
print(input)
eval(input)

load(RESFILE)
load("EBcodProcessedData.RData")
area <- EBarea.s

library(TMB)
library(gridConstruct)

## Get reported density
dyn.load(dynlib("model"))
lpb <- obj$env$last.par.best
obj$env$type <- "Fun" ## No AD !
rep4 <- obj$report(lpb) ## Q4 !
## Swap Q1/Q4 columns:
cn <- colnames(obj$env$data$Apredict)
cnswap <- sub("Quarter4","quarter1",cn)
cnswap <- sub("Quarter1","quarter4",cnswap)
cnswap <- sub("quarter","Quarter",cnswap)
Aswap <- obj$env$data$Apredict[, cnswap]
obj$env$data$Apredict <- Aswap
rep1 <- obj$report(lpb) ## Q1 !
##
times <- levels(obj$env$data$time)
rep <- rep4
timesQ1 <- substring(times,5,6) == ".1"
rep$dens[, timesQ1 ] <- rep1$dens[ , timesQ1 ]

## For each time, attach spatial 'index' to area and plot:
pdf(OUTFILE)
for (i in 1:ncol(rep$dens)) {
    u <- concTransform(log(rep$dens[,i]))
    area$u <- u
    mat <- xtabs(u ~ cututmx + cututmy, data=area)
    not.missing <- xtabs( ~ cututmx + cututmy, data=area)
    mat[!as.logical(not.missing)] <- NA
    image(mat)
    title(levels(obj$env$data$time)[i])
}
dev.off()
