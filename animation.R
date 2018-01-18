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
rep <- obj$report(lpb)

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
