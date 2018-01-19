RESFILE <- "results/cod-Q14-cm0012.RData"
OUTFILE <- "post.RData"

## For scripting
input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
print(input)
eval(input)

load(RESFILE)

library(TMB)
library(gridConstruct)

## Index
logindex1 <- summary(sdr1,"report")
logindex4 <- summary(sdr4,"report")

times <- levels(obj$env$data$time)
rownames(logindex1) <- times
rownames(logindex4) <- times

timesQ1 <- substring(times,5,6) == ".1"
logindex <- logindex4
logindex[timesQ1, ] <- logindex1[timesQ1, ]

## Effects
beta <- as.list(sdr1, "Est")$beta
betaSD <- as.list(sdr1, "Std")$beta

save(logindex, beta, betaSD, file=OUTFILE)
