RESFILE <- "results/cod-Q14-cm0024.RData"
OUTFILE <- "anim.pdf"
QUARTERS <- "1"

## For scripting
input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
print(input)
eval(input)

load(RESFILE)
load("EBcodProcessedData.RData")
area <- EBarea.s

library(TMB)
library(gridConstruct)
library(sp)
library(rgdal)

## Get reported density
dyn.load(dynlib("model"))
lpb <- obj$env$last.par.best
obj$env$data$p[] <- 1
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
grid.size.m = 2500
utmx <- seq(min(area$utm.x)-1,max(area$utm.x)+grid.size.m,by=grid.size.m)
utmy <- seq(min(area$utm.y)-1,max(area$utm.y)+grid.size.m,by=grid.size.m)
area$cututmx = cut(area$utm.x,utmx)
area$cututmy = cut(area$utm.y,utmy)
for (i in 1:ncol(rep$dens)) {
    quarter <- substring(times[i],6,6)
    if(quarter %in% QUARTERS) {
        u <- concTransform(log(rep$dens[,i]))
        area$u <- u
        mat <- xtabs(u ~ cututmx + cututmy, data=area)
        not.missing <- xtabs( ~ cututmx + cututmy, data=area)
        breaks <- seq(0,1,length=11)
        mat[] <- cut(mat,breaks)
        mat[!as.logical(not.missing)] <- NA
        image(utmx, utmy, mat, col=rev(heat.colors(10)))
        points(coastlines, pch=".", col="green")
        ## Trawl points
        LongLatToUTM<-function(x,y,zone=33){
            xy <- data.frame(ID = 1:length(x), X = x, Y = y)
            coordinates(xy) <- c("X", "Y")
            proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
            res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
            return(as.data.frame(res))
        }
        long2UTMzone <- function(long) {
            (floor((long + 180)/6) %% 60) + 1
        }
        lon <- obj$env$data$lon[obj$env$data$time + 1 == i]
        lat <- obj$env$data$lat[obj$env$data$time + 1 == i]
        points(LongLatToUTM(lon,lat)[-1], pch="+", cex=.8)
        title(times[i])
    }
}
dev.off()
