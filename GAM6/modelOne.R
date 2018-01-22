## GAM6 - remove depth effect 
CMGROUP <- 3

## For scripting - replace CMGROUP
input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
print(input)
eval(input)

library(mgcv)
library(DATRAS)
library(sp)
library(rgdal)
library(maps)
library(mapdata)

datafile = "../EBcodProcessedData.RData"
load(datafile)

OUTFILE = paste0("results/",CMGROUP,".RData")

source("fun.R")

dQ14$DepthRS = (dQ14$Depth - 50)/20
EBarea.s$Ship = names(sort(table(dQ14$Ship),TRUE))[1]

dQ14$Gear = factor( dQ14$Gear, levels=names(sort(table(dQ14$Gear),TRUE))) 

mymodel ="Quarter+s(Year,Quarter,bs='re') + Gear + s(utm.x,utm.y,bs='ts',k=144,by=Quarter) + ti(ctime,utm.x,utm.y,d=c(1,2),bs=c('ts','ts'),k=c(16,9)) + s(TimeShotHour,k=6,bs='cc') + offset(log(HaulDur))";

tsel = which( !duplicated(dQ14$ctime))

tselQ1 = which( !duplicated(dQ14$ctime) & dQ14$Quarter == "1" )
tselQ4 = which( !duplicated(dQ14$ctime) & dQ14$Quarter == "4" )


knots=list(TimeShotHour=seq(0,24,length=6))

EBarea.s$DepthRS = (EBarea.s$Depth - 50)/20

print( system.time( myidx <- getSurveyIdxLength(dQ14,CMGROUP,NULL,nBoot=1000,predD=EBarea.s,model=mymodel,knots=knots,mc.cores=1,pred.ctimes=dQ14$ctime[tsel],pred.years=dQ14$Year[tsel],pred.quarter=dQ14$Quarter[tsel]) ) )

plotDepthEffect<-function(model, nBoot=4000){
    pd = EBarea.s[rep(1,100),]
    pd$Year = levels(dQ14$Year)[1]
    pd$DepthRS = seq(min(dQ14$DepthRS),max(dQ14$DepthRS),length.out=100)
    pd$Quarter="1"
    pd2 = pd
    pd2$Quarter="4"
    pd = rbind(pd,pd2)

    depthSel = grepl("DepthRS",names(coef(model)))
    depthSel[1] = TRUE
    Xp.1=predict(model,newdata=pd,type="lpmatrix")[,depthSel];
    library(MASS)
    brp.1=mvrnorm(n=nBoot,coef(model)[depthSel],model$Vp[depthSel,depthSel]);
    OS.pos = matrix(0,nrow(pd),nBoot);
    terms.pos=terms(model)
    if(!is.null(model$offset)){
        off.num.pos <- attr(terms.pos, "offset")
        for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos, 
                                                            "variables")[[i + 1]], pd)
    }
    rep1=exp(Xp.1%*%t(brp.1)+OS.pos);
    
    Q1 = rep1[1:100,]
    Q4 = rep1[101:200,]
    for(i in 1:nBoot) Q1[,i] = Q1[,i] 

    Q1.up = apply(Q1,1,quantile,probs=0.975)
    Q1.lo = apply(Q1,1,quantile,probs=0.025)
    Q4.up = apply(Q4,1,quantile,probs=0.975)
    Q4.lo = apply(Q4,1,quantile,probs=0.025)

    x = pd2$DepthRS*20+50
    plot(x,rowMeans(Q1),type="l",ylim=c(0,max(c(Q1.up,Q4.up))),xlab="Depth",ylab="logEffect")
    lines(x,rowMeans(Q4),type="l",col=2)
    lines(x,Q1.up,lty=2)
    lines(x,Q1.lo,lty=2)
    lines(x,Q4.lo,col=2,lty=2)
    lines(x,Q4.up,col=2,lty=2)
    legend("topright",legend=c("Q1","Q4"),lty=1,col=c(1,2))
    
}

pdf(paste0("results/covariates-",CMGROUP,".pdf"))
plot(myidx$models[[1]],select=6,rug=TRUE) ## time of day
plotDepthEffect(myidx$models[[1]])

gearCoef<-function(m) coef(m)[grepl("Gear",names(coef(m)))]
gearCoefSd<-function(m) sqrt(diag(vcov(m)))[grepl("Gear",names(coef(m)))]

namedCoef<-function(m,name="Gear") {
    list( mu=coef(m)[grepl(name,names(coef(m)))],
         sd=sqrt(diag(vcov(m)))[grepl(name,names(coef(m)))] )
}

gc=gearCoef(myidx$models[[1]])
gc.sd=gearCoefSd(myidx$models[[1]])
dc = namedCoef(myidx$models[[1]],"DepthRS")

sink(file=paste0("results/summary-",CMGROUP,".txt"))
summary(myidx$models[[1]])
gam.check(myidx$models[[1]])
sink()
dev.off()

my.cols=rev(heat.colors(10))

concTransform<-function (x) 
{
    i <- order(x)
    ys <- sort(exp(x))
    p <- ys/sum(ys)
    x[i] <- cumsum(p)
    x
}

plotMapFit<-function(m, pred.ctimes,pred.years,pred.quarter){
    for(ctime in pred.ctimes){
        sel = which( dQ14$ctime == ctime )
        year = dQ14$Year[sel[1]]
        quarter = pred.quarter[ which(pred.ctimes==ctime) ]
        
        fit <- exp(m$gPreds2[[1]][[ which(pred.ctimes==ctime) ]])
        concT = concTransform(log(fit))
 
        plot(EBarea.s$utm.x,EBarea.s$utm.y,col=my.cols[cut(concT,0:length(my.cols)/length(my.cols))],pch=15,cex=0.5,main=paste(year,quarter,colnames(dQ14$N)[CMGROUP]))
        
        points(coastlines,pch=".",col=3)

        points( dQ14$utm.x[sel], dQ14$utm.y[sel], pch=1, cex = sqrt( dQ14$N[,CMGROUP][sel]/(4*colMeans(dQ14$N)[CMGROUP]  )))
        sel2 = which( dQ14$ctime == ctime & dQ14$N[,CMGROUP] == 0)
        points( dQ14$utm.x[sel2], dQ14$utm.y[sel2], pch="+", cex = 1,col="blue" )
        
    }
}


pdf(paste0("results/mapsQ1-cm",CMGROUP,".pdf"))
 plotMapFit(myidx,pred.ctimes=dQ14$ctime[tselQ1],pred.years=dQ14$Year[tselQ1],pred.quarter=dQ14$Quarter[tselQ1])
dev.off()


pdf(paste0("results/mapsQ4-cm",CMGROUP,".pdf"))
 plotMapFit(myidx,pred.ctimes=dQ14$ctime[tselQ4],pred.years=dQ14$Year[tselQ4],pred.quarter=dQ14$Quarter[tselQ4])
dev.off()


out<-list(idx = myidx$idx, gc=gc, gc.sd=gc.sd, dc=dc, lo=myidx$lo, up=myidx$up )
save(out,file=paste0("results/out-cm",CMGROUP,".RData"))

