##
MODELNR <- 1

## For scripting - replace MODELNR
input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
print(input)
eval(input)

library(DATRAS)
library(parallel)
library(mgcv)
load("splitModels.RData")
datafile = "../EBcodProcessedData.RData"
load(datafile)

dQ14$yqf = factor(paste(dQ14$Year, dQ14$Quarter,sep=":"))
Q1files <- dir("../anim", pattern=glob2rx("*Q1*.RData"),full=TRUE)
Q4files <- dir("../anim", pattern=glob2rx("*Q4*.RData"),full=TRUE)

alldens=array(NA,dim=c(nrow(EBarea.s),nlevels(dQ14$yqf),ncol(dQ14$N)))
for(f in 1:length(Q1files)){
    load(Q1files[f])
    alldens[,,f] = dens
}

alldens.split = alldens ##list()

splitIt<-function(ll,ct){
    pd = data.frame(  length_cm=ll, lon = EBarea.s$lon, ctime = ct) ##, haulId=d$haulId)
    ppp = predict(splitModel, newdata=pd,type="response",newdata.guaranteed=TRUE)
    ##cat(ll,"\n")
    ppp
}

first = which.min(dQ14$ctime)
last = which.max(dQ14$ctime)

f2n<-function(x)as.numeric(as.character(x))

firstq = f2n(dQ14[first]$Quarter)
lastq = f2n(dQ14[last]$Quarter)
firsty = f2n(dQ14[first]$Year)
lasty = f2n(dQ14[last]$Year)
q1time = median(dQ14$timeOfYear[ dQ14$Quarter=="1"] )
q4time = median(dQ14$timeOfYear[ dQ14$Quarter=="4"] )
ys = rep(firsty:lasty,each=2)
ys = ys + rep(c(q1time,q4time),length=length(ys))
if(firstq==4) ys=ys[-1]
if(lastq==1) ys=ys[-length(ys)]

cm.b = attr(dQ14,"cm.breaks")
cm.mid = cm.b[-length(cm.b)]+diff(cm.b)/2
cm.mid[ length(cm.mid) ] = 65 ## obs 
cm.mid[ 1 ] = 9


##for(sm in 1:length(splitModels)){
splitModel = splitModels[[MODELNR]]

for(ll in 1:(dim(alldens)[3]) ){
    ##cat("Split model ", sm,"/",length(splitModels), " length ", ll, "/", dim(alldens)[2],"\n")
    for( tt in 1:(dim(alldens)[2])){
        pwest = splitIt( cm.mid[ll], ys[tt] )
        alldens.split[,tt,ll] = pwest            
    }
}
##}

save(alldens.split,file=paste0("alldenssplit-",MODELNR,".RData"))
