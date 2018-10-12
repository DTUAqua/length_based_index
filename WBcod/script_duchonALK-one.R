SETTING <- 1

## For scripting - replace SETTING
input <- parse(text=Sys.getenv("SCRIPT_INPUT"))
print(input)
eval(input)


library(DATRAS)
library(parallel)
source("alkfun.R")
datafile = "../EBcodProcessedData.RData"
load(datafile)
library(mgcv)
load("splitModels.RData")

baselines=c(rep(0,5),rep(1,4))
settings = data.frame( modelnr = c(1:9,1:9,-1), baseline=c(rep(baselines,2),-1), fixedSplit=c(rep(NA,18),13), soft=c(rep(1,9),rep(0,9),0))

dQ14$yqf = factor(paste(dQ14$Year, dQ14$Quarter,sep=":"))
Q1files <- dir("../anim", pattern=glob2rx("*Q1*.RData"),full=TRUE)
Q4files <- dir("../anim", pattern=glob2rx("*Q4*.RData"),full=TRUE)

alldens=array(NA,dim=c(nrow(EBarea.s),nlevels(dQ14$yqf),ncol(dQ14$N)))
for(f in 1:length(Q1files)){
    load(Q1files[f])
    alldens[,,f] = dens
}

allageindicesQ4 = list()
allageindicesQ1 = list()

### Loop through split models
##for(setting in 1:nrow(settings)){
setting = SETTING
    cat("Doing setting ",setting, " out of ",nrow(settings),"\n")
    modelnr = settings$modelnr[setting]
    load(datafile)
    dQ14[[1]]=subset(dQ14[[1]], !is.na(Age))

    if(modelnr>0){
        cat("Splitting age data\n");
        sp2<-read.table("BITS1992_2017_baseline.txt",header=TRUE)
        
        if(settings$baseline[setting]==1) sp2<-read.table("BITS1992_2017_baseline_juv.txt",header=TRUE)
        
        sp2$day <- as.numeric( unlist( lapply( sp2$capture_date, substr, start=7,stop=8) ) )
        sp2$rn <- 1:nrow(sp2)
        oldcoords = sp2[,c("lat","long")]
        sp2 = sp2[,which( !names(sp2)%in%c("lat","long"))]
        
        ## find lat/lons
        DEdata = subset(dQ14, Country=="GFR")
        tmp = merge( dQ14[[1]], dQ14[[2]][,c("Month","Day","StatRec","lat","lon","haul.id")],by="haul.id",all.x=TRUE,suffixes=c("",".y"))
        sp2$hid = paste(sp2$year,sp2$month, sp2$day, sp2$rectangle, sp2$station, floor(sp2$length_cm),sep=":")
        tmp$hid = paste(tmp$Year, tmp$Month, tmp$Day, tmp$StatRec, tmp$StNo, tmp$LngtCm,sep=":" )
        tmp$rn=1:nrow(tmp)
        
        dd = merge(sp2, tmp[,c("hid","lon","lat","haul.id")], by="hid",all.x=TRUE)
        dd$age[dd$age<0]=NA
        dd = dd[ !duplicated(dd$rn), ]
        dd = subset(dd, !is.na(haul.id))
        
        tmp2 = data.frame(haul.id = dd$haul.id, LngtCm = dd$length_cm, Age = dd$age, NoAtALK = 1, stock=dd$stock )
        tmp2$ProbWest=NA
        tmp2$ProbWest[ tmp2$stock=="WEST" ] = 1
        tmp2$ProbWest[ tmp2$stock=="EAST" ] = 0
        
        
        ## replace original German data with Franziska's
        oldim=dim(dQ14[[1]])
        dQ14[[1]] <- subset(dQ14[[1]], !haul.id%in%dd$haul.id | LngtCm<20) 
        dQ14[[1]] <- dQ14[[1]][,c("haul.id","LngtCm","Age","NoAtALK")]
        dQ14[[1]]$stock=NA
        dQ14[[1]]$ProbWest=NA
        
        dQ14[[1]] <- rbind( dQ14[[1]], tmp2)
        dQ14[[1]]$stock=factor(dQ14[[1]]$stock)
        ## ##################
        ## Predict stock for un-stock'ed age data
        ## #######

        dQ14[[1]]$length_cm = dQ14[[1]]$LngtCm
        print(nrow(dQ14[[1]]))
        dQ14[[1]] = merge( dQ14[[1]], dQ14[[2]][,c("Year","Month","Day","Quarter","lon","haul.id")],by="haul.id",all.x=TRUE,suffixes=c("",".y"))
        print(nrow(dQ14[[1]]))
        
        dQ14[[1]]$ctime = as.numeric(as.character(dQ14[[1]]$Year)) + (dQ14[[1]]$Month-1)/12 + (dQ14[[1]]$Day-1)/365
        dQ14[[2]]$ctime = as.numeric(as.character(dQ14[[2]]$Year)) + (dQ14[[2]]$Month-1)/12 + (dQ14[[2]]$Day-1)/365
        
        sel=which(is.na(dQ14[[1]]$stock) & !is.na(dQ14[[1]]$Age))
        
        
        xtabs(NoAtALK ~ Year + Age + Quarter,data=dQ14[[1]]) ##OK
        xtabs(!is.na(Age) ~ Year + Age + Quarter,data=dQ14[[1]]) ##OK
        xtabs(is.na(stock) ~ Year + Age + Quarter,data=dQ14[[1]])
        xtabs(is.na(stock) & !is.na(Age) ~ Year + Age + Quarter,data=dQ14[[1]])
        summary(dQ14[[1]][sel,])
        xtabs( ~ Year + Age + Quarter,data=dQ14[[1]][sel,])
        xtabs(NoAtALK ~ Year + Age + Quarter,data=dQ14[[1]][,])
        xtabs( ~ Year + Age + Quarter,data=dQ14[[1]][sel,])
        
        probw=predict(splitModels[[modelnr]],newdata=dQ14[[1]][sel,],type="response")
        dQ14[[1]][sel,"ProbWest"]=probw
        xtabs(is.na(ProbWest) ~ Year + Age + Quarter,data=dQ14[[1]]) ## all are now NA!
        xtabs(!is.na(ProbWest) ~ Year + Age + Quarter,data=dQ14[[1]]) ## all are now NA!
        summary(dQ14[[1]][is.na(dQ14[[1]]$ProbWest),])

        dQ14[[1]]$NoAtALK = dQ14[[1]]$NoAtALK*dQ14[[1]]$ProbWest ## Weighting
    }##fi modelnr>0

###########
    alldens.prod = alldens
    
    if(modelnr<0){ 
        ## hard fixed split
        fs = settings$fixedSplit[setting]
        
        dQ14 = subset(dQ14,lon<fs)
        sel = which(EBarea.s$lon<fs)
        EBarea.s = EBarea.s[sel,]
        alldens.prod = alldens.prod[sel,,]
    } else {   
        load(paste0("alldenssplit-",modelnr,".RData"))
        if(settings$soft[setting] == 1){
            alldens.prod = alldens * alldens.split
        } else {
            ## hard split
            alldens.prod = alldens * as.numeric(alldens.split>0.5)
        }
    }


    dQ14[[1]] = subset(dQ14[[1]], !is.na(Age) & !is.na(NoAtALK))
    
    
    cms=attr(dQ14,"cm.breaks");
    Q1d = subset(dQ14,Quarter==1)
    Q4d = subset(dQ14,Quarter==4)
    

    if(FALSE){
        lo = matrix(c(1,1,2,2,
                      3,3,4,4),2,4,byrow=TRUE)
        
        nf<-layout(lo)
        layout.show(nf)
        
        
        timestep=5
        ## Abundance unsplit
        concT = surveyIndex:::concTransform(log(alldens[,timestep,16]))
        mycolors = rev(heat.colors(10))
        zFac = cut(concT, 0:length(mycolors)/length(mycolors))
        plot(EBarea.s$utm.x,EBarea.s$utm.y,col=mycolors[zFac],pch=15,cex=0.8,xlab="utm.x",ylab="utm.y",main="East+West Abundance")
        points(coastlines$x,coastlines$y,pch=".",col=3)
        
        ## split
        my.palette2<-colorRampPalette(c("red","white","blue"))(100)
        zFac = cut(alldens.split[sel,timestep,16], 0:length(my.palette2)/length(my.palette2))
        plot(EBarea.s$utm.x,EBarea.s$utm.y,col=my.palette2[zFac],pch=15,cex=0.8,xlab="utm.x",ylab="utm.y",main="East-West Split")
        points(coastlines$x,coastlines$y,pch=".",col=3)
        
        ## product
        west30cm = alldens[,timestep,16] * alldens.split[sel,timestep,16]
        concT = surveyIndex:::concTransform(log(west30cm))
        zFac = cut(concT, 0:length(mycolors)/length(mycolors))
        plot(EBarea.s$utm.x,EBarea.s$utm.y,col=mycolors[zFac],pch=15,cex=0.8,xlab="utm.x",ylab="utm.y",main="West Abundance")
        points(coastlines$x,coastlines$y,pch=".",col=3)
        
        ## product
        east30cm = alldens[,timestep,16] * (1-alldens.split[sel,timestep,16])
        concT = surveyIndex:::concTransform(log(east30cm))
        zFac = cut(concT, 0:length(mycolors)/length(mycolors))
        plot(EBarea.s$utm.x,EBarea.s$utm.y,col=mycolors[zFac],pch=15,cex=0.8,xlab="utm.x",ylab="utm.y",main="East Abundance")
        points(coastlines$x,coastlines$y,pch=".",col=3)

        
        ## ############################
        ## Hard split alternative
        ## ############################

        ## Abundance unsplit
        concT = surveyIndex:::concTransform(log(alldens[,timestep,16]))
        mycolors = rev(heat.colors(10))
        zFac = cut(concT, 0:length(mycolors)/length(mycolors))
        plot(EBarea.s$utm.x,EBarea.s$utm.y,col=mycolors[zFac],pch=15,cex=0.8,xlab="utm.x",ylab="utm.y",main="East+West Abundance")
        points(coastlines$x,coastlines$y,pch=".",col=3)

        ## split
        my.palette2<-colorRampPalette(c("red","white","blue"))(100)
        zFac = cut(alldens.split[sel,timestep,16], 0:length(my.palette2)/length(my.palette2))
        plot(EBarea.s$utm.x,EBarea.s$utm.y,col=my.palette2[zFac],pch=15,cex=0.8,xlab="utm.x",ylab="utm.y",main="East-West Split")
        points(coastlines$x,coastlines$y,pch=".",col=3)

        ## product

        west30cm = alldens[,timestep,16] * as.numeric(alldens.split[sel,timestep,16]>0.5)
        concT = surveyIndex:::concTransform(log(west30cm))
        zFac = cut(concT, 0:length(mycolors)/length(mycolors))
        plot(EBarea.s$utm.x,EBarea.s$utm.y,col=mycolors[zFac],pch=15,cex=0.8,xlab="utm.x",ylab="utm.y",main="West Abundance")
        points(coastlines$x,coastlines$y,pch=".",col=3)

        ## product
        east30cm = alldens[,timestep,16] * as.numeric(alldens.split[sel,timestep,16]<0.5)
        concT = surveyIndex:::concTransform(log(east30cm))
        zFac = cut(concT, 0:length(mycolors)/length(mycolors))
        plot(EBarea.s$utm.x,EBarea.s$utm.y,col=mycolors[zFac],pch=15,cex=0.8,xlab="utm.x",ylab="utm.y",main="East Abundance")
        points(coastlines$x,coastlines$y,pch=".",col=3)

    }

    cat("Fitting ALK Q4\n");
    mf = "" 
    ack=TRUE;
    useBICs=TRUE;
    varCofs=FALSE;
    maxKs=50;
    mc.cores=3
    cmSize=1
    
    Q4d = subset(Q4d,Year!="1991") ## no age samples in 1991
    xtabs(!is.na(Age)~Year+Age+Quarter,data=Q4d[[1]])
    Q4ys = split(Q4d,Q4d$Year)
    Q1ages=1:5
    Q4ages=0:5
    ages=Q4ages

    
    d.ALK= mclapply(Q4ys,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=ack,useBIC=useBICs,varCof=varCofs,maxK=maxKs,splineType=2,mc.cores=mc.cores)

    plotIt=TRUE
    Q1times = seq(3, dim(dens)[2], by=2) ## OBS, skip 1991
    Q4times = seq(4, dim(dens)[2], by=2)
    
    ageindicesQ1 = matrix(NA,length(Q1times),length(Q1ages))
    ageindicesQ4 = matrix(NA,length(Q4times),length(Q4ages))

    nagefun<-function(i,t){
        alldens.prod[i,t,]%*%ALKs[[i]]
    }

    cat("Applying ALK Q4 \n");
    pdf("Q4maps_duchon.pdf")
    for(tt in Q4times){
        cat("Doing time ",which(Q4times==tt)," of ", length(Q4times),"\n")
        curALK=d.ALK[[which(Q4times==tt)]] 
        fakedat <- list(NULL,EBarea.s, N = matrix(NA,nrow=nrow(EBarea.s),ncol=ncol(dQ14$N)) )
        cm.b = attr(dQ14,"cm.breaks")
        cm.mid = cm.b[-length(cm.b)]+diff(cm.b)/2
        cm.mid[ length(cm.mid) ] = 65 ## obs hack
        cm.mid[ 1 ] = 9
        
        attr(fakedat,"cm.breaks")<-cm.mid
        
        attr(curALK,"data") <- fakedat
        
        ALKs=mclapply(1:nrow(EBarea.s),NageByHaul,x=curALK,returnALK=TRUE,mc.cores=3)
        
        nage = sapply(1:nrow(EBarea.s),nagefun,t=tt) 
        ageindicesQ4[which(Q4times==tt),] = rowSums(nage)
        if(plotIt){
            par(mfrow=c(2,3))
            for(aa in 1:6){
                concT = surveyIndex:::concTransform(log(nage[aa,]))
                mycolors = rev(heat.colors(6))
                zFac = cut(concT, 0:length(mycolors)/length(mycolors))
                plot(EBarea.s$utm.x,EBarea.s$utm.y,col=mycolors[zFac],pch=16,cex=0.4,main=paste(levels(dQ14$yqf)[tt],"age",aa))
                points(coastlines$x,coastlines$y,pch=".",col=3)
            }
        } 
    }
    dev.off()
    allageindicesQ4[[ setting ]] = ageindicesQ4
    
    ## #############
    ## Repeat Q1
    ## #############
    cat("Fitting ALK Q1\n");
    Q1d = subset(Q1d,Year!="1991") ## no age samples in 1991
    Q1d[[1]] = subset(Q1d[[1]],Age>0)
    ages=Q1ages
    Q1ys = split(Q1d,Q1d$Year)
    d.ALK= mclapply(Q1ys,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=ack,useBIC=useBICs,varCof=varCofs,maxK=maxKs,splineType=2,mc.cores=1)
    

    pdf("Q1maps_duchon.pdf")
    cat("Applying ALK Q1 \n");
    for(tt in Q1times){
        cat("Doing time ",which(Q1times==tt)," of ", length(Q1times),"\n")
        curALK=d.ALK[[which(Q1times==tt)]] 
        fakedat <- list(NULL,EBarea.s, N = matrix(NA,nrow=nrow(EBarea.s),ncol=ncol(dQ14$N)) )
        cm.b = attr(dQ14,"cm.breaks")
        cm.mid = cm.b[-length(cm.b)]+diff(cm.b)/2
        cm.mid[ length(cm.mid) ] = 65 ## obs hack
        cm.mid[ 1 ] = 9
        
        attr(fakedat,"cm.breaks")<-cm.mid
        
        attr(curALK,"data") <- fakedat
        
        ALKs=mclapply(1:nrow(EBarea.s),NageByHaul,x=curALK,returnALK=TRUE,mc.cores=3)
        
        nage = sapply(1:nrow(EBarea.s),nagefun,t=tt) 
        ageindicesQ1[which(Q1times==tt),] = rowSums(nage)
        if(plotIt){
            par(mfrow=c(2,3))
            for(aa in Q1ages){
                concT = surveyIndex:::concTransform(log(nage[aa,]))
                mycolors = rev(heat.colors(6))
                zFac = cut(concT, 0:length(mycolors)/length(mycolors))
                plot(EBarea.s$utm.x,EBarea.s$utm.y,col=mycolors[zFac],pch=16,cex=0.4,main=paste(levels(dQ14$yqf)[tt],"age",aa))
                points(coastlines$x,coastlines$y,pch=".",col=3)
            }
        } 
    }
    dev.off()
    ##allageindicesQ1[[ setting ]] = ageindicesQ1

##}## rof settings


save(ageindicesQ1,ageindicesQ4,file=paste0("ageindices_duchon-",SETTING,".RData"))
