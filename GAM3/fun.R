getSurveyIdxLength <-
    function(x,lengths,myids,nBoot=1000,mc.cores=2,method="REML",predD,
             model,pred.ctimes,pred.years,pred.quarters,knots=NULL){
        
        
        if(length(model)<length(lengths)) stop(" length(model) < length(lengths)");
        
        x[[1]]$Year=as.factor(x[[1]]$Year);
        x[[2]]$Year=as.factor(x[[2]]$Year);
        models=list()
        gPreds=list() ##last data year's predictions
        gPreds2=list() ## all years predictions
        require(mgcv)
        require(parallel)
        ##yearNum=as.numeric(as.character(x$Year));
        ##yearRange=min(yearNum):max(yearNum);

        ## Choose most frequent gear as reference gear
        gearNames=names(xtabs(~Gear,data=x[[2]]))
        myGear=names(xtabs(~Gear,data=x[[2]]))[which.max(xtabs(~Gear,data=x[[2]]))]
        print(paste("using",myGear,"as reference gear"))
        resMat=matrix(NA,nrow=length(pred.ctimes),ncol=length(lengths));
        upMat=resMat;
        loMat=resMat;
        do.one.a<-function(a){
            ddd=x[[2]]; ddd$dum=1.0;
            ddd$y=ddd$N[,lengths[a]]
            
            formu = as.formula( paste( "y ~",model[a]));
            
            print(system.time(m<-tryCatch.W.E(gam(formu,data=ddd,method=method,family=nb(),knots=knots))$value));

            if(class(m)[2] == "error") {
                print(m)
                stop("Error occured for length ", a, "\n")
            }
            
            ## Calculate total log-likelihood
            totll = logLik(m)
            
            if(is.null(predD)) stop("predD is null") ##predD=subset(ddd,haul.id %in% myids);
            res=numeric(length(pred.ctimes));
            lores=res;
            upres=res;
            gp2=list()

            predD$TimeShotHour=mean(ddd$TimeShotHour)
            predD$Ship=names(which.max(summary(ddd$Ship)))
            predD$timeOfYear=mean(ddd$timeOfYear);
            predD$HaulDur=30.0
            predD$Gear=myGear;
            predD$dum=0;
            
            for(y in 1:length(pred.ctimes)){ 
                ## take care of years with all zeroes
                ##if(all(ddd$y[ddd$Year==y]==0)){
                ##    res[which(as.character(yearRange)==y)]=0;
                ##    upre*s[which(as.character(yearRange)==y)] = 0;
                ##    lores[which(as.character(yearRange)==y)] = 0;
                ##    next;
                ##}
                ## OBS: effects that should be removed should be included here
                predD$Year=pred.years[y]; 
                predD$ctime=pred.ctimes[y] ##as.numeric(as.character(y));
                predD$Quarter=pred.quarters[y]
                
                p=try(predict(m,newdata=predD,newdata.guaranteed=TRUE));
                
                ## take care of failing predictions
                ##if(!is.numeric(p)) {
                ##    res[which(as.character(yearRange)==y)]=0;
                ##    upres[which(as.character(yearRange)==y)] = 0;
                ##    lores[which(as.character(yearRange)==y)] = 0;
                ##    next;
                ##}

                res[y] = sum(exp(p));
                gPred=p
                gp2[[y]]=gPred;
                if(nBoot>10){
                    Xp.1=predict(m,newdata=predD,type="lpmatrix");
                    library(MASS)
                    brp.1=mvrnorm(n=nBoot,coef(m),m$Vp);
                    OS.pos = matrix(0,nrow(predD),nBoot);
                    terms.pos=terms(m)
                    if(!is.null(m$offset)){
                        off.num.pos <- attr(terms.pos, "offset")
                        for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos, 
                                                                            "variables")[[i + 1]], predD)
                    }
                    rep1=exp(Xp.1%*%t(brp.1)+OS.pos);
                    
                    idxSamp = colSums(rep1);
                    upres[y] = quantile(idxSamp,0.975);
                    lores[y] = quantile(idxSamp,0.025);
                }
            } ## rof years
            list(res=res,m=m,lo=lores,up=upres,gp=gPred,ll=totll,gp2=gp2);
        }## end do.one
        noLengths=length(lengths);
        rr=mclapply(1:noLengths,do.one.a,mc.cores=mc.cores);
        logl=0;
        for(a in 1:noLengths){
            resMat[,a]=rr[[a]]$res;
            models[[a]]=rr[[a]]$m;
            loMat[,a]=rr[[a]]$lo;
            upMat[,a]=rr[[a]]$up;
            gPreds[[a]]=rr[[a]]$gp;
            logl=logl+rr[[a]]$ll
            gPreds2[[a]]=rr[[a]]$gp2
        }
        getEdf<-function(m) sum(m$edf)
        totEdf=sum( unlist( lapply(models,getEdf)))
        list(idx=resMat,models=models,lo=loMat,up=upMat,gPreds=gPreds,logLik=logl,edfs=totEdf,gPreds2=gPreds2);
    }
