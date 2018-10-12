fitALK<-function (x, minAge, maxAge, mc.cores = 1, model = c("cra~LngtCm", 
    "cra~poly(LngtCm,2)", "cra~LngtCm+s(lon,lat,bs='ts')")[method], 
    method = 1, autoChooseK = FALSE, useBIC = FALSE, varCof = FALSE, 
    maxK = 100, gamma = 1.4, verbose = FALSE, ...) 
{
    checkSpectrum(x)
    if ((minAge + 1) > maxAge) 
        stop("Invalid age selection.")
    ages = minAge:maxAge
    nAges = length(ages)
    lastAge = ages[nAges]
    ages = ages[-nAges]
    extraVars = unlist(lapply(lapply(model, as.formula), xtraVars, 
        x = x))
    if (length(extraVars) == 0) 
        extraVars = NULL
    x[[1]] = merge(x[[1]], x[[2]][c("lon", "lat", "haul.id", 
        extraVars)], by = "haul.id", all.x = TRUE, sort = FALSE, 
        suffixes = c("", ".y"))
    x[[1]] = subset(x[[1]], !is.na(Year) & !is.na(Age) & !is.na(LngtCm))
    mylapply <- function(...) {
        hasmc = (mc.cores > 1 && require(multicore, quietly = TRUE))
        if (!hasmc) 
            return(lapply(...))
        else return(mclapply(..., mc.cores = mc.cores))
    }
    if (verbose) 
        cat("Fitting model...")
    models = mylapply(ages, fitALKone, ages = ages, AL = x[[1]], 
        model = model, gamma = gamma, autoChooseK = autoChooseK, 
        useBIC = useBIC, varCof = varCof, maxK = maxK, verbose = verbose, 
        ...)
    class(models) <- "ALKmodel"
    attr(models, "data") <- x
    attr(models, "ALKformula") <- model
    attr(models, "ages") <- ages
    models
}

fitALKone<-function (a, ages, AL, model, gamma, autoChooseK = FALSE, useBIC = FALSE, 
    varCof = FALSE, maxK = 100, verbose = FALSE, splineType = 1, ...) 
{
    ##cat("flaf\n");
    if (length(model) > 1) {
        idx = which(ages == a)
        f <- as.formula(model[idx])
    }
    else {
        f <- as.formula(model)
    }
    require(mgcv, quietly = TRUE)
    myd = subset(AL, Age >= a)
    myd$cra = as.factor(myd$Age > a)
    if (autoChooseK) {
        uniqueCovs = length(unique(paste(myd$lon, myd$lat)))
        if (!varCof) {
            k = min(maxK, uniqueCovs - 1)
            f = as.formula(paste("cra~LngtCm+s(lon,lat,k=", k, 
                ",bs='ts')"))
            if(splineType==2){
                khalf = floor(sqrt(k))
                f = as.formula(paste("cra~LngtCm+te(lon,lat,d=c(1,1),bs=c('ds','ds'),k=c(", khalf,",",khalf,
                "),m=list(c(1,0),c(1,0)))")) ## Duchon spline 1st order
                 
            }
            if (uniqueCovs < 10) 
                f = as.formula("cra~LngtCm")
        }
        else {
            k = min(maxK, uniqueCovs/2 - 1)
            f = as.formula(paste("cra~s(lon,lat,by=LngtCm,k=", 
                k, ",bs='ts')+s(lon,lat,k=", k, ",bs='ts')"))
            if(splineType==2) f = as.formula(paste("cra~LngtCm+s(lon,lat,k=", k, 
                ",bs=c('ds','ds'),m=list(c(1,0),c(1,0)),by=LngtCm)+s(lon,lat,k=", k, 
                ",bs=c('ds','ds'),m=list(c(1,0),c(1,0))) ")) ## Duchon spline 1st order 
            if (uniqueCovs < 10) 
                f = as.formula("cra~LngtCm")
        }
        if (useBIC) 
            gamma = log(sum(myd$NoAtALK))/2
    }
    ##print(f);
    ##browser()
    m <- tryCatch.W.E(gam(f, data = myd, family = "binomial", 
        weights = NoAtALK, gamma = gamma, ...))$value
    if (class(m)[2] == "error") {
        print(m)
        print(summary(myd))
        stop("Error occured for age ", a, "\\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\\n")
    }
    if (verbose) {
        print(summary(m))
    }
    return(m)
}
