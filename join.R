files <- dir("post", full=TRUE)

getEst <- function(x, what) {
    load(x)
    eval(parse(text=what))
}

beta <- sapply(files, getEst, "beta")
betaSD <- sapply(files, getEst, "betaSD")

logindex <- sapply(files, getEst, "logindex[,1]")
logindexSD <- sapply(files, getEst, "logindex[,2]")

## Gear
i <- grep("Gear", rownames(beta))
matplot( t(beta[i, ] ), type="b")
matplot( t(betaSD[i, ] ), type="b")

## Top punkt
tpRS <- -.5 * beta["Quarter1:DepthRS",] / beta["Quarter1:I(DepthRS^2)",]
tp1 <- (tpRS * 20 + 50)

tpRS <- -.5 * beta["Quarter4:DepthRS",] / beta["Quarter4:I(DepthRS^2)",]
tp4 <- (tpRS * 20 + 50)

matplot(cbind(tp1,tp4), type="b")
