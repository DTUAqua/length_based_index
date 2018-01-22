files <- dir("post", full=TRUE)

getEst <- function(x, what) {
    load(x)
    eval(parse(text=what))
}

## All fixed effect estimates (Gear, timeofday, Depth) and uncertainties
beta <- sapply(files, getEst, "beta")
betaSD <- sapply(files, getEst, "betaSD")

## Log index and uncertainties
logindex <- sapply(files, getEst, "logindex[,1]")
logindexSD <- sapply(files, getEst, "logindex[,2]")

## Example: Gear
i <- grep("Gear", rownames(beta))
matplot( t(beta[i, ] ), type="b")
matplot( t(betaSD[i, ] ), type="b")

## Parabola vertex
tpRS <- -.5 * beta["Quarter1:DepthRS",] / beta["Quarter1:I(DepthRS^2)",]
tp1 <- (tpRS * 20 + 50)

tpRS <- -.5 * beta["Quarter4:DepthRS",] / beta["Quarter4:I(DepthRS^2)",]
tp4 <- (tpRS * 20 + 50)

matplot(cbind(tp1,tp4), type="b")
