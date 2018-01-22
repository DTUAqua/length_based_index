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
mat <- t(beta[i, ] )
matplot( mat, type="b", main="Gear effect (relative to TVL)", xlab="Size group", ylab="log(Effect)")
legend("bottomright", colnames(mat), pch=c(1L:9L, 0L, letters, LETTERS),  col=1:6)

mat <- t(betaSD[i, ] )
matplot( mat, type="b", main="Gear effect Std. Dev.", xlab="Size group", ylab="SD(log(Effect))")
legend("topright", colnames(mat), pch=c(1L:9L, 0L, letters, LETTERS),  col=1:6)


## Parabola vertex
tpRS <- -.5 * beta["Quarter1:DepthRS",] / beta["Quarter1:I(DepthRS^2)",]
tp1 <- (tpRS * 20 + 50)
tpRS <- -.5 * beta["Quarter4:DepthRS",] / beta["Quarter4:I(DepthRS^2)",]
tp4 <- (tpRS * 20 + 50)

## Parabola 'sd'
sdRS <- 1 / sqrt(-2 * beta["Quarter1:I(DepthRS^2)",])
sd1 <- (sdRS * 20)
sdRS <- 1 / sqrt(-2 * beta["Quarter4:I(DepthRS^2)",])
sd4 <- (sdRS * 20)

mat <- cbind(tp1,tp4)
matsd <- cbind(sd1, sd4)

## Add sds to preferred depth:
## mat <- cbind(mat, mat-matsd, mat+matsd)

matplot(mat, type="b", main="Preferred depth", xlab="Size group", ylab="Depth (m)", col=1:2)
legend("topleft", c("Quarter 1", "Quarter 4"), pch=c(1L:9L, 0L, letters, LETTERS),  col=1:2)
