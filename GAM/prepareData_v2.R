library(DATRAS)
library(sp)
library(rgdal)
library(maps)
library(mapdata)

datafile="../EBcod.RData"

if(!file.exists(datafile)){
    
    sti="~/Documents/DATRAS";
    years=1991:2018
    genus="Gadus"
    bfamily="morhua";
    
    dAll <- readExchangeDir(paste(sti,"/exchange/BITS/",sep=""),strict=FALSE)
    
    dAll<-addSpatialData(dAll,"~/Documents/shapefiles/ICES_areas.shp")
    
    dAll<-subset(dAll,Species==paste(genus,bfamily),Year %in% years,HaulVal %in% c("V","A","C","N"),StdSpecRecCode%in%c(1,3))
    
    
    mytab=xtabs(~Gear,dAll[["HH"]])
    
    goodGears=names(mytab[mytab>120])
    
    dAll<-subset(dAll, Gear %in% goodGears)
    
    save(dAll,file=datafile)
} else load(datafile)

LongLatToUTM<-function(x,y,zone){
 xy <- data.frame(ID = 1:length(x), X = x, Y = y)
 coordinates(xy) <- c("X", "Y")
 proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
 res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
 return(as.data.frame(res))
}

long2UTMzone <- function(long) {
    (floor((long + 180)/6) %% 60) + 1
}

## convert to utm spatial coordinates
myzone = long2UTMzone( median(dAll$lon) )
utmcoords = LongLatToUTM(dAll$lon, dAll$lat,myzone)
dAll[[2]]$utm.x = utmcoords$X
dAll[[2]]$utm.y = utmcoords$Y

## add spectrum (OBS last group is NOT a plus group, do we want that instead?)
##dAll<-addSpectrum(dAll,cm.breaks=seq(5,100,by=1))
## Margit wants: 5 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80 85 90 95 100 105 110 115 120
##dAll=addSpectrum(dAll,cm.breaks=c(0,10,14,16,18,20:41,43,45,47,50,55,61,200))
dAll=addSpectrum(dAll,cm.breaks=c(0,10,12,14,16,18,20:41,42,44,46,48,50,52,54,56,58,60,200))


## make "continuous" time variable, discretized to quarterly to prevent too much time wiggliness
dAll[[2]]$ctime = as.numeric( as.factor(as.numeric(dAll$Year)+(as.numeric(as.character(dAll$Quarter) )-1)/4)   ) ##dAll$timeOfYear
##dAll[[2]]$ctime = as.numeric( cut( dAll$ctime, seq( floor(min(dAll$ctime)), ceiling(max(dAll$ctime)), by=1/4)))


dQ14 = subset(dAll, Quarter %in% c("1","4"))

hist(dQ14$HaulDur[dQ14$HaulVal=="N"])
dQ14$HaulVal[is.na(dQ14$HaulDur)]

## some "No oxygen" hauls have not been performed (HaulDur=0 or NA)
## Set HaulDur to a positive number for these such that they are included
dQ14$HaulDur[ dQ14$HaulDur<5 | is.na(dQ14$HaulDur) ] = 10


## Get depth database :Download link
## http://data.bshc.pro/ogc/bsbd-0.9.3?SERVICE=wcs&VERSION=1.0.0&REQUEST=GetCoverage&coverage=bsbd&CRS=EPSG:4326&bbox=7.274582044420753,53.09571677875055,24.780872897306068,58.88722383749357&width=2174&height=1289&format=XYZ
bathy <- read.table("~/Desktop/bathy/bsbd-0.9.3.xyz")
bathy$V3[ bathy$V3 > 0 ] = 0 ## not interested in land

depthLookup<-function(lon,lat, bathy){ 
   dLon = diff(as.numeric(rownames(bathy)[1:2]))
   dLat = diff(as.numeric(colnames(bathy)[1:2]))
   lon0 =  as.numeric(rownames(bathy)[1])
   lat0 =  as.numeric(colnames(bathy)[1])
   ## lon = lon0 + (row-1)*dLon
   row = round(lon/dLon - lon0/dLon + 1)
   col = round(lat/dLat - lat0/dLat + 1)
   if(row>0 && col>0 && row <= nrow(bathy) && col<=ncol(bathy))
       return(bathy[row,col])
   return(NA)
}

## make it into a matrix
bathymatrix <- xtabs( V3 ~ V1 + V2, data=bathy,addNA=TRUE)

## look up depth by coordinates
dQ14$depth2 = sapply( 1:length(dQ14$Depth), function(i) -depthLookup(dQ14$lon[i],dQ14$lat[i],bathymatrix) )

## replace Depth NA's with lookup
sel=which(is.na(dQ14$Depth) )
dQ14$Depth[sel] = dQ14$depth2[sel] 
depthModel=lm( Depth ~ depth2, data=dQ14[[2]])

## start with 26th cm group (most abundant)
##dQ14$y = dQ14$N[,CMGROUP]

knots=list(TimeShotHour=seq(0,24,length=6))

## include a bit of kattegat to improve edge estimates
dQ14 = subset(dQ14,!ICES_SUB %in% c(20,21) | lat < 56.5) 


## Add ICES areas
names(bathy)<-c("lon","lat","Depth")
require(maptools)
shp<-readShapeSpatial("~/Documents/shapefiles/ICES_areas.shp")
coordinates(bathy) <- ~lon + lat
xtra <- over(bathy, shp)
bathy$ICES_SUB = xtra$ICES_SUB

## Extract EBcod area grid with depths (25:32)
EBarea<-subset(bathy, ICES_SUB %in% as.character(22:32))

gridutm = LongLatToUTM(EBarea$lon,EBarea$lat,zone=myzone)
EBarea$utm.x = gridutm$X
EBarea$utm.y = gridutm$Y
EBarea=as.data.frame(EBarea)

my.palette<-colorRampPalette(c("darkblue","mediumblue","lightblue1"))
my.palette.vec=my.palette(100);

EBarea$ctime=median(dQ14$ctime)
EBarea$HaulDur=30
EBarea$TimeShotHour=median(dQ14$TimeShotHour)
EBarea$Quarter="1"
EBarea$Gear=tail(names(sort(table(dQ14$Gear))),1)
EBarea$Depth = - EBarea$Depth
EBarea=subset(EBarea,Depth>5)

grid.size.m = 2500
EBarea$cututmx = cut(EBarea$utm.x,seq(min(EBarea$utm.x)-1,max(EBarea$utm.x)+grid.size.m,by=grid.size.m))
EBarea$cututmy = cut(EBarea$utm.y,seq(min(EBarea$utm.y)-1,max(EBarea$utm.y)+grid.size.m,by=grid.size.m))
EBarea$gridcell = factor(paste(EBarea$cututmy,EBarea$cututmx,sep=":"))

## reduce grid resolution to speed up predictions
meanDepth=aggregate(Depth ~ gridcell, data=EBarea,FUN=mean)
EBarea$meanDepth = meanDepth$Depth[ EBarea$gridcell ]
meanutmx = aggregate( utm.x ~ gridcell, data=EBarea,FUN=mean)
meanutmy = aggregate( utm.y ~ gridcell, data=EBarea,FUN=mean)
EBarea$meanutmx = meanutmx$utm.x[ EBarea$gridcell ]
EBarea$meanutmy = meanutmy$utm.y[ EBarea$gridcell ]
EBarea.s = EBarea[!duplicated(EBarea$gridcell),]

## replace depth with corrected depth
EBarea.s$Depth = predict( depthModel, newdata=data.frame(depth2=EBarea.s$meanDepth))
EBarea.s$utm.x = EBarea.s$meanutmx
EBarea.s$utm.y = EBarea.s$meanutmy

## cut out north eastern most samples (very few)
plot(dQ14$lon,dQ14$lat)
abline(h=58.15,col=2)
abline(v=22,col=2)
dQ14 = subset( dQ14, lon < 22 & lat < 58.15)


## Cut out north eastern areas with no samples
EBarea.s = subset( EBarea.s, utm.x < max(dQ14$utm.x) )
EBarea.s = subset( EBarea.s, utm.y < max(dQ14$utm.y) )

## Almost no data from 2017 Q4, omit
xtabs(~ Year + Quarter, dQ14[[2]])
##dQ14 = subset( dQ14, !( Year=="2017" & Quarter=="4"))

DepthQs = quantile(dQ14$Depth,probs=c(0.005,0.995))

coastlines = map('worldHires',xlim=range(EBarea.s$lon),ylim=range(EBarea.s$lat),plot=FALSE) 
nonNA = !is.na(coastlines$x)
coastlinesproj = LongLatToUTM(coastlines$x[nonNA],coastlines$y[nonNA],myzone)
coastlines$x[nonNA] = coastlinesproj$X[nonNA]
coastlines$y[nonNA] = coastlinesproj$Y[nonNA]

## Cut out depths outside measured range [ No... instead force parabolic depth effect ]
##EBarea.s = subset(EBarea.s, Depth>=DepthQs[1] & Depth<=DepthQs[2])

dQ14$timeOfDay = 1
dQ14$timeOfDay[  dQ14$TimeShotHour<7 | dQ14$TimeShotHour>17  ] = 2
dQ14$timeOfDay[ dQ14$TimeShotHour>10 & dQ14$TimeShotHour<14  ] = 3

##WBEBarea.s<-EBarea.s
##EBarea.s<-subset(EBarea.s, ICES_SUB %in% as.character(25:32))
pvec=array(NA,dim=c(ncol(dQ14$N),nrow(EBarea.s),length(unique(dQ14$ctime))))
pvecs = list()

## EB cod area 25+
pvecEB.1 = as.numeric(EBarea.s$ICES_SUB %in% as.character(25:32))
for(lg in 1:ncol(dQ14$N)){
    for(ct in 1:length(unique(dQ14$ctime))){
        pvec[lg,,ct] = pvecEB.1
    }
}
pvecs[[1]]<-pvec
### EB cod 13 degrees
pvecEB.2 = as.numeric(EBarea.s$lon>13)
for(lg in 1:ncol(dQ14$N)){
    for(ct in 1:length(unique(dQ14$ctime))){
        pvec[lg,,ct] = pvecEB.2
    }
}
pvecs[[2]]<-pvec
## WB cod (13 degrees)
pvecEB.3 = as.numeric(EBarea.s$lon<=13)
for(lg in 1:ncol(dQ14$N)){
    for(ct in 1:length(unique(dQ14$ctime))){
        pvec[lg,,ct] = pvecEB.3
    }
}
pvecs[[3]]<-pvec

#### WB cod (soft)
load("~/Documents/DTUAqua/length_based_index/WBcod/alldenssplit-10.RData")
pvec=array(NA,dim=c(ncol(dQ14$N),nrow(EBarea.s),length(unique(dQ14$ctime))))
for(lg in 1:ncol(dQ14$N)){
    for(ct in 1:length(unique(dQ14$ctime))){
        pvec[lg,,ct] = alldens.split[,ct,lg]
    }
}
pvecs[[4]]<-pvec



## WB cod (penkowa)
sel = which(as.numeric(as.character(EBarea.s$ICES_SUB))>24)
sel2 = which(as.numeric(as.character(EBarea.s$ICES_SUB))<24)
for(lg in 1:ncol(dQ14$N)){
    for(ct in 1:length(unique(dQ14$ctime))){
        pvec[lg,sel,ct] = 0
        pvec[lg,sel2,ct] = 1
    }
}
pvecs[[5]]<-pvec

## EB cod (penkowa)
pvecs[[6]]<-1-pvec

save(dQ14,EBarea.s,pvecs,coastlines,file="../EBcodProcessedData.RData")

###
my.palette<-colorRampPalette(c("darkblue","mediumblue","lightblue1"))
my.palette.vec=my.palette(100);
plot(EBarea.s$utm.x,EBarea.s$utm.y,col=rev(my.palette.vec)[cut(EBarea.s$Depth,100)],pch=16,cex=0.5)

plot(EBarea.s$utm.x,EBarea.s$utm.y,col=rev(my.palette.vec)[cut(pvecs[[5]][10,,10],100)],pch=16,cex=0.5)


##plot(WBEBarea.s$utm.x,WBEBarea.s$utm.y,col=rev(my.palette.vec)[cut(WBEBarea.s$Depth,100)],pch=16,cex=0.5)

###
