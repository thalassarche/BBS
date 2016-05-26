######################
### INITIAL SET UP ###
######################

# Clear workspace
rm(list = ls())

# Load packages
require(XLConnect)   # for reading Excel files
require(dplyr)       # for joins
require(coda)
require(arm)
require(lme4)        # For mixed effect models

########################
### PREPARE DATA     ###
########################

setwd("~/Dropbox/FILE_BACKUP/Projects/OregonForestsBBS/data/BBS_data")
CA <- read.csv("Califor.csv", header=TRUE)
OR <- read.csv("Oregon.csv", header=TRUE)
WA <- read.csv("Washing.csv", header=TRUE)
weather <- read.csv("weather.csv", header=TRUE)
nwfp_routes <- read.csv("BBS_Request_Sauer/nwfp_routes.csv", header=TRUE)
rd0 <- readWorksheetFromFile("custom/routes.xlsx", sheet=1, header=TRUE)

# define folder for model outputs
setwd("~/Dropbox/FILE_BACKUP/Projects/OregonForestsBBS/output")

# combine the three relevant states
CAORWA <- rbind(CA, OR, WA)

# add the zeroes for species not recorded on a route
roy <- CAORWA[,c("countrynum", "statenum", "Route", "RPID", "Year")]
uroy <- unique(roy[c("countrynum", "statenum", "Route", "RPID", "Year")]) # all route-year combinations counted
uroy$id <- seq(1,nrow(uroy),1)
spp <- CAORWA[,"Aou"]
uspp <- unique(spp)                                # all species counted in these states
allcomb <- expand.grid(seq(1,nrow(uroy),1), uspp)  # table with a key plus each species
names(allcomb) <- c("id", "Aou")
routespp <- left_join(allcomb, uroy)               # joins on id field
datz <- left_join(routespp, CAORWA)
datz[is.na(datz)] <- 0                             # fill in zeroes for species not recorded on transects

# merge with weather data (automagically works out what variables to use)
dat1 <- left_join(datz, weather)

# make Route a three-digit number, and AOU species code a five digit number
dat1$Route <- formatC(dat1$Route, width = 3, format = "d", flag = "0")
dat1$Aou <- formatC(dat1$Aou, width = 5, format = "d", flag = "0")

# add new column with StateRoute, then that plus Observer ID
dat1$StateRoute <- as.numeric(paste(dat1$statenum, dat1$Route, sep=""))
dat1$RouteObs <- paste(dat1$StateRoute, dat1$ObsN, sep="_")

# extract data for the NWFP area transects, and where RunType=1, and where (weather) data complete
# RunType "1" is acceptable data
# RunType "0" is unacceptable data for USGS trend analyses purposes
dat <- left_join(nwfp_routes, dat1)
dat <- dat[dat$RunType==1,]
dat <- dat[complete.cases(dat),]   # Only removes three rows

# add data on ownership
dat <- left_join(dat, rd0)

# define single breakpoint if using two segments (1990/1991, i.e. 1991 splits data evenly)
breakp <- 1991

# select and extract only the useful columns
seldat <- dat[,c("StateRoute","Year", "period", "Aou", "SpeciesTotal", "ObsN", "RouteObs", "RouteName",
                 "Elevation", "Latitude", "Longitude", "ForestProp", "TreeCover",
                 "PropFederal", "PropIndustry", "PropPrivate", "HasFederal", "HasPrivate")]

# Extract data for each focal species 
# see: Phalan_SpeciesTimesRoutes.xlsx
focalsp <- "04330"
RUHU1 <- seldat[seldat$Aou==focalsp,]

# Function to remove transects with all zero counts for a species
# and single counts for an observer-route combination (of which, v few)
# This second part disabled by setting to 0 (not 1) in calculation of c
removezero <- function(d){
  a <- aggregate(d$SpeciesTotal, list(d$StateRoute), sum)
  names(a) <- c("StateRoute", "AllYearTotal")
  updated <- left_join(d, a)
  b <- updated[updated$AllYearTotal>0,]
  routesums2 <- aggregate(b$SpeciesTotal, list(b$StateRoute, b$ObsN), length)
  names(routesums2) <- c("StateRoute", "ObsN", "ObsRouteCounts")
  updated2 <- left_join(b, routesums2)
  c <- updated2[updated2$ObsRouteCounts>0,]
  return(c)
}

# Select years
# e.g., remove 2014, to compare with BBS online route summaries
RUHU <- RUHU1[RUHU1$Year<breakp,]          
dat <- removezero(RUHU)

# generate list of routes
routelist <- unique(dat$StateRoute)

###############################
### RUN MODEL FOR ONE ROUTE ###
###############################

m1 = NULL
rdat <- dat[dat$StateRoute==69047,]  # Select a route using routelist[i] or StateRoute#
head(rdat)

# Define model components
nobservers <- length(unique(rdat$ObsN))
count <- rdat$SpeciesTotal
ncounts <- length(count)
year=rdat$Year-1967
nyears <- length(year)
fixedyear <- round(nyears/2-1)  # In BBS code, fixedyear is 20
obser=as.numeric(factor(rdat$ObsN))
nchains <- 10
year2 <- year-fixedyear

# Prepare the data
jags.dat <- list(ncounts=ncounts, nobservers=nobservers, nyears=nyears,
                 fixedyear=fixedyear, count=count, year2=year2, obser=obser)

# Run model as Poisson with log link (or, first, negative binomial)
#m1 <- glmer.nb(count ~ year2 + (1|obser) -1, data=jags.dat)
m1 <- glmer(count ~ year2 + (1|obser) - 1, data=jags.dat, poisson(link = "log"))
summary(m1)

# Predictions from model
mean.obs <- mean(unlist(ranef(m1)$obser))
beta <- fixef(m1)                      # coefficient for year
sd.obs <- sd(unlist(ranef(m1)$obser))
fitted <- rep(NA, length(nyears))
for (t in 1:nyears) {
  fitted[t] <- exp(mean.obs + beta*year2[t] + 0.5*sd.obs*sd.obs)
  fitted <- fitted
}

#fitted <- exp(predict(m1, newdata=NULL, re.form=NA)) # predictions
B <- 100*((fitted[nyears]/fitted[1])^(1/(nyears-1))-1)
B

#################################################
### DEFINE AND RUN MODELS FOR MULTIPLE ROUTES ###
#################################################

##### Define functions #####

# Function for when there is one observer
singleobs <- function(jags.dat, year2, nyears, year, focalsp, routelist, nobservers, count, out) {
  m1 <- NULL   # reset m1
  try(m1 <- glm(count ~ year2, data=jags.dat, poisson(link = "log")))
  
  alpha <- as.numeric(coef(m1)[1])
  beta <- as.numeric(coef(m1)[2])                           # coefficient for year
  sd.yr <- coef(summary(m1))[2,2] * sqrt(length(year2))     # standard deviation of year (may be inaccurate)
  vari <- sd.yr*sd.yr                                       # variance
  fitted <- rep(0,length(year2))
  for (t in 1:nyears) {
    fitted[t] <- exp(alpha + beta*year2[t])  # Predictions from model
  }
  B <- 100*((fitted[nyears]/fitted[1])^(1/(nyears-1))-1)            # Annual % trend
  
  out[i,] <<- c(focalsp, routelist[i], min(year), max(year), nyears, nobservers, mean(count), B, vari)
}

# Function for when there are multiple observers
multiobs <- function(jags.dat, year2, nyears, year, focalsp, routelist, nobservers, count, out) {
  m1 <- NULL   # reset m1
  try(m1 <- glmer(count ~ year2 + (1|obser) - 1, data=jags.dat, poisson(link = "log")))
  
  mean.obs <- mean(unlist(ranef(m1)$obser))
  beta <- fixef(m1)                               # coefficient for year
  sd.obs <- as.numeric(as.data.frame(VarCorr(m1))[5])         # not: sd(unlist(ranef(m1)$obser))
  sd.yr <- se.fixef(m1) * sqrt(length(year2))     # standard deviation of year (may be inaccurate)
  vari <- sd.obs*sd.obs + sd.yr*sd.yr             # sum of variances
  fitted <- rep(0,length(year2))
  for (t in 1:nyears) {
    fitted[t] <- exp(mean.obs + beta*year2[t] + 0.5*sd.obs*sd.obs)  # Predictions from model
  }
  B <- 100*((fitted[nyears]/fitted[1])^(1/(nyears-1))-1)            # Annual % trend
  
  out[i,] <<- c(focalsp, routelist[i], min(year), max(year), nyears, nobservers, mean(count), B, vari)
}

# Set up dataframe to store results
N <- length(unique(dat$StateRoute))
cols <- c("species", "StateRoute", "startyear", "endyear", "nyears",
          "nobservers", "meancount", "B", "variance")
out <- data.frame(matrix(0, nrow=N, ncol=length(cols)))
colnames(out) <- cols

# Record the time
time1 <- Sys.time()

# Assumes model already defined from previous section
for (i in 1:length(routelist)) {
  
  # Select route
  temp <- dat[dat$StateRoute==routelist[i],]  # Select a route
  
  # Define model components
  nobservers <- length(unique(temp$ObsN))
  count <- temp$SpeciesTotal ;   ncounts <- length(count)
  year=temp$Year-1967        ;   nyears <- length(year)
  fixedyear <- round(nyears/2-1)  # In BBS code, fixedyear is 20
  obser=as.numeric(factor(temp$ObsN))
  year2 <- year-fixedyear
  oneobs <- nobservers==1      # Check if one observer, or more
  
  # Prepare the data
  jags.dat <- list(ncounts=ncounts, nobservers=nobservers, nyears=nyears,
                   fixedyear=fixedyear, count=count, year=year, obser=obser)
  
  ifelse (oneobs=="TRUE", try(singleobs(jags.dat, year2, nyears, year, focalsp, routelist, nobservers, count, out)),
          try(multiobs(jags.dat, year2, nyears, year, focalsp, routelist, nobservers, count, out)))
}

# Output the time taken
Sys.time()-time1

# Coerce back to numeric (the NAs coerce it to character)
for (i in 1:ncol(out)) {
  out[,i] <- as.numeric(out[,i])
}

# Subset routes with meancount of at least 0.25 (1 bird/four years)
# and where average observer contributes at least four years of data
out_clean <- subset(out, meancount>=0.25 & (nyears/nobservers)>=4)
# It may be better to produce a column for each of these and then select/weight using those

# Save out table (or out_clean)
write.table(out,file="RUHU_trends_poisson.csv", sep=",")
#out <- read.csv("RUHU_trends_poisson.csv", sep=",")

#######################
### EXAMINE RESULTS ###
#######################

par(mfrow=c(1,1))

# Visualise data
plot(count ~ year2, col=obser)
x <- year2
y <- count
lines(x, fitted)

# Compare with trends from BBS
wdir <- getwd()
setwd("/Users/benphalan/Dropbox/FILE_BACKUP/Projects/OregonForestsBBS/misc/Analyses/2016_01_28/data")
bbs0 <- readWorksheetFromFile("SppBBS_Trends04330RUHU_20151203.xlsx", sheet=1, header=TRUE)
setwd(wdir)

out$StateRoute <- as.numeric(out$StateRoute)
RUHU.cf <- left_join(out_clean, bbs0, by = c("StateRoute"="RouteNum"))
RUHU.cf <- data.frame(RUHU.cf)
RUHU.cf[,c(2,3,4,5,6,7,8,12,21,23,25)]   # Need to update

x <- as.numeric(RUHU.cf$B)
y <- as.numeric(RUHU.cf$EEQ_TrendEst)
x.sd <- as.numeric(RUHU.cf$variance)
y.sd <- sqrt(as.numeric(RUHU.cf$EEQ_Variance))
plot(x,y, xlab="B", ylab="EEQ Trend")#, xlim=c(-30,30), ylim=c(-30,75))
abline(0,1, lty=3, lwd=2)
abline(v=0, lty=3, lwd=1)
abline(h=0, lty=3, lwd=1)
arrows(x, y-y.sd, x, y+y.sd, length=0.05, angle=90, code=3, col="grey")
arrows(x-x.sd, y, x+x.sd, y, length=0.05, angle=90, code=3, col="grey")
points(x,y)
abline(lm(y~x, na.action=na.exclude))    # plot regression line
summary(lm(y~x, na.action=na.exclude))   # check R2 of trend comparison

# Data visualising
d <- x-y
plot(d ~ RUHU.cf$meancount)
require(ggplot2)
cs <- RUHU.cf$meancount
bs <- RUHU.cf$B
ratioyo <- RUHU.cf$nyears/RUHU.cf$nobservers
p <- qplot(ratioyo, bs, colour=cs)
p + scale_colour_gradient(limits=c(0,0.25))

### FOOTNOTES
# BBS model from: ftp://ftpext.usgs.gov/pub/er/md/laurel/Sauer/Population%20Change%20Estimation%20Webinar/
# JAGS e.g. adapted from POLS 506 Bonus lecture here: http://jee3.web.rice.edu/teaching.htm

