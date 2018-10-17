#--Create package------------------------------------------------------------------

if(F){
  #Packages
  .rs.restartR()
  packages  <- c("devtools", "roxygen2")
  if(length(packages[!(packages %in% installed.packages()[,"Package"])])){
    install.packages(packages[!(packages %in% installed.packages()[,"Package"])])}
  suppressWarnings(loaded <- lapply(packages, require, character.only = TRUE)); names(loaded) <- packages; loaded
  
  #PackageName
  pn <- "VerDi"
  
  #DescriptionFile
  custdes <- list("Package" = "VerDi",
                  "Title" = "Modelling the vertical distribution of aquatic organisms",
                  "Maintainer" = "Robert Herrmann <robertherrmann@mail.de>",
                  "Author" = "Robert Herrmann",
                  "Version" = "0.1",
                  "Description" =  "Predicting the vertical distribution of pelagic and semi-pelagic fish and zooplankton. Model approach is based on Evaluation Functions (rule-based systems). Required knowledge concerns the species' habitat preferences represented as environmental factors, i.e abiotic or biotic parameter thresholds (e.g. water temperature, irradiation, food availability, etc.). Further functions (e.g. modelling underwater light regimes based on remote sensing) are included...",
                  "Depends"= "R (>= 3.5.1)",
                  "License" = "What license is it under?",
                  "Encoding" = "UTF-8",
                  "LazyData" = "true",
                  "RoxygenNote" = "6.1.0",
                  "Import" = "data.table")
  
  #CreatePackage
  setwd(paste(Sys.getenv("HOME"), "Dropbox/GitHub", sep = "/"))
  create(pn, custdes)
  
  #CreateDocumentation
  setwd(paste("./", pn, sep = ""))
  getwd()
  document()
  
  #TestingPackage
  setwd("..")
  install(pn)
  library(VerDi)
  help(package = as.character(pn))
  
  #WritingManual
  check(pn, check_dir = getwd(), manual = TRUE)
  
  #LoadFromGitHub
  .rs.restartR()
  remove.packages(pn)
  install_github(paste("herrmannrobert", pn, sep = "/"))
  library(VerDi)
  help(package = as.character(pn))
}

#----------------------------------------------------------------------------------
#--Demographic noise including different options-----------------------------------
#----------------------------------------------------------------------------------

#' @title Demographic noise
#' @description Function generates random walk (demographic noise) based on a normal distribution (ND) 
#' or uniform distribution (UD).
#' @param par numeric vector or R object.
#' @param hs numeric, half saturation parameter, Default: 0.1
#' @param mv numeric, mean when op = "ND", Default: 0
#' @param sd numeric, standard deviation when op = "ND", Default: 0.1
#' @param min numeric, minimum value when op = "UD", Default: -0.2
#' @param max numeric, maximum value when op = "UD", Default: 0.2
#' @param op character, either 'ND' (normal distribution) or 'UD' (uniform distribution), Default: 'ND'
#' @details Is used to add demograhic noise (random walk) to e.g. the population's vertical distribution. 
#' Further details in rnorm() and runif().
#' @export
#' @examples 
#' \dontrun{
#' 
#' #Example shows difference between options (op)
#' #Variables (hs, mv, sd, min, max) are set to default (see description above)
#' 
#' rw <- data.frame()
#' 
#' set.seed(123)
#' for(i in 1:100){
#'  rw[i,1] <- RWalk(0)
#'  rw[i,2] <- RWalk(0, op = "UD")
#'  
#'  if(i == 100){colnames(rw) <- c("Normally distributed random walk", "Uniformly distributed random walk")
#'  
#'   par(mfrow = c(2, 2))
#'   plot(rw[,1], type = "l", main = colnames(rw)[1], ylab = "Random walk", xlab = "", ylim = c(-0.4, 0.4))
#'   plot(rw[,2], type = "l", main = colnames(rw)[2], ylab = "", xlab = "", ylim = c(-0.4, 0.4))
#'   hist(rw[,1], breaks = 100, main = "", xlab = "", xlim = c(-0.4, 0.4), ylim = c(0, 14))
#'   lines(density(rw[,1], adjust=2), lwd = 2)
#'   hist(rw[,2], breaks = 100, main = "", ylab = "", xlab = "", xlim = c(-0.4, 0.4), ylim = c(0, 14))
#'   lines(density(rw[,2], adjust=2), lwd = 2)
#'   }
#'   
#'  }
#' }
#' @export

RWalk <- function(par, hs = 0.1, mv = 0, sd = 0.1, min = -0.2, max = 0.2, op = "ND"){
  
  if(op == "ND"){
    (rnorm(1, mv, sd) + par) * abs(rnorm(1, mv, sd) + par) / (hs + (rnorm(1, mv, sd) + par)^2)
  }
  else{
    if(op == "UD"){
      (runif(1, min, max) + par) * abs(runif(1, min, max) + par) / (hs + (runif(1, min, max) + par)^2)
    }
  }
} 

#----------------------------------------------------------------------------------
#---------Evaluation function------------------------------------------------------
#----------------------------------------------------------------------------------

#' @title Evaluation function
#' @description Function computes the organism's vertical movement based on a set of single environmental factors. 
#' @param hydro      data frame including set of environmental parameters.
#' @param faclist    data frame specifying the species specific thresholds that induce a vertical movement.
#' @param nInd       numeric, number of individuals that are included in modelled environment.
#' @param mSpeed     numeric, maximum vertical distance that can be reached per time interval, Default: 1.
#' @param tStep      numeric, number of time increments, Default: 1.
#' @param dimList    numeric, intended return. See 'Details'.
#' @param rwalk      numeric vector indicating the demographic noise of the modelled population. A pre-defined random walk is included by default. See 'Details'.
#' @param rWalkOff   logical, if TRUE, demographic noise is turned off. If FLASE, a pre-defined random walk is included (default). See 'Details'.
#' 
#' @details The arguments hydro and faclist should be of class data.frame. Note that colnames of hydro and faclist 
#' have to be equal. Required layouts are given in 'Examples'.
#' 
#' @details If the argument nInd is provided, the model output includes the position of each individual for each time increament. The individuals are randomly 
#' distributed over the entire water column when the calculations begin. Note that computation time is strongly correlated with the amount of individuals included 
#' in the model environment.
#' 
#' @details Note that the arguments mSpeed and tStep are linked, i.e. mSpeed is the vertical distance each individual can reach in a time interval given 
#' by tStep. For example, if you aim to predict the vertical ditribution for an entire day with one-minute-intervals, than tStep is set to 1440 minutes, and 
#' consequently, mSpeed must be converted to meter per minute. For further details see 'Examples'.
#' 
#' @details The argument dimList provides three diffrent options of what specific data type is returned. Further details need to be written. See 'Examples'.
#'
#' @details Function makes use of RWalk. Default: c("ND", 0, 0.1, 0.1). If rWalk is set to NULL, the random walk is turned off and each individual would remain at 
#' the given water depth as long as the ambient environment stays sufficient enough.
#'
#' @details Further details in 'Examples'.
#' @examples 
#' \dontrun{
#'
#' #-1-----Example for static hydrography (stat.hyd) and given factor list (fl)
#' #-1.1---Generating hydrography commonly found in the Bornholm Basin (Baltic Sea)
#' set.seed(123456)
#' stat.hyd <- HydroGen(data.frame(Temp = c(17,13,4,4.8,8),
#'                                 Ox = c(6,5,3,2,0),
#'                                 Sal = c(7,8,9,18,20)), by = 25)
#'                                
#' #-1.2---Generating factor list for any kind of organism
#' fac  <- c("UL", "LL", "MOV", "RI")
#' temp <- c(  10,    8,    -1,  0.33)
#' ox   <- c(   1,  0.5,     1,  0.33)
#' sal  <- c(  10,   11,    -1,  0.33)
#' fl   <- data.frame(rbind(temp, ox, sal))
#' colnames(fl) <- fac
#' rownames(fl) <- colnames(stat.hyd)
#'
#' #-1.3.1-Executing EvalFunc(defaults) and calculating free vertival range of organisms
#' output.a <- EvalFunc(hydro = stat.hyd, faclist = fl)
#' FVR.a    <- which(output.a[,"PredMovByAllPar",] == 0)
#'
#' #-1.3.2-Alternatively use EvalFunc(dimList = 2) for extracting "PredMovByAllPar"
#' output.b <- EvalFunc(hydro = stat.hyd, faclist = fl, dimList = 2)
#' FVR.b    <- which(output.b$PredMovByAllPar == 0)
#'
#' #-1.4---Plotting
#' par(mfrow = c(1,1))
#' df <- as.data.frame(output.a)
#' plot(x = 0, ylim = c(nrow(df), 0), xlim = c(min(df[,1:3]), max(df[,1:3])), 
#'      type = "n", ylab = "Water depth [m]", main = "Hydrography and free vertical range (shaded area)",
#'      xlab = "Temperature (T) / Oxygen (O) / Salinty (S)")
#' rect(0, max(FVR.a), 20, min(FVR.a), angle = 45, density = 6, col = "grey70")
#' for(i in 1:3){
#'   c <- c("red", "blue", "darkgreen")
#'   lines(df[,i], as.numeric(rownames(df)), lty = i, lwd = 2, col = c[i])
#'   posY <- round(nrow(df) * 0.9, 0)
#'   points(df[posY,i], posY, pch = 21, bg = "white", cex = 4)
#'   text(df[posY,i], posY, substr(colnames(df)[i], 1, 1))
#' }
#'
#' #-2-----Example for dynamic hydrography (DynHyd) including intra daily light regimes 
#' #-2.1---Generating interpolated hydrography profiles for CTDs by use of InterPro()
#' set.seed(123)
#' hyd.a <- HydroGen(data.frame(Temp = c(12,7,5,7,8,9), Ox = c(8,7,5,3,2,0)), by = 20)
#' hyd.b <- HydroGen(data.frame(Temp = c(10,7,4.8,7,8), Ox = c(8,7,6,5,0)), by = 25)
#' hyd.c <- HydroGen(data.frame(Temp = c(11,4,8), Ox = c(8,4,0)), by = 50)
#' hyd.d <- HydroGen(data.frame(Temp = c(11,7,4,6,8), Ox = c(8,7,6,4,0)), by = 25)
#' hyds <- list(hyd.a, hyd.b, hyd.c, hyd.d)
#' tp <- matrix(c(1,2,2,3,3,4), ncol = 3)
#' dyn.hyd <- list()
#' for(i in 1:3){dyn.hyd <- c(dyn.hyd, InterPro(data = list(hyds[[tp[1,i]]], hyds[[tp[2,i]]]), 
#'                                              steps = 59, opList = TRUE))}
#'
#' #-2.2---Generating light regimes (19:00 - 22:00) by use of date-at-location specific regression 
#' light <- LightBornholm(depth   = c(0:99), ATTk = -0.16, lx = TRUE,
#'                        daytime = seq(19,22,0.01666667), min = 2)
#'
#' for(i in 1:180){dyn.hyd[[i]] <- cbind(dyn.hyd[[i]], light[,i], light[,i])
#' colnames(dyn.hyd[[i]]) <- c(colnames(hyd.a), "UpperLight", "LowerLight")}
#'
#' hyd <- array(unlist(dyn.hyd), dim = c(dim(dyn.hyd[[1]]), length(dyn.hyd)), 
#'              dimnames = list(NULL,colnames(dyn.hyd[[1]]), NULL))
#'
#' #-2.3---Generating factor list for any kind of organism (e.g. European sprat)
#' fac      <- c("UL",  "LL", "MOV", "RI")
#' temp     <- c(  4,     5,      0,  0.5)
#' ox       <- c(   1,   0.5,     1,  0.5)
#' lightUp  <- c( 500,    10,    -1,  0.2)
#' lightLo  <- c( 0.1, 0.005,     1,  0.2)
#' fl       <- data.frame(rbind(temp, ox, lightUp, lightLo))
#' colnames(fl) <- fac
#' rownames(fl) <-  colnames(hyd)
#'
#' #-2.3---Executing EvalFunc() for given hydro (including light regimes) and factor list
#' df     <- EvalFunc(hydro = hyd, faclist = fl)
#' head(df[,,1],10)
#'
#' #-2.4---Plotting
#' lists  <- EvalFunc(hydro = hyd, faclist = fl, dimList = 2)
#' PVD    <- lists$PredVerDi
#' image  <- as.raster(PVD)
#' time   <- format(seq(from = as.POSIXct("2001-06-04 19:00"), 
#'                      to = as.POSIXct("2001-06-04 22:00"), 
#'                      by = "min"), "%H:%M")[c(1,31,61, 91,121,151,181)]
#'
#' AS <- function(m) t(m)[,nrow(m):1]
#' FS <- function(x) (x-min(x))/(max(x)-min(x)) 
#'
#' layout(matrix(c(1,2,3), 3, 1), 
#'        widths=c(1,1,1), heights=c(1,1,2))
#'
#' #-2.4.1-Temperature profile
#' par(mar = c(0.5,4.1,2.1,2.1))
#' colfunc <- colorRampPalette(c("blue1", "lightblue", "limegreen", "yellowgreen", 
#'                               "yellow", "orange", "red"))
#' image(AS(lists$Temp), useRaster=TRUE, col = colfunc(500), axes = FALSE, ylab = "Water depth [m]")
#' legend(-0.04, 1.05, "TEMP",  bty="n", cex = 1.8, xjust = 0, yjust = 1)
#' abline(v = seq(0,1,length.out = 7), h = seq(0,1,0.2), col = "grey70", lty = 3)
#' axis(2, at = seq(0,1,0.2), labels = seq(100,0,-20), las = 2); box()
#'
#' #-2.4.2-Oxygen profile
#' par(mar = c(2,4.1,0.5,2.1))
#' colfunc <- colorRampPalette(c("slategray4", "slategray", "slategray3", "slategray2", 
#'                               "slategray1", "white", "orangered1"))
#' image(AS(lists$Ox), useRaster=TRUE, col = rev(colfunc(500)), axes = FALSE, ylab = "Water depth [m]")
#' legend(-0.04, 1.05, "OX",  bty="n", cex = 1.8, xjust = 0, yjust = 1)
#' abline(v = seq(0,1,length.out = 7), h = seq(0,1,0.2), col = "grey70", lty = 3)
#' axis(2, at = seq(0,1,0.2), labels = seq(100,0,-20), las = 2); box()
#'
#' #-2.4.2-Predicted vertical distribution of sprat - high probabilty (bright araes) vs. low prob (dark areas)
#' par(mar = c(5.1,4.1,0.5,2.1))
#' colfunc <- colorRampPalette(c("black", "white"))
#' image(AS(lists$PredVerDi), useRaster=TRUE, col = colfunc(500), axes = FALSE, ylab = "Water depth [m]", 
#'       xlab = "Day time [UTC+2]")
#' legend(-0.04, 1.03, "MODEL",  bty="n", cex = 1.8, xjust = 0, yjust = 1)
#' abline(v = seq(0,1,length.out = 7), h = seq(0,1,0.2), col = "grey70", lty = 3)
#' axis(2, at = seq(0,1,0.2), labels = seq(100,0,-20), las = 2)
#' axis(1, at = seq(0,1,length.out = 7), labels = time); box()
#'
#' #-2.4.2-Adding movement pattern for 60 individuals of European sprat
#' set.seed(12345)
#' nIndivi <- 60
#' testInd <- EvalFunc(hydro = hyd, faclist = fl, nInd = nIndivi, mSpeed = 3, rWalk = c("ND", 0, 0.1, 0.03))
#' colfunc <- colorRampPalette(c("blue", "limegreen", "yellow", "orange", "red"))
#' for(i in 1:nIndivi){
#'   lines(FS(c(1:length(testInd[i,]))), 1-testInd[i,]/nrow(lists$PredVerDi), 
#'         col = colfunc(nIndivi)[i], cex = 1.2)
#' }
#' 
#' par(mfrow = c(1,1))
#' }
#' @export

EvalFunc  <- function(hydro, faclist, nInd, mSpeed = 1, tStep = 1, rWalk = c("ND", 0, 0.1, 0.1), dimList){
  
  #Check class of object hydro and transform to array
  if(is.data.frame(hydro) || is.matrix(hydro)){
    hydro <- array(unlist(list(hydro)), dim = c(dim(hydro), 1), dimnames = list(NULL, colnames(hydro), NULL))
  } else if(is.list(hydro)){hydro <- array(unlist(hydro), dim = c(dim(hydro[[1]]), length(hydro)), 
                                           dimnames = list(NULL, colnames(hydro[[1]]), NULL))}
  
  #Execute TransFunc() and sum that up
  PreMov <- array(numeric(), dim = c(nrow(hydro), nrow(faclist)+2, dim(hydro)[3]))
  for(i in 1:ncol(hydro)){
    PreMov[,i,] <- TransFunc(par = hydro[,i,], 
                             UL = faclist$UL[i], 
                             LL = faclist$LL[i],
                             MOV = faclist$MOV[i],
                             RI = faclist$RI[i])
  }
  
  rSums <- matrix(t(apply(PreMov, 1, colSums, na.rm = TRUE)))
  
  FeatScale <- function(x){(x-min(x))/(max(x)-min(x))}
  PVD       <- abs(FeatScale(matrix(t(apply(abs(PreMov), 1, colSums, na.rm = TRUE))))-1)
  
  PreMov[,ncol(PreMov)-1,] <- rSums
  PreMov[,ncol(PreMov),]   <- PVD
  
  #Returns array or list including hydro, parameter dependent individual movement and PVD
  if(missing(nInd)){
    hydArray <- array(numeric(), dim = c(nrow(hydro), ncol(PreMov)+ncol(hydro), dim(hydro)[3]))
    for(i in 1:dim(hydro)[3]){
      hydArray[,,i] <- cbind(hydro[,,i], PreMov[,,i])
    }
    n <- matrix(c(rep("PredMovBy", nrow(faclist)), rownames(faclist)), nrow(faclist))
    dimnames(hydArray) <- list(NULL, c(colnames(hydro), apply(n, 1, paste, collapse=""),"PredMovByAllPar", "PredVerDi"), NULL)
    
    if(!missing(dimList)){
      split.along.dim <- function(a, n){
        setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                        array, dim = dim(a)[-n], dimnames(a)[-n]),
                 dimnames(a)[[n]])
      }
      hydList <- split.along.dim(hydArray, n = dimList)
      return(hydList)
    }
    return(hydArray)
  }
  
  #Calculates movements of individual organisms - step by step
  else{
    
    if(dim(hydro)[3] > 1){tStep <- dim(hydro)[3]}
    if(length(nInd) == 1){
      indlist <- matrix(runif(nInd, 1, nrow(hydro)), nrow = nInd, ncol = tStep + 1)
    }else{indlist <- matrix(nInd, nrow = length(nInd), ncol = tStep + 1)}
    
    for(i in 1:tStep){
      vom <- PreMov[indlist[,i], ncol(PreMov)-1, ifelse(dim(hydro)[3] == 1, 1, i)]
      if(is.null(rWalk)){
        mov <- vom
      }else{mov <- sapply(vom, RWalk, op = as.character(rWalk[1]), 
                          mv = as.numeric(rWalk[2]), 
                          sd = as.numeric(rWalk[3]),
                          hs = as.numeric(rWalk[4]))
      }
      indlist[,i+1] <- round(indlist[,i] + mov * mSpeed * -1, 0)
      indlist[,i+1][indlist[,i+1] < 0] <- 0
      indlist[,i+1][indlist[,i+1] > nrow(hydro)] <- nrow(hydro)
    }
    return(indlist)
  }
}

#----------------------------------------------------------------------------------
#--Calculating the general areal overlap among two species-------------------------
#----------------------------------------------------------------------------------

#' @title General areal overlap
#' @description Function calculates interspecific overlap among two populations.
#' @param data matrix or data frame including two colums - one for each population's distribution.
#' @details The general areal overlap is defined as the area occupied by both populations devided by 
#' the area occupied by either one or the other population.
#' @examples 
#' #Example shows how two populations overlap in one of five water layers
#' 
#' df <- cbind(c(1,1,1,0,0),
#'             c(0,0,1,1,1))
#' 
#' GeneralArealOverlap(df)
#' @export

GeneralArealOverlap <- function(data){
  
  a <- table(ifelse(data[,1] >= 0.5, 1, 0))["1"]
  b <- table(ifelse(data[,2] >= 0.5, 1, 0))["1"]
  c <- table(ifelse(data[,1] >= 0.5 & data[,2] >= 0.5, 1, 0))["1"]
  d <- c/(a+b-c)
  d[is.na(d)] <- 0
  return(d)
  
}

#----------------------------------------------------------------------------------
#-------Hydrography generator------------------------------------------------------
#----------------------------------------------------------------------------------

#' @title Hydrography generator
#' @description Function generates customized hydrography for test purposes.
#' @param data data frame including environmental parameters.
#' @param by numeric, interval between each data point, Default: 10.
#' @details Generates a customized hydrography. More details in 'Examples'.
#' @examples 
#' #Example generates hydro including temperature, salinity and oxygen
#' 
#' data <- data.frame(Temperature = c(17, 15, 10,   8,  4, 4.5,   6,  8),
#'                    Oxygen =      c(6,   5,  4,   3,  2,   1, 1.5,  0),
#'                    Salinity =    c(7, 7.5,  8, 8.5,  9, 9.5,  14, 20))
#'                    
#' HydroGen(data = data, by = 10)
#' @export

HydroGen <- function(data, by = 10){
  
  hydro <- matrix(ncol = ncol(data), nrow = nrow(data)*by-by)
  for(j in 1:ncol(data)){
    rows <- nrow(data)-1
    asd <- NULL
    for(i in c(1:rows)){
      asd <- c(asd, sort(runif(by, min = min(data[i:(i+1), j]), max = max(data[i:(i+1), j])), 
                         decreasing = data[i, j] > data[i+1, j]))
    }
    hydro[,j] <- asd  
  }
  hydro <- as.data.frame(hydro)
  colnames(hydro) <- colnames(data)
  return(hydro)
}

#----------------------------------------------------------------------------------
#-------Interpolating between two CTD profiles-------------------------------------
#----------------------------------------------------------------------------------

#' @title Hydrography interpolator
#' @description Function computes linear interpolation between two hydrographies.  
#' @param data list of two hydrography profiles.
#' @param steps numeric, number of profiles to be interpolated.
#' @param opList logical, if TRUE, returns list. If FALSE, returns array.
#' @details To be written. Further details in 'Examples'.
#' @examples 
#' \dontrun{
#' 
#' #Example for two customized CTDs along south-to-north transect (20 nautical miles)
#' #Computation of linearly interpolated hydrographies for every nautical mile
#'
#' #Generating two customized hydrographies
#' set.seed(12345)
#' hydNorth <- HydroGen(data.frame(Temperature = c(17, 10, 4, 7, 9),
#'                                Oxygen = c(6, 5, 3, 2, 0),
#'                                Salinity = c(7, 8, 9, 18, 20)), by = 20)
#'
#' hydSouth <- HydroGen(data.frame(Temperature = c(20, 17, 3, 7, 8),
#'                                 Oxygen = c(7, 6, 4, 1, 0),
#'                                 Salinity = c(6, 8, 9, 16, 20)), by = 20)
#'
#' #Generating list including hydNorth and hydSouth
#' dl  <- list(hydNorth, hydSouth)
#'
#' #Calculating nineteen interpolated profiles
#' hydro <- InterPro(data = dl, steps = 19, opList = TRUE)
#'
#' #Generating plot including temperature and salinity profiles from southern towards northern station 
#' par(mfrow = c(1,2))
#' for(j in c(1,3)){
#'   plot(x = 0, xlim = c(-20, 20), ylim = c(85,0), ylab = ifelse(j == 1, "Water depth", ""), xlab = "", xaxt='n', type = "n")
#'   points(c(-18, 18), c(84, 84), pch = 21, cex = 4)
#'   text(c(-18, 18), c(84, 84), c("S", "N"), cex = 1.5)
#'   title(main = c("Interpolated temperature profiles\nfrom one station to another", "", 
#'                  "Interpolated salinity profiles\nfrom one station to another")[j])
#'   
#'   for(i in 20:1){
#'     lines(hydro[[i]][,j]-i, c(1:80), col = c("red3", "", "green4")[j])
#'   }
#'  }
#' }
#' @export

InterPro <- function(data, steps, opList = FALSE){
  
  if(nrow(data[[1]]) != nrow(data[[2]])){
    nr <- abs(nrow(data[[1]]) - nrow(data[[2]]))
    ifelse(nrow(data[[1]]) < nrow(data[[2]]),
           data[[1]] <- rbind(data[[1]], matrix(rep(c(tail(data[[1]], 1)), nr), 
                                                ncol = ncol(data[[1]]), byrow = TRUE)),
           data[[2]] <- rbind(data[[2]], matrix(rep(c(tail(data[[2]], 1)), nr), 
                                                ncol = ncol(data[[2]]), byrow = TRUE)))
  }
  
  dataOne    <- data[[1]]
  dataTwo    <- data[[2]]
  sum        <- (dataOne - dataTwo) / steps * -1
  hydroList  <- list(dataOne)
  
  for(i in 1:steps){
    hydroList[[i+1]]  <- hydroList[[i]] + sum
  }
  
  ifelse(opList == TRUE,
         return(hydroList),
         hydroArray <- array(unlist(hydroList), dim = c(dim(hydroList[1]), length(hydroList))))
  
  return(hydroArray)
}

#----------------------------------------------------------------------------------
#-------Light issues---------------------------------------------------------------
#----------------------------------------------------------------------------------

#' @title Light regime generator
#' @description Function generates semi-customized light profile typically found in the Bornholm Basin (Central Baltic Sea) in early May.  
#' @param depth numeric vector, depth profile which is used by light model, Default: c(1:100).
#' @param daytime numeric vector, time sequence in decimals, e.g. every single minute over 24 hours (default).
#' @param ATTk numeric, attenuation coefficient.
#' @param min numeric, minimum illumination at shallowest depth.
#' @param lx logical, if TRUE, output is in lux [lx], if FALSE, computations are in W/m2 (default).
#' @details Function makes use of light model developed by (...). Further details in 'Examples'.
#' @examples 
#' \dontrun{
#' 
#' #The minimum illumination at >0m of depth is assumed to be 2 lx
#' minVal <- 2
#' 
#' #Attenuation coefficiant is set to -0.16, which can be considered as relatively 'clear'
#' atteCo <- -0.16
#'
#' #Computing intra daily under water light regime [lx] down to 100m of depth in Bornholm Basin
#' testmatrix <- LightBornholm(depth = c(1:100), ATTk = atteCo, lx = TRUE, min = minVal)
#'
#' #Assigning time of day to each column
#' ts <- format(seq(from = as.POSIXct("2017-08-19 0:00"), length.out = 1440, by = "min"), "%H:%M")
#' colnames(testmatrix) <- ts
#' 
#' #Printing underwater light regime between 20 and 40 meters of depth at 00:00, 06:00, 12:00 and 18:00
#' testmatrix[20:40, c("00:00", "06:00", "12:00", "18:00")]
#' }
#' @export

LightBornholm <- function(depth, daytime = seq(0,24,0.01666667), ATTk, min, lx){
  
  a <- exp(ATTk * depth)    
  if(lx){b <- 393.199 * sin(2*pi*daytime/30.3702+5.15129) * (4.6 / 0.01953)
  }else{b <- 393.199 * sin(2*pi*daytime/30.3702+5.15129)}
  c <- ifelse(b < min, min, b)
  return(a %o% c)
}

#------------------------------------------------------------------------------------
#---------Calculating gradients------------------------------------------------------
#------------------------------------------------------------------------------------

#' @title Calculating gradients
#' @description Function calculates in which direction and at what rate a certain parameter changes over depth.
#' @param depth numeric vector of depth indication.
#' @param par   numeric vector of single environmental parameter (e.g. water temperature).
#' @param range numeric integer indicating the range for which the delta should be calculated, Default: 5.
#' @examples 
#' #Generating customized temperature and salinity profile
#' hydro <- HydroGen(data.frame(temp = c(17,10,6,8,9),
#'                              salt = c(7,8,9,18,20)),
#'                              by = 25)
#'                                 
#' #Calculating salinity and temperature gradient on a 3-meter-range
#' hydro$tempGrad <- parGrad(depth = seq_along(hydro$temp), par = hydro$temp, range = 3)
#' hydro$saltGrad <- parGrad(depth = seq_along(hydro$salt), par = hydro$salt, range = 3)
#' 
#' head(hydro, 15)
#' @rdname parGrad
#' @export

parGrad <- function(depth, par, range = 5) {
  
  suppressPackageStartupMessages(library(data.table))
  
  if(is.unsorted(depth)) {
    stop("Data must be sorted by depth...")}                
  else{
    ifelse(is.na((shift(par,type = "lead") -          
                    shift(par)) *                         
                   (range/abs(shift(depth, type = "lead") - 
                                shift(depth)))), 0, 
           (shift(par, type = "lead") -                   
              shift(par)) *                           
             (range/abs(shift(depth, type = "lead") - 
                          shift(depth))) )
  }
}

#----------------------------------------------------------------------------------
#--Calculating the specific areal overlap among two species------------------------
#----------------------------------------------------------------------------------

#' @title Specific areal overlap
#' @description Function calculates fraction of areas occupied by one population on areas occupied by 
#' another population.
#' @param data matrix or data frame including two colums - one for each population's distribution.
#' @param col  numeric integer indicating the column in which the species if interest in given, Default = 1.
#' @details The specific areal overlap is defined as the area occupied by both populations devided by 
#' the area occupied by one of these populations.
#' @examples 
#' #Example shows water column occupied by clupeid fish and cod
#' 
#' df <- data.frame(Clupeids = c(0,1,1,1,1,0),
#'                       Cod = c(0,0,0,1,1,1))
#' 
#' #Calculating specific areal overlap for clupeids and cod, respectively
#' 
#' clu <- SpecificArealOverlap(df, col = 1) * 100
#' cod <- SpecificArealOverlap(df, col = 2) * 100
#' 
#' paste(clu, "% of the water layers occupied by clupeid fish are also occupied by cod.", sep = "")
#' paste(cod, "% of the water layers occupied by cod are also occupied by clupeid fish.", sep = "")
#' @export

SpecificArealOverlap <- function(data, col = 1){
  
  colTwo <- ifelse(col == 1, 2, 1)
  a <- table(ifelse(data[,col] >= 0.5, 1, 0))["1"]
  b <- table(ifelse(data[,col] >= 0.5 & data[,colTwo] >= 0.5, 1, 0))["1"]
  c <- b/a
  c[is.na(c)] <- 0
  return(c)
  
}

#----------------------------------------------------------------------------------
#---------Transition function------------------------------------------------------
#----------------------------------------------------------------------------------

#' @title Transition function
#' @description Function computes organism's vertical movement based on one single environmental factor. 
#' @param par numeric vector of single environmental parameter.
#' @param UL  numeric, species specific parameter threshold located in upper water column.
#' @param LL  numeric, species specific parameter threshold located in lower water column.
#' @param MOV numeric, direction of movement when threshold is reached, Default: -1.
#' @param RI  numeric, parameter weight, Default: 1.
#' @details Code generates transition function and applies that to the vertical profile of a 
#' single environmental factor (e.g. water temperature). Further details in example below...
#' @examples 
#' \dontrun{
#' 
#' #Example for induced species' movement by ambient salinity and oxygen content (par)
#' #Generating customized salinity and oxygen profile by executing HydroGen()
#' set.seed(1231)
#' hydFic <- HydroGen(data = data.frame(Salinity = c(7,8,9,18,20),
#'                                      Oxygen = c(6,5,4,3,0)), by = 25)
#'  
#' #Generating a bunch of randomly distributed individuals (n = 20)
#' ind <- data.frame(runif(20, 1, nrow(hydFic)))
#'  
#' #Executing transition functions by applying TransFunc()
#' #Salinity: Maximum strength of downward movement (MOV = -1) when par <= upper limit (UL)
#' #Oxygen: Maximum strength of upward movement (MOV = 1) when par <= lower limit (LL)
#' 
#' movSal <- TransFunc(par = hydFic$Salinity, UL = 8, LL = 11)
#' movOx  <- TransFunc(par = hydFic$Oxygen, UL = 2, LL = 1, MOV = 1)
#'  
#' #Assigning strength and direction of movement to each individual
#' ind$movSal <- movSal[ind[,1]]
#' ind$movOx  <- movOx[ind[,1]]
#' colnames(ind) <- c("CurrentDepthOfIndividual", "PredictedMovementBy:Sal", "PredictedMovementBy:Ox")
#' head(ind, 10)
#'  
#' #Setting graphical paramaters
#' sig        <- c(25, 21, 24)
#' pa         <- cbind(hydFic$Salinity, hydFic$Oxygen)
#' col        <- c("limegreen", "blue")
#' ind$simSal <- sig[sign(movSal[ind[,1]])+2]
#' ind$simOx  <- sig[sign(movOx[ind[,1]])+2]
#' lsal       <- c(tail(which(hydFic$Salinity <= 8), 1), tail(which(hydFic$Salinity <= 11), 1))
#' lox        <- c(head(which(hydFic$Oxygen <= 2), 1), head(which(hydFic$Oxygen <= 1), 1))
#' limline    <- cbind(lsal, lox)
#' tit        <- c("Salinity profile including", "Oxygen profile including", 
#'                 "predicted individual movement\nand its strength")
#' xl         <- c("Practical salinity [PSU]", "Oxygen content [ml/l]")
#' yl         <- c("Water depth [m]", "")
#'  
#' #Plotting
#' par(mfrow = c(1,2),
#'     mar = c(5.1, 4.1, 4.1, 2.1))
#' for(i in 1:2){
#'   plot(pa[,i], c(1:length(pa[,i])), ylim = c(length(pa[,i]),1), xlab = xl[i], ylab = yl[i],
#'        type = "l", col = col[i], lwd = 2, main = paste(tit[i], tit[3], sep = "\n"))
#'   abline(h = limline[,i], col = "gray60", lty = 2)
#'   text(c(19, 5.5)[i], c(limline[1,i]+2, limline[2,i]+2), c("UL", "LL"))
#'   points(seq(min(pa[,i]), max(pa[,i]), length.out = 20), ind[,1], 
#'          pch = ind[,i+3], cex = (ind[,i+1] * -1.3) ^ 4 + 2, bg = "grey")
#'   }
#' }
#' @export

TransFunc <- function(par, UL, LL, MOV = -1, RI = 1){
  
  #AdditionalFunction
  parGradExtention <- function(par){
    depth <- seq_along(par)
    parGrad(depth, par = par, range = 5)
  }
  
  #PreCalculations
  shape <- sign(UL - LL)                      
  a <- 1 / (UL - LL)                                             
  b <- -1 * a * LL
  
  x <- pmax(shape * par, shape * LL) * shape
  x <- pmin(shape *    x, shape * UL) * shape
  factor <- a * x + b
  
  if(MOV == 0 & is.matrix(factor)){
    pargrad <- apply(par, 2, parGradExtention)
  }else{
    pargrad <- parGrad(seq_along(par), par = par, range = 5)
  }
  
  #FinalCalculations
  if(MOV == 0){
    return(((sign(pargrad) * -1) * factor) * RI)
  }
  
  if(MOV != 0){
    if(MOV == -1){(MOV * factor) * RI} 
    else{
      if(MOV == 1 & UL > LL){(MOV - factor) * RI} 
      else{(1 - factor) * RI}
    }
  }
}