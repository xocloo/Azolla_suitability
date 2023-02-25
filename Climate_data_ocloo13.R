setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature") #this is the folder with all the GeoTiffs

library(bbmle)
library(RColorBrewer)


azolla = read.csv("Azolla_temp.csv")
fili = subset(azolla, Species == "A. filiculoides")

Room = function(t, Pm,  Topt, a, b){
  if(t == Topt){return(Pm)}else{
    if(t < Topt){return(Pm*exp(-a*(t - Topt)^2))}else{
      return(Pm*exp(-b*(t - Topt)^2))
    }
  }
}

Room_LL = function(Pm, Topt, a, b, sd, data){
  predictions = mapply(FUN = Room, t=data$Temp, MoreArgs = list(Pm=Pm, Topt=Topt, a=a, b=b))
  sum(dnorm(data$rgr, mean=predictions, sd=sd, log=T))
}

Room_LL(Pm=0.30, Topt=30, a=0.1, b=0.1, sd=0.1, data=subset(fili, Study==54))

Room_ILL = function(Pm_mean, Topt, a, b, sd_Pm, sd_obs, n=50, data){
  sd_Pm = abs(sd_Pm)
  sd_obs = abs(sd_obs)
  qs = seq(from= 0.001, to = 0.999, length.out = n) #select quantiles
  Pms = qnorm(qs, mean=Pm_mean, sd=sd_Pm) # obtain random effect values
  if(min(c(Pms, Topt, a, b, sd_Pm, sd_obs)) <= 0){return(-1e6)}
  ds = dnorm(Pms, mean=Pm_mean, sd=sd_Pm) # get the density of these values
  probs = ds/sum(ds)
  delta_x = qs[2] - qs[1]
  LL = numeric()
  ### This band is just for plotting
  # plot(rgr ~ Temp, data=data, pch=21, bg="black", ylim=c(0, 0.75))
  temps =  seq(from=5, to=50, by=0.1)
  # ys = mapply(FUN=Room, t=temps, MoreArgs = list(a=a, Pm=Pm_mean, Topt=Topt, b=b))
  # lines(ys ~ temps, col="black")
  ### end band
  for(i in 1:length(Pms)){
    LL[i] = Room_LL(Pms[i], Topt, a, b, sd_obs, data) # total LL for a given Pm
    ### This band is just for plotting
    ys = mapply(FUN=Room, t=temps, MoreArgs = list(a=a, Pm=Pms[i], Topt=Topt, b=b))
    lines(ys ~ temps, col="red")
    ### end band
  }
  LL_max = max(LL)
  relative_probs = sum(exp(LL - LL_max))
  (LL_max + log(relative_probs))
}

Room_ILL(Pm_mean = 0.4, Topt=24, a = 0.007, b = 0.007, sd_Pm = 0.1, n=10, sd_obs=0.07, data=subset(fili, Study==54))

study_list = unique(fili$Study)

Azolla_total_ILL = function(Pm_mean, Topt, a, b, sd_Pm, sd_obs, n=50, data=fili, NLL=TRUE){
  study_list = unique(data$Study)
  total_ILL = 0
  plot(rgr ~ Temp, data=data, pch=21, bg="black", ylim=c(0, 0.75))
  for(i in 1:length(study_list)){
    ILL = Room_ILL(Pm_mean, Topt, a, b, sd_Pm, sd_obs, n, data=subset(data, Study == study_list[i]))
    total_ILL = total_ILL + ILL
  }
  print(total_ILL)
  ifelse(NLL==T, -total_ILL, total_ILL)
}

Azolla_total_ILL(Pm_mean = 0.25, Topt=24, a = 0.007, b = 0.007, sd_Pm = 0.05, sd_obs=0.07)

# hard coded in data for filiculoides
m2 = mle2(minuslogl = Azolla_total_ILL, start=list(Pm_mean = 0.19, Topt=23.5, a = 0.005, b = 0.005, sd_Pm = 0.05, sd_obs=0.04),
          control=list(parscale=c(Pm_mean = 0.19, Topt=23.5, a = 0.005, b = 0.005, sd_Pm = 0.05, sd_obs=0.04), trace=2))
m2


Azolla_total_ILL(Pm_mean = coef(m2)["Pm_mean"], Topt=coef(m2)["Topt"], a = coef(m2)["a"], b = coef(m2)["b"], 
                 sd_Pm = coef(m2)["sd_Pm"], sd_obs=coef(m2)["sd_obs"])


pinnata = subset(azolla, Species == "A. pinnata")
  
Azolla_total_ILL = function(Pm_mean, Topt, a, b, sd_Pm, sd_obs, n=50, data=pinnata, NLL=TRUE){
  study_list = unique(data$Study)
  total_ILL = 0
  plot(rgr ~ Temp, data=data, pch=21, bg="black", ylim=c(0, 0.75))
  for(i in 1:length(study_list)){
    ILL = Room_ILL(Pm_mean, Topt, a, b, sd_Pm, sd_obs, n, data=subset(data, Study == study_list[i]))
    total_ILL = total_ILL + ILL
  }
  ifelse(NLL==T, -total_ILL, total_ILL)
}

Azolla_total_ILL(Pm_mean = 0.25, Topt=24, a = 0.007, b = 0.007, sd_Pm = 0.05, sd_obs=0.07)

# hard coded in data for pinnata
m3 = mle2(minuslogl = Azolla_total_ILL, start=list(Pm_mean = 0.25, Topt=27, a = 0.004, b = 0.004, sd_Pm = 0.01, sd_obs=0.01),
          control=list(parscale=c(Pm_mean = 0.25, Topt=27, a = 0.004, b = 0.004, sd_Pm = 0.01, sd_obs=0.001), maxit=1e5, trace=T))
m3

Azolla_total_ILL(Pm_mean = coef(m3)["Pm_mean"], Topt=coef(m3)["Topt"], a = coef(m3)["a"], b = coef(m3)["b"], 
                 sd_Pm = coef(m3)["sd_Pm"], sd_obs=coef(m3)["sd_obs"])

azolla_fit = function(Tmin=5, Tmax=50, df=fili, fit=m2){
  temps =  seq(from=Tmin, to=50, by=0.1)
  ys = mapply(FUN=Room, t=temps, MoreArgs = list(a=coef(fit)["a"], Pm=coef(fit)["Pm_mean"], Topt=coef(fit)["Topt"], b=coef(fit)["b"]))
  plot(rgr ~ Temp, data=df, pch=21, bg="black", ylim=c(0, 0.4), xlim=c(Tmin, Tmax))
  lines(ys ~ temps, col="red", lwd = 5)
  
}
azolla_fit(Tmin=5, Tmax=50, df=fili, fit=m2) # For filiculoides
azolla_fit(Tmin=5, Tmax=50, df=pinnata, fit=m3) #For pinnata


#You will get back the dimensions of those data sets. Each line of 
#code should give you back 2 numbers. The first number should be the 
#larger number in both cases, and it will be the number of observations for both data sets
dim(fili) #149
dim(pinnata) #40



##############CURRENT average plot monthly for Senegal at 30s resolution#######
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/Current tavg_30s")
#install.packages("raster")
library(raster)
library(sp)
#install.packages("rasterVis")
library("rasterVis")
#install.packages("rgdal")
library(rgdal)
library(RColorBrewer)

meanJan <- raster("wc2.1_30s_tavg_01.tif") #mean temp for January so on so forth..
meanFeb <- raster("wc2.1_30s_tavg_02.tif")
meanMar <- raster("wc2.1_30s_tavg_03.tif")
meanApr <- raster("wc2.1_30s_tavg_04.tif")
meanMaj <- raster("wc2.1_30s_tavg_05.tif")
meanJun <- raster("wc2.1_30s_tavg_06.tif")
meanJul <- raster("wc2.1_30s_tavg_07.tif")
meanAvg <- raster("wc2.1_30s_tavg_08.tif")
meanSep <- raster("wc2.1_30s_tavg_09.tif")
meanOkt <- raster("wc2.1_30s_tavg_10.tif")
meanNov <- raster("wc2.1_30s_tavg_11.tif")
meanDec <- raster("wc2.1_30s_tavg_12.tif")
#plot(meanJan) #will simply get the basic R plot of a GeoTiff raster
#plot(meanAvg)

##zooming in 
mtStack <- stack(meanJan, meanFeb, meanMar, meanApr, meanMaj, meanJun, meanJul, meanAvg, meanSep, meanOkt, meanNov, meanDec)
plot(meanJan) #plot the raster file for any month
extentA <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack <- (crop(mtStack, extentA)) #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack) <- month.abb #adds month abbreviations
mapTheme2 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack, main = "Current tavg Temperature", par.settings = mapTheme2 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 
library(RColorBrewer)
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))


levelplot(r2, main = "Monthly tavg Relative Growth Rate of Azolla filiculoides", par.settings = mapTheme) 


## Pinnata monthly
Pm <- 0.252 #inputing parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2, main = "Current tavg Relative Growth Rate of Azolla pinnata", par.settings = mapTheme) 


#####CURRENT MAX SENEGAL 30 s#######
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/current tmax_30s/wc2.1_30s_tmax")
meanJan <- raster("wc2.1_30s_tmax_01.tif") #mean temp for January so on so forth..
meanFeb <- raster("wc2.1_30s_tmax_02.tif")
meanMar <- raster("wc2.1_30s_tmax_03.tif")
meanApr <- raster("wc2.1_30s_tmax_04.tif")
meanMaj <- raster("wc2.1_30s_tmax_05.tif")
meanJun <- raster("wc2.1_30s_tmax_06.tif")
meanJul <- raster("wc2.1_30s_tmax_07.tif")
meanAvg <- raster("wc2.1_30s_tmax_08.tif")
meanSep <- raster("wc2.1_30s_tmax_09.tif")
meanOkt <- raster("wc2.1_30s_tmax_10.tif")
meanNov <- raster("wc2.1_30s_tmax_11.tif")
meanDec <- raster("wc2.1_30s_tmax_12.tif")
#plot(meanJan) #will simply get the basic R plot of a GeoTiff raster
#plot(meanAvg)

##zooming in 
mtStack <- stack(meanJan, meanFeb, meanMar, meanApr, meanMaj, meanJun, meanJul, meanAvg, meanSep, meanOkt, meanNov, meanDec)
plot(meanJan) #plot the raster file for any month
extentA <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack <- (crop(mtStack, extentA)) #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack) <- month.abb #adds month abbreviations
mapTheme2 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack, main = "Current tx Temperature", par.settings = mapTheme2 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 
library(RColorBrewer)
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2, main = "Current tx Relative Growth Rate of Azolla filiculoides", par.settings = mapTheme) 


## Pinnata monthly
Pm <- 0.252 #inputing parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2, main = "Current tx RGR of A. pinnata", par.settings = mapTheme) 




##CURRENT MIN SENEGAL 30 S####
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/current tmin_30s/wc2.1_30s_tmin")
meanJan <- raster("wc2.1_30s_tmin_01.tif") #mean temp for January so on so forth..
meanFeb <- raster("wc2.1_30s_tmin_02.tif")
meanMar <- raster("wc2.1_30s_tmin_03.tif")
meanApr <- raster("wc2.1_30s_tmin_04.tif")
meanMaj <- raster("wc2.1_30s_tmin_05.tif")
meanJun <- raster("wc2.1_30s_tmin_06.tif")
meanJul <- raster("wc2.1_30s_tmin_07.tif")
meanAvg <- raster("wc2.1_30s_tmin_08.tif")
meanSep <- raster("wc2.1_30s_tmin_09.tif")
meanOkt <- raster("wc2.1_30s_tmin_10.tif")
meanNov <- raster("wc2.1_30s_tmin_11.tif")
meanDec <- raster("wc2.1_30s_tmin_12.tif")
#plot(meanJan) #will simply get the basic R plot of a GeoTiff raster
#plot(meanAvg)

##zooming in 
mtStack <- stack(meanJan, meanFeb, meanMar, meanApr, meanMaj, meanJun, meanJul, meanAvg, meanSep, meanOkt, meanNov, meanDec)
plot(meanJan) #plot the raster file for any month
extentA <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack <- (crop(mtStack, extentA)) #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack) <- month.abb #adds month abbreviations
mapTheme2 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack, main = "Current tn Temperature", par.settings = mapTheme2 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 
library(RColorBrewer)
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2, main = "Current tn RGR of A. filiculoides", par.settings = mapTheme) 


## Pinnata monthly
Pm <- 0.252 #inputing parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2, main = "Current tn RGR of Azolla pinnata", par.settings = mapTheme) 















###################### CURRENT AFRICA Monthly(TN/TX/TAVG) and ANNUAL Temperature ##############
library('maptools')
##pinnata
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/current tavg_5min/wc2.1_5m_tavg")
dir()
meanJan <- raster("wc2.1_5m_tavg_01.tif") #mean temp for January so on so forth..
meanFeb <- raster("wc2.1_5m_tavg_02.tif")
meanMar <- raster("wc2.1_5m_tavg_03.tif")
meanApr <- raster("wc2.1_5m_tavg_04.tif")
meanMaj <- raster("wc2.1_5m_tavg_05.tif")
meanJun <- raster("wc2.1_5m_tavg_06.tif")
meanJul <- raster("wc2.1_5m_tavg_07.tif")
meanAvg <- raster("wc2.1_5m_tavg_08.tif")
meanSep <- raster("wc2.1_5m_tavg_09.tif")
meanOkt <- raster("wc2.1_5m_tavg_10.tif")
meanNov <- raster("wc2.1_5m_tavg_11.tif")
meanDec <- raster("wc2.1_5m_tavg_12.tif")
#plot(meanJan) #will simply get the basic R plot of a GeoTiff raster
#plot(meanAvg)

##zooming in 
mtStack <- stack(meanJan, meanFeb, meanMar, meanApr, meanMaj, meanJun, meanJul, meanAvg, meanSep, meanOkt, meanNov, meanDec)
extentA <- extent(-20,60,-35,40)

##cutting, stack and drawing level plot
cutStack <- (crop(mtStack, extentA)) #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack) <- month.abb #adds month abbreviations
mapTheme2 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack, main = "Africa Current tavg Temperature", par.settings = mapTheme2 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 
library(RColorBrewer)
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2, main = "Current Monthly tavg RGR of Azolla filiculoides", par.settings = mapTheme) 


## Pinnata monthly
Pm <- 0.252 #inputing parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2, main = "Current tavg RGR of Azolla pinnata", par.settings = mapTheme) 


#####CURRENT MAX Africa 5 min resolution#######
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/Current tmax_5m/wc2.1_5m_tmax")
meanJan <- raster("wc2.1_5m_tmax_01.tif") #mean temp for January so on so forth..
meanFeb <- raster("wc2.1_5m_tmax_02.tif")
meanMar <- raster("wc2.1_5m_tmax_03.tif")
meanApr <- raster("wc2.1_5m_tmax_04.tif")
meanMaj <- raster("wc2.1_5m_tmax_05.tif")
meanJun <- raster("wc2.1_5m_tmax_06.tif")
meanJul <- raster("wc2.1_5m_tmax_07.tif")
meanAvg <- raster("wc2.1_5m_tmax_08.tif")
meanSep <- raster("wc2.1_5m_tmax_09.tif")
meanOkt <- raster("wc2.1_5m_tmax_10.tif")
meanNov <- raster("wc2.1_5m_tmax_11.tif")
meanDec <- raster("wc2.1_5m_tmax_12.tif")
#plot(meanJan) #will simply get the basic R plot of a GeoTiff raster
#plot(meanAvg)

##zooming in 
mtStack <- stack(meanJan, meanFeb, meanMar, meanApr, meanMaj, meanJun, meanJul, meanAvg, meanSep, meanOkt, meanNov, meanDec)
extentA <- extent(-20,60,-35,40)

##cutting, stack and drawing level plot
cutStack <- (crop(mtStack, extentA)) #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack) <- month.abb #adds month abbreviations
mapTheme2 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack, main = "Africa Current tx Temperature", par.settings = mapTheme2 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2, main = "Current tx RGR of Azolla filiculoides", par.settings = mapTheme) 


## Pinnata monthly
Pm <- 0.252 #inputing parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2, main = "Current tx RGR of Azolla pinnata", par.settings = mapTheme) 




##CURRENT MIN AFRICA  5MINN####
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/Current tmin_5m/wc2.1_5m_tmin")
meanJan <- raster("wc2.1_5m_tmin_01.tif") #mean temp for January so on so forth..
meanFeb <- raster("wc2.1_5m_tmin_02.tif")
meanMar <- raster("wc2.1_5m_tmin_03.tif")
meanApr <- raster("wc2.1_5m_tmin_04.tif")
meanMaj <- raster("wc2.1_5m_tmin_05.tif")
meanJun <- raster("wc2.1_5m_tmin_06.tif")
meanJul <- raster("wc2.1_5m_tmin_07.tif")
meanAvg <- raster("wc2.1_5m_tmin_08.tif")
meanSep <- raster("wc2.1_5m_tmin_09.tif")
meanOkt <- raster("wc2.1_5m_tmin_10.tif")
meanNov <- raster("wc2.1_5m_tmin_11.tif")
meanDec <- raster("wc2.1_5m_tmin_12.tif")
#plot(meanJan) #will simply get the basic R plot of a GeoTiff raster
#plot(meanAvg)

##zooming in 
mtStack <- stack(meanJan, meanFeb, meanMar, meanApr, meanMaj, meanJun, meanJul, meanAvg, meanSep, meanOkt, meanNov, meanDec)
extentA <- extent(-20,60,-35,40)

##cutting, stack and drawing level plot
cutStack <- (crop(mtStack, extentA)) #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack) <- month.abb #adds month abbreviations
mapTheme2 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack, main = "Africa Current tn Temperature", par.settings = mapTheme2 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2, main = "Current tn RGR of A. filiculoides", par.settings = mapTheme) 


## Pinnata monthly
Pm <- 0.252 #inputing parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344


r <- raster(nrow= nrow(cutStack[[1]]), ncol= ncol(cutStack[[1]]), xmn= extent(cutStack[[1]])[1], xmx=extent(cutStack[[1]])[2],
            ymn=extent(cutStack[[1]])[3], ymx=extent(cutStack[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2 <- brick(r,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack)){
  r2[[ii]][which(cutStack[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack[[ii]][which(cutStack[[ii]][] <= Topt)] - Topt)^2))
  r2[[ii]][which(cutStack[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack[[ii]][which(cutStack[[ii]][] > Topt)] - Topt)^2))
}
names(r2) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster
mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2, main = "Current Africa tn RGR of A. pinnata", par.settings = mapTheme) 


##BIOCLIM DATA FOR AFRICA ANNUAL TEMPERATURE
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/Current_bioclim_5min/wc2.1_5m_bio")
annualTemp <- raster("wc2.1_5m_bio_1.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
mapTheme2 <- rasterTheme(region = rev(brewer.pal(9,"RdYlBu")))
cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Current Annual Temperature of Africa", par.settings = mapTheme2) #plot the monthly mean temperatures

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r5 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r5[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r5[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r5, main = "Annual Bio1 Africa RGR of A. pinnata", par.settings = mapTheme) 
Annual_bio_Africa_current_Apinnata <- r5
Annual_bio_Africa_current_Apinnata

##annual filiculoides
annualTemp <- raster("wc2.1_5m_bio_1.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Annual Temperature") #plot the monthly mean temperatures

Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495

r6 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r6[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r6[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r6, main = "Annual Bio1 RGR of A. filiculoides", par.settings = mapTheme) 












########################################################################################################
##############################################################################################################
##############################################################################################################
######################FUTURE PROJECTIONS MONTHLY SENEGAL at 30s reoslution#############################
##############################################################################################################
##############################################################################################################
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/gs26tx50")
dir()

##### Senegal: pull the monthly mean max temperatures for scenario max mean temp RCP 2.6

meanJan_RCP26_tx <- raster("gs26tx501.tif") #mean temp for January so on so forth..
meanFeb_RCP26_tx <- raster("gs26tx502.tif")
meanMar_RCP26_tx <- raster("gs26tx503.tif")
meanApr_RCP26_tx <- raster("gs26tx504.tif")
meanMaj_RCP26_tx <- raster("gs26tx505.tif")
meanJun_RCP26_tx <- raster("gs26tx506.tif")
meanJul_RCP26_tx <- raster("gs26tx507.tif")
meanAvg_RCP26_tx <- raster("gs26tx508.tif")
meanSep_RCP26_tx <- raster("gs26tx509.tif")
meanOkt_RCP26_tx <- raster("gs26tx5010.tif")
meanNov_RCP26_tx <- raster("gs26tx5011.tif")
meanDec_RCP26_tx <- raster("gs26tx5012.tif")

plot(meanJan_RCP26_tx)
##zooming in 
mtStack_RCP26_tx <- stack(meanJan_RCP26_tx, meanFeb_RCP26_tx, meanMar_RCP26_tx, meanApr_RCP26_tx, meanMaj_RCP26_tx, meanJun_RCP26_tx, meanJul_RCP26_tx, meanAvg_RCP26_tx, meanSep_RCP26_tx, meanOkt_RCP26_tx, meanNov_RCP26_tx, meanDec_RCP26_tx)
plot(meanJan_RCP26_tx) #plot the raster file for any month
extentA_RCP26_tx <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP26_tx <- (crop(mtStack_RCP26_tx, extentA_RCP26_tx))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP26_tx) <- month.abb #adds month abbreviations
mapTheme2_RCP26_tx <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP26_tx, main = "Monthly Temperature_RCP26_tx", par.settings = mapTheme2_RCP26_tx ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP26_tx <- raster(nrow= nrow(cutStack_RCP26_tx[[1]]), ncol= ncol(cutStack_RCP26_tx[[1]]), xmn= extent(cutStack_RCP26_tx[[1]])[1], xmx=extent(cutStack_RCP26_tx[[1]])[2],
            ymn=extent(cutStack_RCP26_tx[[1]])[3], ymx=extent(cutStack_RCP26_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP26_tx <- brick(r_RCP26_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP26_tx)){
  r2_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP26_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP26_tx <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP26_tx, main = "RGR A_filiculoides_RCP26_mean max temperatures", par.settings = mapTheme_RCP26_tx) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP26_tx <- raster(nrow= nrow(cutStack_RCP26_tx[[1]]), ncol= ncol(cutStack_RCP26_tx[[1]]), xmn= extent(cutStack_RCP26_tx[[1]])[1], xmx=extent(cutStack_RCP26_tx[[1]])[2],
                  ymn=extent(cutStack_RCP26_tx[[1]])[3], ymx=extent(cutStack_RCP26_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP26_tx <- brick(r_RCP26_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP26_tx)){
  r2_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP26_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP26 <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP26_tx, main = "RGR of A.pinnata_RCP26_mean max temperatures", par.settings = mapTheme_RCP26) 


############################################# Senegal mean min temperatures RCP2.6#####################

## pull the monthly mean max temperatures for scenario RCP 2.6
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/gs26tn50")

meanJan_RCP26_tn <- raster("gs26tn501.tif") #mean temp for January so on so forth..
meanFeb_RCP26_tn <- raster("gs26tn502.tif")
meanMar_RCP26_tn <- raster("gs26tn503.tif")
meanApr_RCP26_tn <- raster("gs26tn504.tif")
meanMaj_RCP26_tn <- raster("gs26tn505.tif")
meanJun_RCP26_tn <- raster("gs26tn506.tif")
meanJul_RCP26_tn <- raster("gs26tn507.tif")
meanAvg_RCP26_tn <- raster("gs26tn508.tif")
meanSep_RCP26_tn <- raster("gs26tn509.tif")
meanOkt_RCP26_tn <- raster("gs26tn5010.tif")
meanNov_RCP26_tn <- raster("gs26tn5011.tif")
meanDec_RCP26_tn <- raster("gs26tn5012.tif")

plot(meanJan_RCP26_tn)
##zooming in 
mtStack_RCP26_tn <- stack(meanJan_RCP26_tn, meanFeb_RCP26_tn, meanMar_RCP26_tn, meanApr_RCP26_tn, meanMaj_RCP26_tn, meanJun_RCP26_tn, meanJul_RCP26_tn, meanAvg_RCP26_tn, meanSep_RCP26_tn, meanOkt_RCP26_tn, meanNov_RCP26_tn, meanDec_RCP26_tn)
plot(meanJan_RCP26_tn) #plot the raster file for any month
extentA_RCP26_tn <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP26_tn <- (crop(mtStack_RCP26_tn, extentA_RCP26_tn))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP26_tn) <- month.abb #adds month abbreviations
mapTheme2_RCP26 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP26_tn, main = "Monthly Temperature_RCP26_min_temp_tn", par.settings = mapTheme2_RCP26 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP26_tn <- raster(nrow= nrow(cutStack_RCP26_tn[[1]]), ncol= ncol(cutStack_RCP26_tn[[1]]), xmn= extent(cutStack_RCP26_tn[[1]])[1], xmx=extent(cutStack_RCP26_tn[[1]])[2],
                     ymn=extent(cutStack_RCP26_tn[[1]])[3], ymx=extent(cutStack_RCP26_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP26_tn <- brick(r_RCP26_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP26_tn)){
  r2_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP26_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP26_tn <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP26_tn, main = "Monthly RGR A.filiculoides_RCP26_tn_min temperatures", par.settings = mapTheme_RCP26_tn) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP26_tn <- raster(nrow= nrow(cutStack_RCP26_tn[[1]]), ncol= ncol(cutStack_RCP26_tn[[1]]), xmn= extent(cutStack_RCP26_tn[[1]])[1], xmx=extent(cutStack_RCP26_tn[[1]])[2],
                     ymn=extent(cutStack_RCP26_tn[[1]])[3], ymx=extent(cutStack_RCP26_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP26_tn <- brick(r_RCP26_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 


for(ii in 1:nlayers(cutStack_RCP26_tn)){
  r2_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP26_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP26 <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP26_tn, main = "Monthly RGR of A.pinnata_RCP26_tn_min temperatures", par.settings = mapTheme_RCP26) 





########################################################################################################
########################## FUTURE PROJECTIONS RCP 4.5#################################################
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/gs45tx50")
dir()

##### Senegal: pull the monthly mean max temperatures for scenario max mean temp RCP 4.5

meanJan_RCP45_tx <- raster("gs45tx501.tif") #mean temp for January so on so forth..
meanFeb_RCP45_tx <- raster("gs45tx502.tif")
meanMar_RCP45_tx <- raster("gs45tx503.tif")
meanApr_RCP45_tx <- raster("gs45tx504.tif")
meanMaj_RCP45_tx <- raster("gs45tx505.tif")
meanJun_RCP45_tx <- raster("gs45tx506.tif")
meanJul_RCP45_tx <- raster("gs45tx507.tif")
meanAvg_RCP45_tx <- raster("gs45tx508.tif")
meanSep_RCP45_tx <- raster("gs45tx509.tif")
meanOkt_RCP45_tx <- raster("gs45tx5010.tif")
meanNov_RCP45_tx <- raster("gs45tx5011.tif")
meanDec_RCP45_tx <- raster("gs45tx5012.tif")

plot(meanJan_RCP45_tx)
##zooming in 
mtStack_RCP45_tx <- stack(meanJan_RCP45_tx, meanFeb_RCP45_tx, meanMar_RCP45_tx, meanApr_RCP45_tx, meanMaj_RCP45_tx, meanJun_RCP45_tx, meanJul_RCP45_tx, meanAvg_RCP45_tx, meanSep_RCP45_tx, meanOkt_RCP45_tx, meanNov_RCP45_tx, meanDec_RCP45_tx)
plot(meanJan_RCP45_tx) #plot the raster file for any month
extentA_RCP45_tx <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP45_tx <- (crop(mtStack_RCP45_tx, extentA_RCP45_tx))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP45_tx) <- month.abb #adds month abbreviations
mapTheme2_RCP45_tx <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP45_tx, main = "Monthly Temperature_RCP45_tx", par.settings = mapTheme2_RCP26_tx ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP45_tx <- raster(nrow= nrow(cutStack_RCP45_tx[[1]]), ncol= ncol(cutStack_RCP45_tx[[1]]), xmn= extent(cutStack_RCP45_tx[[1]])[1], xmx=extent(cutStack_RCP45_tx[[1]])[2],
                     ymn=extent(cutStack_RCP45_tx[[1]])[3], ymx=extent(cutStack_RCP45_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP45_tx <- brick(r_RCP45_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP45_tx)){
  r2_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP45_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP45_tx <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP45_tx, main = "RGR A_filiculoides_RCP45_mean max temperatures", par.settings = mapTheme_RCP45_tx) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP45_tx <- raster(nrow= nrow(cutStack_RCP45_tx[[1]]), ncol= ncol(cutStack_RCP45_tx[[1]]), xmn= extent(cutStack_RCP45_tx[[1]])[1], xmx=extent(cutStack_RCP45_tx[[1]])[2],
                     ymn=extent(cutStack_RCP45_tx[[1]])[3], ymx=extent(cutStack_RCP45_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP45_tx <- brick(r_RCP45_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP45_tx)){
  r2_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP45_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP45_tx <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP26_tx, main = "RGR of A.pinnata_RCP45_maxtemperatures_tx", par.settings = mapTheme_RCP26) 


############################################# Senegal mean min temperatures RCP4.5#####################

## pull the monthly mean max temperatures for scenario RCP 4.5
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/gs45tn50")

meanJan_RCP45_tn <- raster("gs45tn501.tif") #mean temp for January so on so forth..
meanFeb_RCP45_tn <- raster("gs45tn502.tif")
meanMar_RCP45_tn <- raster("gs45tn503.tif")
meanApr_RCP45_tn <- raster("gs45tn504.tif")
meanMaj_RCP45_tn <- raster("gs45tn505.tif")
meanJun_RCP45_tn <- raster("gs45tn506.tif")
meanJul_RCP45_tn <- raster("gs45tn507.tif")
meanAvg_RCP45_tn <- raster("gs45tn508.tif")
meanSep_RCP45_tn <- raster("gs45tn509.tif")
meanOkt_RCP45_tn <- raster("gs45tn5010.tif")
meanNov_RCP45_tn <- raster("gs45tn5011.tif")
meanDec_RCP45_tn <- raster("gs45tn5012.tif")


plot(meanJan_RCP45_tn)
##zooming in 
mtStack_RCP45_tn <- stack(meanJan_RCP45_tn, meanFeb_RCP45_tn, meanMar_RCP45_tn, meanApr_RCP45_tn, meanMaj_RCP45_tn, meanJun_RCP45_tn, meanJul_RCP45_tn, meanAvg_RCP45_tn, meanSep_RCP45_tn, meanOkt_RCP45_tn, meanNov_RCP45_tn, meanDec_RCP45_tn)
plot(meanJan_RCP45_tn) #plot the raster file for any month
extentA_RCP45_tn <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP45_tn <- (crop(mtStack_RCP45_tn, extentA_RCP45_tn))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP45_tn) <- month.abb #adds month abbreviations
mapTheme2_RCP45 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP45_tn, main = "Monthly Temperature_RCP45_min_temp_tn", par.settings = mapTheme2_RCP45 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP45_tn <- raster(nrow= nrow(cutStack_RCP45_tn[[1]]), ncol= ncol(cutStack_RCP45_tn[[1]]), xmn= extent(cutStack_RCP45_tn[[1]])[1], xmx=extent(cutStack_RCP45_tn[[1]])[2],
                     ymn=extent(cutStack_RCP45_tn[[1]])[3], ymx=extent(cutStack_RCP45_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP45_tn <- brick(r_RCP45_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP45_tn)){
  r2_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP45_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP45_tn <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP45_tn, main = "Monthly RGR A.filiculoides_RCP45_min temp tn", par.settings = mapTheme_RCP45_tn) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP45_tn <- raster(nrow= nrow(cutStack_RCP45_tn[[1]]), ncol= ncol(cutStack_RCP45_tn[[1]]), xmn= extent(cutStack_RCP45_tn[[1]])[1], xmx=extent(cutStack_RCP45_tn[[1]])[2],
                     ymn=extent(cutStack_RCP45_tn[[1]])[3], ymx=extent(cutStack_RCP45_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP45_tn <- brick(r_RCP45_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP45_tn)){
  r2_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP45_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP45_tn <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP45_tn, main = "Monthly RGR A.pinnata_RCP45_min temp tn", par.settings = mapTheme_RCP45_tn) 


###### RCP 60########
########################################################################################################
########################## FUTURE PROJECTIONS RCP 6.0#################################################
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/gs60tx50")
dir()

##### Senegal: pull the monthly mean max temperatures for scenario max mean temp RCP 6.0

meanJan_RCP60_tx <- raster("gs60tx501.tif") #mean temp for January so on so forth..
meanFeb_RCP60_tx <- raster("gs60tx502.tif")
meanMar_RCP60_tx <- raster("gs60tx503.tif")
meanApr_RCP60_tx <- raster("gs60tx504.tif")
meanMaj_RCP60_tx <- raster("gs60tx505.tif")
meanJun_RCP60_tx <- raster("gs60tx506.tif")
meanJul_RCP60_tx <- raster("gs60tx507.tif")
meanAvg_RCP60_tx <- raster("gs60tx508.tif")
meanSep_RCP60_tx <- raster("gs60tx509.tif")
meanOkt_RCP60_tx <- raster("gs60tx5010.tif")
meanNov_RCP60_tx <- raster("gs60tx5011.tif")
meanDec_RCP60_tx <- raster("gs60tx5012.tif")

plot(meanJan_RCP60_tx)
##zooming in 
mtStack_RCP60_tx <- stack(meanJan_RCP60_tx, meanFeb_RCP60_tx, meanMar_RCP60_tx, meanApr_RCP60_tx, meanMaj_RCP60_tx, meanJun_RCP60_tx, meanJul_RCP60_tx, meanAvg_RCP60_tx, meanSep_RCP60_tx, meanOkt_RCP60_tx, meanNov_RCP60_tx, meanDec_RCP60_tx)
plot(meanJan_RCP60_tx) #plot the raster file for any month
extentA_RCP60_tx <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP60_tx <- (crop(mtStack_RCP60_tx, extentA_RCP60_tx))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP60_tx) <- month.abb #adds month abbreviations
mapTheme2_RCP60_tx <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP60_tx, main = "Monthly Temperature_RCP60_tx", par.settings = mapTheme2_RCP60_tx ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP60_tx <- raster(nrow= nrow(cutStack_RCP60_tx[[1]]), ncol= ncol(cutStack_RCP60_tx[[1]]), xmn= extent(cutStack_RCP60_tx[[1]])[1], xmx=extent(cutStack_RCP60_tx[[1]])[2],
                     ymn=extent(cutStack_RCP60_tx[[1]])[3], ymx=extent(cutStack_RCP60_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP60_tx <- brick(r_RCP60_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP60_tx)){
  r2_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP60_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP60_tx_fil <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP60_tx, main = "RGR A_filiculoides_RCP60_mean max temp_tx", par.settings = mapTheme_RCP60_tx_fil) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP60_tx <- raster(nrow= nrow(cutStack_RCP60_tx[[1]]), ncol= ncol(cutStack_RCP60_tx[[1]]), xmn= extent(cutStack_RCP60_tx[[1]])[1], xmx=extent(cutStack_RCP60_tx[[1]])[2],
                     ymn=extent(cutStack_RCP60_tx[[1]])[3], ymx=extent(cutStack_RCP60_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP60_tx <- brick(r_RCP60_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP60_tx)){
  r2_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP60_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature


mapTheme_RCP60_tx_pin <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP60_tx, main = "RGR of A.pinnata_RCP60_maxtemperatures_tx", par.settings = mapTheme_RCP60_tx_pin) 


############################################# Senegal mean min temperatures RCP6.0#####################

## pull the monthly mean max temperatures for scenario RCP 6.0
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/gs60tn50")

meanJan_RCP60_tn <- raster("gs60tn501.tif") #mean temp for January so on so forth..
meanFeb_RCP60_tn <- raster("gs60tn502.tif")
meanMar_RCP60_tn <- raster("gs60tn503.tif")
meanApr_RCP60_tn <- raster("gs60tn504.tif")
meanMaj_RCP60_tn <- raster("gs60tn505.tif")
meanJun_RCP60_tn <- raster("gs60tn506.tif")
meanJul_RCP60_tn <- raster("gs60tn507.tif")
meanAvg_RCP60_tn <- raster("gs60tn508.tif")
meanSep_RCP60_tn <- raster("gs60tn509.tif")
meanOkt_RCP60_tn <- raster("gs60tn5010.tif")
meanNov_RCP60_tn <- raster("gs60tn5011.tif")
meanDec_RCP60_tn <- raster("gs60tn5012.tif")


plot(meanJan_RCP60_tn)
##zooming in 
mtStack_RCP60_tn <- stack(meanJan_RCP60_tn, meanFeb_RCP60_tn, meanMar_RCP60_tn, meanApr_RCP60_tn, meanMaj_RCP60_tn, meanJun_RCP60_tn, meanJul_RCP60_tn, meanAvg_RCP60_tn, meanSep_RCP60_tn, meanOkt_RCP60_tn, meanNov_RCP60_tn, meanDec_RCP60_tn)
plot(meanJan_RCP60_tn) #plot the raster file for any month
extentA_RCP60_tn <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP60_tn <- (crop(mtStack_RCP60_tn, extentA_RCP60_tn))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP60_tn) <- month.abb #adds month abbreviations
mapTheme2_RCP60 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP60_tn, main = "Monthly Temperature_RCP60_min_temp_tn", par.settings = mapTheme2_RCP45 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP60_tn <- raster(nrow= nrow(cutStack_RCP60_tn[[1]]), ncol= ncol(cutStack_RCP60_tn[[1]]), xmn= extent(cutStack_RCP60_tn[[1]])[1], xmx=extent(cutStack_RCP60_tn[[1]])[2],
                     ymn=extent(cutStack_RCP60_tn[[1]])[3], ymx=extent(cutStack_RCP60_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP60_tn <- brick(r_RCP60_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP60_tn)){
  r2_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP60_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP60_tn_fil <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP60_tn, main = "Monthly RGR A.filiculoides_RCP60_min temp tn", par.settings = mapTheme_RCP60_tn_fil) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP60_tn <- raster(nrow= nrow(cutStack_RCP60_tn[[1]]), ncol= ncol(cutStack_RCP60_tn[[1]]), xmn= extent(cutStack_RCP60_tn[[1]])[1], xmx=extent(cutStack_RCP60_tn[[1]])[2],
                     ymn=extent(cutStack_RCP60_tn[[1]])[3], ymx=extent(cutStack_RCP60_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP60_tn <- brick(r_RCP60_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP60_tn)){
  r2_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP60_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 
mapTheme_RCP60_tn_pin <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP60_tn, main = "Monthly RGR A.pinnata_RCP60_min temp tn", par.settings = mapTheme_RCP60_tn_pin) 




###### RCP 8.5########
########################################################################################################
########################## FUTURE PROJECTIONS RCP 6.0#################################################
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/gs85tx50")
dir()

##### Senegal: pull the monthly mean max temperatures for scenario max mean temp RCP 8.5

meanJan_RCP85_tx <- raster("gs85tx501.tif") #mean temp for January so on so forth..
meanFeb_RCP85_tx <- raster("gs85tx502.tif")
meanMar_RCP85_tx <- raster("gs85tx503.tif")
meanApr_RCP85_tx <- raster("gs85tx504.tif")
meanMaj_RCP85_tx <- raster("gs85tx505.tif")
meanJun_RCP85_tx <- raster("gs85tx506.tif")
meanJul_RCP85_tx <- raster("gs85tx507.tif")
meanAvg_RCP85_tx <- raster("gs85tx508.tif")
meanSep_RCP85_tx <- raster("gs85tx509.tif")
meanOkt_RCP85_tx <- raster("gs85tx5010.tif")
meanNov_RCP85_tx <- raster("gs85tx5011.tif")
meanDec_RCP85_tx <- raster("gs85tx5012.tif")

plot(meanJan_RCP85_tx)
##zooming in 
mtStack_RCP85_tx <- stack(meanJan_RCP85_tx, meanFeb_RCP85_tx, meanMar_RCP85_tx, meanApr_RCP85_tx, meanMaj_RCP85_tx, meanJun_RCP85_tx, meanJul_RCP85_tx, meanAvg_RCP85_tx, meanSep_RCP85_tx, meanOkt_RCP85_tx, meanNov_RCP85_tx, meanDec_RCP85_tx)
plot(meanJan_RCP85_tx) #plot the raster file for any month
extentA_RCP85_tx <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP85_tx <- (crop(mtStack_RCP85_tx, extentA_RCP85_tx))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP85_tx) <- month.abb #adds month abbreviations
mapTheme2_RCP85_tx <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP85_tx, main = "Monthly Temperature_RCP85_tx", par.settings = mapTheme2_RCP85_tx ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP85_tx <- raster(nrow= nrow(cutStack_RCP85_tx[[1]]), ncol= ncol(cutStack_RCP85_tx[[1]]), xmn= extent(cutStack_RCP85_tx[[1]])[1], xmx=extent(cutStack_RCP85_tx[[1]])[2],
                     ymn=extent(cutStack_RCP85_tx[[1]])[3], ymx=extent(cutStack_RCP85_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP85_tx <- brick(r_RCP85_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP85_tx)){
  r2_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP85_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP85_tx_fil <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP85_tx, main = "RGR A_filiculoides_RCP85_mean max temp_tx", par.settings = mapTheme_RCP85_tx_fil) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP85_tx <- raster(nrow= nrow(cutStack_RCP85_tx[[1]]), ncol= ncol(cutStack_RCP85_tx[[1]]), xmn= extent(cutStack_RCP85_tx[[1]])[1], xmx=extent(cutStack_RCP85_tx[[1]])[2],
                     ymn=extent(cutStack_RCP85_tx[[1]])[3], ymx=extent(cutStack_RCP85_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP85_tx <- brick(r_RCP85_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP85_tx)){
  r2_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP85_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP85_tx_pin <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP85_tx, main = "RGR of A.pinnata_RCP85_maxtemperatures_tx", par.settings = mapTheme_RCP85_tx_pin) 


############################################# Senegal mean min temperatures RCP8.5#####################

## pull the monthly mean max temperatures for scenario RCP 8.5
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/gs85tn50")

meanJan_RCP85_tn <- raster("gs85tn501.tif") #mean temp for January so on so forth..
meanFeb_RCP85_tn <- raster("gs85tn502.tif")
meanMar_RCP85_tn <- raster("gs85tn503.tif")
meanApr_RCP85_tn <- raster("gs85tn504.tif")
meanMaj_RCP85_tn <- raster("gs85tn505.tif")
meanJun_RCP85_tn <- raster("gs85tn506.tif")
meanJul_RCP85_tn <- raster("gs85tn507.tif")
meanAvg_RCP85_tn <- raster("gs85tn508.tif")
meanSep_RCP85_tn <- raster("gs85tn509.tif")
meanOkt_RCP85_tn <- raster("gs85tn5010.tif")
meanNov_RCP85_tn <- raster("gs85tn5011.tif")
meanDec_RCP85_tn <- raster("gs85tn5012.tif")


plot(meanJan_RCP85_tn)
##zooming in 
mtStack_RCP85_tn <- stack(meanJan_RCP85_tn, meanFeb_RCP85_tn, meanMar_RCP85_tn, meanApr_RCP85_tn, meanMaj_RCP85_tn, meanJun_RCP85_tn, meanJul_RCP85_tn, meanAvg_RCP85_tn, meanSep_RCP85_tn, meanOkt_RCP85_tn, meanNov_RCP85_tn, meanDec_RCP85_tn)
plot(meanJan_RCP85_tn) #plot the raster file for any month
extentA_RCP85_tn <- extent(-18,-11,12,16.5) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP85_tn <- (crop(mtStack_RCP85_tn, extentA_RCP85_tn))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP85_tn) <- month.abb #adds month abbreviations
mapTheme2_RCP85 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP85_tn, main = "Monthly Temperature_RCP85_min_temp_tn", par.settings = mapTheme2_RCP85 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP85_tn <- raster(nrow= nrow(cutStack_RCP85_tn[[1]]), ncol= ncol(cutStack_RCP85_tn[[1]]), xmn= extent(cutStack_RCP85_tn[[1]])[1], xmx=extent(cutStack_RCP85_tn[[1]])[2],
                     ymn=extent(cutStack_RCP85_tn[[1]])[3], ymx=extent(cutStack_RCP85_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP85_tn <- brick(r_RCP85_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP85_tn)){
  r2_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP85_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP85_tn_fil <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP85_tn, main = "Monthly RGR A.filiculoides_RCP85_min temp tn", par.settings = mapTheme_RCP85_tn_fil) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP85_tn <- raster(nrow= nrow(cutStack_RCP85_tn[[1]]), ncol= ncol(cutStack_RCP85_tn[[1]]), xmn= extent(cutStack_RCP85_tn[[1]])[1], xmx=extent(cutStack_RCP85_tn[[1]])[2],
                     ymn=extent(cutStack_RCP85_tn[[1]])[3], ymx=extent(cutStack_RCP85_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP85_tn <- brick(r_RCP85_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP85_tn)){
  r2_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP85_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 


mapTheme_RCP85_tn_pin <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP85_tn, main = "Monthly RGR A.pinnata_RCP85_min temp tn", par.settings = mapTheme_RCP85_tn_pin) 




















################################################################################################
############################################################################################################
############################################################################################################
##############AFRICA FUTURE PROJECTIONS for all RCP at  5min resolution !##########################################################################

setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs26tx50")
dir()

##### Africa: pull the monthly mean max temperatures for scenario max mean temp RCP 2.6

meanJan_RCP26_tx <- raster("gs26tx501.tif") #mean temp for January so on so forth..
meanFeb_RCP26_tx <- raster("gs26tx502.tif")
meanMar_RCP26_tx <- raster("gs26tx503.tif")
meanApr_RCP26_tx <- raster("gs26tx504.tif")
meanMaj_RCP26_tx <- raster("gs26tx505.tif")
meanJun_RCP26_tx <- raster("gs26tx506.tif")
meanJul_RCP26_tx <- raster("gs26tx507.tif")
meanAvg_RCP26_tx <- raster("gs26tx508.tif")
meanSep_RCP26_tx <- raster("gs26tx509.tif")
meanOkt_RCP26_tx <- raster("gs26tx5010.tif")
meanNov_RCP26_tx <- raster("gs26tx5011.tif")
meanDec_RCP26_tx <- raster("gs26tx5012.tif")

plot(meanJan_RCP26_tx)
##zooming in 
mtStack_RCP26_tx <- stack(meanJan_RCP26_tx, meanFeb_RCP26_tx, meanMar_RCP26_tx, meanApr_RCP26_tx, meanMaj_RCP26_tx, meanJun_RCP26_tx, meanJul_RCP26_tx, meanAvg_RCP26_tx, meanSep_RCP26_tx, meanOkt_RCP26_tx, meanNov_RCP26_tx, meanDec_RCP26_tx)
plot(meanJan_RCP26_tx) #plot the raster file for any month
extentA_RCP26_tx <- extent(-20,60,-35,40) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP26_tx <- (crop(mtStack_RCP26_tx, extentA_RCP26_tx))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP26_tx) <- month.abb #adds month abbreviations
mapTheme2_RCP26_tx <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP26_tx, main = "Africa Monthly Temperature_RCP26_tx", par.settings = mapTheme2_RCP26_tx ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP26_tx <- raster(nrow= nrow(cutStack_RCP26_tx[[1]]), ncol= ncol(cutStack_RCP26_tx[[1]]), xmn= extent(cutStack_RCP26_tx[[1]])[1], xmx=extent(cutStack_RCP26_tx[[1]])[2],
                     ymn=extent(cutStack_RCP26_tx[[1]])[3], ymx=extent(cutStack_RCP26_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP26_tx <- brick(r_RCP26_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP26_tx)){
  r2_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP26_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP26_tx <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP26_tx, main = "Future Africa RGR A_filiculoides_RCP26_tx", par.settings = mapTheme_RCP26_tx) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP26_tx <- raster(nrow= nrow(cutStack_RCP26_tx[[1]]), ncol= ncol(cutStack_RCP26_tx[[1]]), xmn= extent(cutStack_RCP26_tx[[1]])[1], xmx=extent(cutStack_RCP26_tx[[1]])[2],
                     ymn=extent(cutStack_RCP26_tx[[1]])[3], ymx=extent(cutStack_RCP26_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP26_tx <- brick(r_RCP26_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP26_tx)){
  r2_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP26_tx[[ii]][which(cutStack_RCP26_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP26_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP26_pin <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2_RCP26_tx, main = "Future Africa RGR of A.pinnata_RCP26_tx", par.settings = mapTheme_RCP26_pin) 


############################################# Africa mean min temperatures RCP2.6#####################

## pull the monthly mean max temperatures for scenario RCP 2.6
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs26tn50")
dir()
meanJan_RCP26_tn <- raster("gs26tn501.tif") #mean temp for January so on so forth..
meanFeb_RCP26_tn <- raster("gs26tn502.tif")
meanMar_RCP26_tn <- raster("gs26tn503.tif")
meanApr_RCP26_tn <- raster("gs26tn504.tif")
meanMaj_RCP26_tn <- raster("gs26tn505.tif")
meanJun_RCP26_tn <- raster("gs26tn506.tif")
meanJul_RCP26_tn <- raster("gs26tn507.tif")
meanAvg_RCP26_tn <- raster("gs26tn508.tif")
meanSep_RCP26_tn <- raster("gs26tn509.tif")
meanOkt_RCP26_tn <- raster("gs26tn5010.tif")
meanNov_RCP26_tn <- raster("gs26tn5011.tif")
meanDec_RCP26_tn <- raster("gs26tn5012.tif")

plot(meanJan_RCP26_tn)
##zooming in 
mtStack_RCP26_tn <- stack(meanJan_RCP26_tn, meanFeb_RCP26_tn, meanMar_RCP26_tn, meanApr_RCP26_tn, meanMaj_RCP26_tn, meanJun_RCP26_tn, meanJul_RCP26_tn, meanAvg_RCP26_tn, meanSep_RCP26_tn, meanOkt_RCP26_tn, meanNov_RCP26_tn, meanDec_RCP26_tn)
plot(meanJan_RCP26_tn) #plot the raster file for any month
extentA_RCP26_tn <- extent(-20,60,-35,40) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP26_tn <- (crop(mtStack_RCP26_tn, extentA_RCP26_tn))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP26_tn) <- month.abb #adds month abbreviations
mapTheme2_RCP26 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP26_tn, main = "Africa Monthly Temperature_RCP26_min_temp_tn", par.settings = mapTheme2_RCP26 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP26_tn <- raster(nrow= nrow(cutStack_RCP26_tn[[1]]), ncol= ncol(cutStack_RCP26_tn[[1]]), xmn= extent(cutStack_RCP26_tn[[1]])[1], xmx=extent(cutStack_RCP26_tn[[1]])[2],
                     ymn=extent(cutStack_RCP26_tn[[1]])[3], ymx=extent(cutStack_RCP26_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP26_tn <- brick(r_RCP26_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP26_tn)){
  r2_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP26_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP26_tn <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP26_tn, main = "Future Africa RGR A.filiculoides_RCP26_tn", par.settings = mapTheme_RCP26_tn) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP26_tn <- raster(nrow= nrow(cutStack_RCP26_tn[[1]]), ncol= ncol(cutStack_RCP26_tn[[1]]), xmn= extent(cutStack_RCP26_tn[[1]])[1], xmx=extent(cutStack_RCP26_tn[[1]])[2],
                     ymn=extent(cutStack_RCP26_tn[[1]])[3], ymx=extent(cutStack_RCP26_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP26_tn <- brick(r_RCP26_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 


for(ii in 1:nlayers(cutStack_RCP26_tn)){
  r2_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP26_tn[[ii]][which(cutStack_RCP26_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP26_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP26 <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2_RCP26_tn, main = "Future Africa  RGR of A.pinnata_RCP26_tn", par.settings = mapTheme_RCP26) 


#BIOCLIM DATA AFRICA at 5m resolution for RCP 2.6#
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs26bi50")
dir()
annualTemp <- raster("gs26bi501.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
mapTheme2 <- rasterTheme(region = rev(brewer.pal(9,"RdYlBu")))

cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "FUTURE Annual Temperature of Africa", par.settings = mapTheme2) #plot the monthly mean temperatures

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r5 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r5[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r5[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r5, main = "FUTURE Annual bio1 RGR RCR2.6 of A. pinnata", par.settings = mapTheme) 


##annual filiculoides
annualTemp <- raster("gs26bi501.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Africa Annual Temperature") #plot the monthly mean temperatures

Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495

r6 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 
r6[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r6[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r6, main = "FUTURE Annual bio 1 RGR RCP2.6 of A. filiculoides", par.settings = mapTheme) 






########################################################################################################
########################## FUTURE PROJECTIONS RCP 4.5#################################################

setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs45tx50")


##### Africa: pull the monthly mean max temperatures for scenario max mean temp RCP 4.5

meanJan_RCP45_tx <- raster("gs45tx501.tif") #mean temp for January so on so forth..
meanFeb_RCP45_tx <- raster("gs45tx502.tif")
meanMar_RCP45_tx <- raster("gs45tx503.tif")
meanApr_RCP45_tx <- raster("gs45tx504.tif")
meanMaj_RCP45_tx <- raster("gs45tx505.tif")
meanJun_RCP45_tx <- raster("gs45tx506.tif")
meanJul_RCP45_tx <- raster("gs45tx507.tif")
meanAvg_RCP45_tx <- raster("gs45tx508.tif")
meanSep_RCP45_tx <- raster("gs45tx509.tif")
meanOkt_RCP45_tx <- raster("gs45tx5010.tif")
meanNov_RCP45_tx <- raster("gs45tx5011.tif")
meanDec_RCP45_tx <- raster("gs45tx5012.tif")

plot(meanJan_RCP45_tx)
##zooming in 
mtStack_RCP45_tx <- stack(meanJan_RCP45_tx, meanFeb_RCP45_tx, meanMar_RCP45_tx, meanApr_RCP45_tx, meanMaj_RCP45_tx, meanJun_RCP45_tx, meanJul_RCP45_tx, meanAvg_RCP45_tx, meanSep_RCP45_tx, meanOkt_RCP45_tx, meanNov_RCP45_tx, meanDec_RCP45_tx)
plot(meanJan_RCP45_tx) #plot the raster file for any month
extentA_RCP45_tx <- extent(-20,60,-35,40) #x, y lims 
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP45_tx <- (crop(mtStack_RCP45_tx, extentA_RCP45_tx))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP45_tx) <- month.abb #adds month abbreviations
mapTheme2_RCP45_tx <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP45_tx, main = "Africa Monthly Temperature_RCP45_tx", par.settings = mapTheme2_RCP26_tx ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP45_tx <- raster(nrow= nrow(cutStack_RCP45_tx[[1]]), ncol= ncol(cutStack_RCP45_tx[[1]]), xmn= extent(cutStack_RCP45_tx[[1]])[1], xmx=extent(cutStack_RCP45_tx[[1]])[2],
                     ymn=extent(cutStack_RCP45_tx[[1]])[3], ymx=extent(cutStack_RCP45_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP45_tx <- brick(r_RCP45_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP45_tx)){
  r2_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP45_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP45_tx <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP45_tx, main = "Future Africa RGR A_filiculoides_RCP45_tx", par.settings = mapTheme_RCP45_tx) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP45_tx <- raster(nrow= nrow(cutStack_RCP45_tx[[1]]), ncol= ncol(cutStack_RCP45_tx[[1]]), xmn= extent(cutStack_RCP45_tx[[1]])[1], xmx=extent(cutStack_RCP45_tx[[1]])[2],
                     ymn=extent(cutStack_RCP45_tx[[1]])[3], ymx=extent(cutStack_RCP45_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP45_tx <- brick(r_RCP45_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP45_tx)){
  r2_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP45_tx[[ii]][which(cutStack_RCP45_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP45_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP45_tx <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP26_tx, main = "Future Africa RGR of A.pinnata_RCP45_tx", par.settings = mapTheme_RCP26) 


############################################# Africa mean min temperatures RCP4.5#####################

## pull the monthly mean max temperatures for scenario RCP 4.5
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs45tn50")

meanJan_RCP45_tn <- raster("gs45tn501.tif") #mean temp for January so on so forth..
meanFeb_RCP45_tn <- raster("gs45tn502.tif")
meanMar_RCP45_tn <- raster("gs45tn503.tif")
meanApr_RCP45_tn <- raster("gs45tn504.tif")
meanMaj_RCP45_tn <- raster("gs45tn505.tif")
meanJun_RCP45_tn <- raster("gs45tn506.tif")
meanJul_RCP45_tn <- raster("gs45tn507.tif")
meanAvg_RCP45_tn <- raster("gs45tn508.tif")
meanSep_RCP45_tn <- raster("gs45tn509.tif")
meanOkt_RCP45_tn <- raster("gs45tn5010.tif")
meanNov_RCP45_tn <- raster("gs45tn5011.tif")
meanDec_RCP45_tn <- raster("gs45tn5012.tif")


plot(meanJan_RCP45_tn)
##zooming in 
mtStack_RCP45_tn <- stack(meanJan_RCP45_tn, meanFeb_RCP45_tn, meanMar_RCP45_tn, meanApr_RCP45_tn, meanMaj_RCP45_tn, meanJun_RCP45_tn, meanJul_RCP45_tn, meanAvg_RCP45_tn, meanSep_RCP45_tn, meanOkt_RCP45_tn, meanNov_RCP45_tn, meanDec_RCP45_tn)
plot(meanJan_RCP45_tn) #plot the raster file for any month
extentA_RCP45_tn <- extent(-20,60,-35,40) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP45_tn <- (crop(mtStack_RCP45_tn, extentA_RCP45_tn))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP45_tn) <- month.abb #adds month abbreviations
mapTheme2_RCP45 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP45_tn, main = "Africa Monthly Temperature_RCP45_min_temp_tn", par.settings = mapTheme2_RCP45 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP45_tn <- raster(nrow= nrow(cutStack_RCP45_tn[[1]]), ncol= ncol(cutStack_RCP45_tn[[1]]), xmn= extent(cutStack_RCP45_tn[[1]])[1], xmx=extent(cutStack_RCP45_tn[[1]])[2],
                     ymn=extent(cutStack_RCP45_tn[[1]])[3], ymx=extent(cutStack_RCP45_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP45_tn <- brick(r_RCP45_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP45_tn)){
  r2_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP45_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP45_tn <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2_RCP45_tn, main = "Future Africa  RGR A.filiculoides_RCP45_tn", par.settings = mapTheme_RCP45_tn) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP45_tn <- raster(nrow= nrow(cutStack_RCP45_tn[[1]]), ncol= ncol(cutStack_RCP45_tn[[1]]), xmn= extent(cutStack_RCP45_tn[[1]])[1], xmx=extent(cutStack_RCP45_tn[[1]])[2],
                     ymn=extent(cutStack_RCP45_tn[[1]])[3], ymx=extent(cutStack_RCP45_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP45_tn <- brick(r_RCP45_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP45_tn)){
  r2_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP45_tn[[ii]][which(cutStack_RCP45_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP45_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP45_tn <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP45_tn, main = "Future Africa RGR A.pinnata_RCP45_tn", par.settings = mapTheme_RCP45_tn) 



#BIOCLIM DATA AFRICA at 5m resolution for RCP 4.5#
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs45bi50")
dir()
annualTemp <- raster("gs45bi501.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
mapTheme2 <- rasterTheme(region = rev(brewer.pal(9,"RdYlBu")))

cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "FUTURE Annual Temperature of Africa", par.settings = mapTheme2) #plot the monthly mean temperatures

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r5 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r5[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r5[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r5, main = "FUTURE Africa Annual bio1 RGR RCP 4.5 of A. pinnata", par.settings = mapTheme) 


##annual filiculoides
annualTemp <- raster("gs45bi501.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Africa Annual Temperature") #plot the monthly mean temperatures

Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495

r6 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r6[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r6[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r6, main = "FUTURE Africa Bio 1 Annual RGR RCP 4.5 of A. filiculoides", par.settings = mapTheme) 


###### RCP 60########
########################################################################################################
########################## FUTURE PROJECTIONS RCP 6.0#################################################
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs60tx50")

dir()

##### Africa: pull the monthly mean max temperatures for scenario max mean temp RCP 6.0

meanJan_RCP60_tx <- raster("gs60tx501.tif") #mean temp for January so on so forth..
meanFeb_RCP60_tx <- raster("gs60tx502.tif")
meanMar_RCP60_tx <- raster("gs60tx503.tif")
meanApr_RCP60_tx <- raster("gs60tx504.tif")
meanMaj_RCP60_tx <- raster("gs60tx505.tif")
meanJun_RCP60_tx <- raster("gs60tx506.tif")
meanJul_RCP60_tx <- raster("gs60tx507.tif")
meanAvg_RCP60_tx <- raster("gs60tx508.tif")
meanSep_RCP60_tx <- raster("gs60tx509.tif")
meanOkt_RCP60_tx <- raster("gs60tx5010.tif")
meanNov_RCP60_tx <- raster("gs60tx5011.tif")
meanDec_RCP60_tx <- raster("gs60tx5012.tif")

plot(meanJan_RCP60_tx)
##zooming in 
mtStack_RCP60_tx <- stack(meanJan_RCP60_tx, meanFeb_RCP60_tx, meanMar_RCP60_tx, meanApr_RCP60_tx, meanMaj_RCP60_tx, meanJun_RCP60_tx, meanJul_RCP60_tx, meanAvg_RCP60_tx, meanSep_RCP60_tx, meanOkt_RCP60_tx, meanNov_RCP60_tx, meanDec_RCP60_tx)
plot(meanJan_RCP60_tx) #plot the raster file for any month
extentA_RCP60_tx <- extent(-20,60,-35,40) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP60_tx <- (crop(mtStack_RCP60_tx, extentA_RCP60_tx))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP60_tx) <- month.abb #adds month abbreviations
mapTheme2_RCP60_tx <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP60_tx, main = "Africa Monthly Temperature_RCP60_tx", par.settings = mapTheme2_RCP60_tx ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP60_tx <- raster(nrow= nrow(cutStack_RCP60_tx[[1]]), ncol= ncol(cutStack_RCP60_tx[[1]]), xmn= extent(cutStack_RCP60_tx[[1]])[1], xmx=extent(cutStack_RCP60_tx[[1]])[2],
                     ymn=extent(cutStack_RCP60_tx[[1]])[3], ymx=extent(cutStack_RCP60_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP60_tx <- brick(r_RCP60_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP60_tx)){
  r2_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP60_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP60_tx_fil <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP60_tx, main = "Future Africa RGR A_filiculoides_RCP60_tx", par.settings = mapTheme_RCP60_tx_fil) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP60_tx <- raster(nrow= nrow(cutStack_RCP60_tx[[1]]), ncol= ncol(cutStack_RCP60_tx[[1]]), xmn= extent(cutStack_RCP60_tx[[1]])[1], xmx=extent(cutStack_RCP60_tx[[1]])[2],
                     ymn=extent(cutStack_RCP60_tx[[1]])[3], ymx=extent(cutStack_RCP60_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP60_tx <- brick(r_RCP60_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP60_tx)){
  r2_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP60_tx[[ii]][which(cutStack_RCP60_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP60_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature


mapTheme_RCP60_tx_pin <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP60_tx, main = "Future Africa RGR of A.pinnata_RCP60_tx", par.settings = mapTheme_RCP60_tx_pin) 


############################################# Africa mean min temperatures RCP6.0#####################

## pull the monthly mean max temperatures for scenario RCP 6.0
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs60tn50")

meanJan_RCP60_tn <- raster("gs60tn501.tif") #mean temp for January so on so forth..
meanFeb_RCP60_tn <- raster("gs60tn502.tif")
meanMar_RCP60_tn <- raster("gs60tn503.tif")
meanApr_RCP60_tn <- raster("gs60tn504.tif")
meanMaj_RCP60_tn <- raster("gs60tn505.tif")
meanJun_RCP60_tn <- raster("gs60tn506.tif")
meanJul_RCP60_tn <- raster("gs60tn507.tif")
meanAvg_RCP60_tn <- raster("gs60tn508.tif")
meanSep_RCP60_tn <- raster("gs60tn509.tif")
meanOkt_RCP60_tn <- raster("gs60tn5010.tif")
meanNov_RCP60_tn <- raster("gs60tn5011.tif")
meanDec_RCP60_tn <- raster("gs60tn5012.tif")


plot(meanJan_RCP60_tn)
##zooming in 
mtStack_RCP60_tn <- stack(meanJan_RCP60_tn, meanFeb_RCP60_tn, meanMar_RCP60_tn, meanApr_RCP60_tn, meanMaj_RCP60_tn, meanJun_RCP60_tn, meanJul_RCP60_tn, meanAvg_RCP60_tn, meanSep_RCP60_tn, meanOkt_RCP60_tn, meanNov_RCP60_tn, meanDec_RCP60_tn)
plot(meanJan_RCP60_tn) #plot the raster file for any month
extentA_RCP60_tn <- extent(-20,60,-35,40) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP60_tn <- (crop(mtStack_RCP60_tn, extentA_RCP60_tn))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP60_tn) <- month.abb #adds month abbreviations
mapTheme2_RCP60 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP60_tn, main = "Africa Monthly Temperature_RCP60_min_temp_tn", par.settings = mapTheme2_RCP45 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP60_tn <- raster(nrow= nrow(cutStack_RCP60_tn[[1]]), ncol= ncol(cutStack_RCP60_tn[[1]]), xmn= extent(cutStack_RCP60_tn[[1]])[1], xmx=extent(cutStack_RCP60_tn[[1]])[2],
                     ymn=extent(cutStack_RCP60_tn[[1]])[3], ymx=extent(cutStack_RCP60_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP60_tn <- brick(r_RCP60_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP60_tn)){
  r2_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP60_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP60_tn_fil <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP60_tn, main = "Africa Monthly RGR A.filiculoides_RCP60_tn", par.settings = mapTheme_RCP60_tn_fil) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP60_tn <- raster(nrow= nrow(cutStack_RCP60_tn[[1]]), ncol= ncol(cutStack_RCP60_tn[[1]]), xmn= extent(cutStack_RCP60_tn[[1]])[1], xmx=extent(cutStack_RCP60_tn[[1]])[2],
                     ymn=extent(cutStack_RCP60_tn[[1]])[3], ymx=extent(cutStack_RCP60_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP60_tn <- brick(r_RCP60_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP60_tn)){
  r2_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP60_tn[[ii]][which(cutStack_RCP60_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP60_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 
mapTheme_RCP60_tn_pin <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP60_tn, main = "future Africa  RGR A.pinnata_RCP60_tn", par.settings = mapTheme_RCP60_tn_pin) 


#BIOCLIM DATA AFRICA at 5m resolution for RCP 4.5#
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs60bi50")
dir()
annualTemp <- raster("gs60bi501.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
mapTheme2 <- rasterTheme(region = rev(brewer.pal(9,"RdYlBu")))

cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Current Annual Temperature of Africa", par.settings = mapTheme2) #plot the monthly mean temperatures

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r5 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r5[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r5[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r5, main = "Future Africa Annual bio1 RGR RCP 6 of A. pinnata", par.settings = mapTheme) 


##annual filiculoides
annualTemp <- raster("gs60bi501.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Africa Annual Temperature") #plot the monthly mean temperatures

Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495

r6 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r6[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r6[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r6, main = "Future Africa Annual bio 1 RGR RCP 6 of A. filiculoides", par.settings = mapTheme) 



###### RCP 8.5########
########################################################################################################
########################## FUTURE PROJECTIONS RCP 8.5#################################################
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs85tx50")
dir()

##### Africa: pull the monthly mean max temperatures for scenario max mean temp RCP 8.5

meanJan_RCP85_tx <- raster("gs85tx501.tif") #mean temp for January so on so forth..
meanFeb_RCP85_tx <- raster("gs85tx502.tif")
meanMar_RCP85_tx <- raster("gs85tx503.tif")
meanApr_RCP85_tx <- raster("gs85tx504.tif")
meanMaj_RCP85_tx <- raster("gs85tx505.tif")
meanJun_RCP85_tx <- raster("gs85tx506.tif")
meanJul_RCP85_tx <- raster("gs85tx507.tif")
meanAvg_RCP85_tx <- raster("gs85tx508.tif")
meanSep_RCP85_tx <- raster("gs85tx509.tif")
meanOkt_RCP85_tx <- raster("gs85tx5010.tif")
meanNov_RCP85_tx <- raster("gs85tx5011.tif")
meanDec_RCP85_tx <- raster("gs85tx5012.tif")

plot(meanJan_RCP85_tx)
##zooming in 
mtStack_RCP85_tx <- stack(meanJan_RCP85_tx, meanFeb_RCP85_tx, meanMar_RCP85_tx, meanApr_RCP85_tx, meanMaj_RCP85_tx, meanJun_RCP85_tx, meanJul_RCP85_tx, meanAvg_RCP85_tx, meanSep_RCP85_tx, meanOkt_RCP85_tx, meanNov_RCP85_tx, meanDec_RCP85_tx)
plot(meanJan_RCP85_tx) #plot the raster file for any month
extentA_RCP85_tx <- extent(-20,60,-35,40) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP85_tx <- (crop(mtStack_RCP85_tx, extentA_RCP85_tx))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP85_tx) <- month.abb #adds month abbreviations
mapTheme2_RCP85_tx <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP85_tx, main = "Africa Monthly Temperature_RCP85_tx", par.settings = mapTheme2_RCP85_tx ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP85_tx <- raster(nrow= nrow(cutStack_RCP85_tx[[1]]), ncol= ncol(cutStack_RCP85_tx[[1]]), xmn= extent(cutStack_RCP85_tx[[1]])[1], xmx=extent(cutStack_RCP85_tx[[1]])[2],
                     ymn=extent(cutStack_RCP85_tx[[1]])[3], ymx=extent(cutStack_RCP85_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP85_tx <- brick(r_RCP85_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP85_tx)){
  r2_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP85_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP85_tx_fil <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP85_tx, main = "Future Africa RGR A_filiculoides_RCP85_tx", par.settings = mapTheme_RCP85_tx_fil) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP85_tx <- raster(nrow= nrow(cutStack_RCP85_tx[[1]]), ncol= ncol(cutStack_RCP85_tx[[1]]), xmn= extent(cutStack_RCP85_tx[[1]])[1], xmx=extent(cutStack_RCP85_tx[[1]])[2],
                     ymn=extent(cutStack_RCP85_tx[[1]])[3], ymx=extent(cutStack_RCP85_tx[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP85_tx <- brick(r_RCP85_tx,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP85_tx)){
  r2_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP85_tx[[ii]][which(cutStack_RCP85_tx[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP85_tx) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP85_tx_pin <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP85_tx, main = "Future Africa RGR of A.pinnata_RCP85_tx", par.settings = mapTheme_RCP85_tx_pin) 


############################################# Africa mean min temperatures RCP8.5#####################

## pull the monthly mean max temperatures for scenario RCP 8.5
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs85tn50")

library(raster)

meanJan_RCP85_tn <- raster("gs85tn501.tif") #mean temp for January so on so forth..
meanFeb_RCP85_tn <- raster("gs85tn502.tif")
meanMar_RCP85_tn <- raster("gs85tn503.tif")
meanApr_RCP85_tn <- raster("gs85tn504.tif")
meanMaj_RCP85_tn <- raster("gs85tn505.tif")
meanJun_RCP85_tn <- raster("gs85tn506.tif")
meanJul_RCP85_tn <- raster("gs85tn507.tif")
meanAvg_RCP85_tn <- raster("gs85tn508.tif")
meanSep_RCP85_tn <- raster("gs85tn509.tif")
meanOkt_RCP85_tn <- raster("gs85tn5010.tif")
meanNov_RCP85_tn <- raster("gs85tn5011.tif")
meanDec_RCP85_tn <- raster("gs85tn5012.tif")


plot(meanJan_RCP85_tn)
##zooming in 
mtStack_RCP85_tn <- stack(meanJan_RCP85_tn, meanFeb_RCP85_tn, meanMar_RCP85_tn, meanApr_RCP85_tn, meanMaj_RCP85_tn, meanJun_RCP85_tn, meanJul_RCP85_tn, meanAvg_RCP85_tn, meanSep_RCP85_tn, meanOkt_RCP85_tn, meanNov_RCP85_tn, meanDec_RCP85_tn)
plot(meanJan_RCP85_tn) #plot the raster file for any month
extentA_RCP85_tn <- extent(-20,60,-35,40) #x, y lims
#plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutStack_RCP85_tn <- (crop(mtStack_RCP85_tn, extentA_RCP85_tn))/10 #temperatures are already standard so don't have to divide by 10
#cutStack

names(cutStack_RCP85_tn) <- month.abb #adds month abbreviations
mapTheme2_RCP85 <- rasterTheme(region = (brewer.pal(9,"YlOrRd")))
levelplot(cutStack_RCP85_tn, main = "Future Africa Temperature_RCP85_min_temp_tn", par.settings = mapTheme2_RCP85 ) #plot the monthly mean temperature

##filiculoides monthly
Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495


r_RCP85_tn <- raster(nrow= nrow(cutStack_RCP85_tn[[1]]), ncol= ncol(cutStack_RCP85_tn[[1]]), xmn= extent(cutStack_RCP85_tn[[1]])[1], xmx=extent(cutStack_RCP85_tn[[1]])[2],
                     ymn=extent(cutStack_RCP85_tn[[1]])[3], ymx=extent(cutStack_RCP85_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP85_tn <- brick(r_RCP85_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP85_tn)){
  r2_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP85_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 

mapTheme_RCP85_tn_fil <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP85_tn, main = "Future Africa RGR A.filiculoides_RCP85_tn", par.settings = mapTheme_RCP85_tn_fil) 

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r_RCP85_tn <- raster(nrow= nrow(cutStack_RCP85_tn[[1]]), ncol= ncol(cutStack_RCP85_tn[[1]]), xmn= extent(cutStack_RCP85_tn[[1]])[1], xmx=extent(cutStack_RCP85_tn[[1]])[2],
                     ymn=extent(cutStack_RCP85_tn[[1]])[3], ymx=extent(cutStack_RCP85_tn[[1]])[4]) #one blank layer. Numbers gives  4 values, xmin, xmax, ymin, and ymax, assigning each of those to the new mins/maxes individually
r2_RCP85_tn <- brick(r_RCP85_tn,nl=12) # repeated it 12 times to make a blank brick
##rasterBrick (which is what cutStack is) is  a multi-layered raster object, so  
# cutStack one had 12 layers, one for each months data. Similar to  list can 
#have multiple data frames within in. Similarly, with lists you do [[x]] to get to the specific data frame, 
#then [y,z] to move within that data frame. For bricks, [[x]] gets you to the xth layer, then [y,z] or [y] 
#moves within the raster layer. 

for(ii in 1:nlayers(cutStack_RCP85_tn)){
  r2_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] <= Topt)] <- Pm * exp(-a * ((cutStack_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] <= Topt)] - Topt)^2))
  r2_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] > Topt)] <- Pm * exp(-b * ((cutStack_RCP85_tn[[ii]][which(cutStack_RCP85_tn[[ii]][] > Topt)] - Topt)^2))
}
names(r2_RCP85_tn) <- month.abb #populated the values for each layer based on if the values were < or > the optimum temperature

#the exact indices where the temperature is less than the optimum and using the equation to output the expected growth rate 
#in the exact same place on the empty raster

#install.packages("RColorBrewer") 


mapTheme_RCP85_tn_pin <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r2_RCP85_tn, main = "Future Africa  RGR A.pinnata_RCP85_tn", par.settings = mapTheme_RCP85_tn_pin) 


#BIOCLIM DATA AFRICA at 5m resolution for RCP 4.5#
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs85bi50")
dir()
annualTemp <- raster("gs85bi501.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
mapTheme2 <- rasterTheme(region = rev(brewer.pal(9,"RdYlBu")))

cutAnnual <- (crop(annualTemp, extentA)/10) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Future Annual Temperature of Africa", par.settings = mapTheme2) #plot the monthly mean temperatures

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r5 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r5[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r5[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
                    
levelplot(r5, main = "Future Africa Annual bio 1 RGR RCP 8.5 of A. pinnata", par.settings = mapTheme) 
#writeRaster(r6, filename = "Africa_annual_8.5_pinnata.tif", format= "GTiff") #convert map to Raster file


##annual filiculoides
annualTemp <- raster("gs85bi501.tif")
plot (annualTemp)
extentA <- extent(-20,60,-35,40)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutAnnual <- (crop(annualTemp, extentA)/10) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Future Africa Annual Temperature") #plot the monthly mean temperatures

Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495

r6 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r6[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r6[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r6, main = "Future Africa Annual bio 1 RGR RCP 8.5 of A. filiculoides", par.settings = mapTheme) 

Africa_annual_8.5_pinnata <- r6   

###################
difference <- (Africa_annual_8.5_pinnata - Annual_bio_Africa_current_Apinnata, main = "difference in pinnata" par.settings = mapTheme) ##play with the ylimits #put 3 colors and middle color should be white. 
#getting an averagge of africa and the regions. 
summary(difference)
levelplot(difference)

str(Annual_bio_Africa_current_Apinnata)
