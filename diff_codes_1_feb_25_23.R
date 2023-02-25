library(raster)
library(rasterVis)
library(RColorBrewer)


##BIOCLIM DATA FOR Current AFRICA ANNUAL TEMPERATURE
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/Current_bioclim_5min/wc2.1_5m_bio")
annualTemp <- raster("wc2.1_5m_bio_1.tif")
plot (annualTemp)
extentA <- extent(-20,52,-35,37.5)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
mapTheme2 <- rasterTheme(region = rev(brewer.pal(9,"RdYlBu")))
cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10

setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/manuscript/PLOS 1/submission 2/Figures_Tiff_files")
tiff("S1.1_r.tiff", res = 600, width = 10, height = 8, units = "in")
levelplot(cutAnnual, main = "Current Annual Temperature of Africa", par.settings = mapTheme2, margin = FALSE) #plot the monthly mean temperatures
dev.off()


#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r5 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r5[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r5[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
tiff("S1.2_r.tiff", res = 600, width = 10, height = 8, units = "in")
levelplot(r5, main = "Current Annual Bio1 Africa RGR of A. pinnata", par.settings = mapTheme, margin = FALSE) 
Current_Africa_annual_pinnata <- r5  
dev.off()


##annual filiculoides
annualTemp <- raster("wc2.1_5m_bio_1.tif")
plot (annualTemp)
extentA <- extent(-20,52,-35,37.5)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutAnnual <- (crop(annualTemp, extentA)) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Current Annual Temperature", margin = FALSE) #plot the monthly mean temperatures

Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495

r6 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r6[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r6[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/manuscript/PLOS 1/submission 2/Figures_Tiff_files")
tiff("S1.3_r.tiff", res = 600, width = 10, height = 8, units = "in")
levelplot(r6, main = "Current Annual Bio1 RGR of A. filiculoides", par.settings = mapTheme, margin = FALSE) 
dev.off()
Current_Africa_annual_filiculoides <- r6   




########FUTURE##########



#BIOCLIM DATA AFRICA at 5m resolution for RCP 4.5#
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs45bi50")
dir()
annualTemp <- raster("gs45bi501.tif")
plot (annualTemp)
extentA <- extent(-20,52,-35,37.5)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
mapTheme2 <- rasterTheme(region = rev(brewer.pal(9,"RdYlBu")))

cutAnnual <- (crop(annualTemp, extentA)/10) #temperatures are already standard so don't have to divide by 10
tiff("S2.300_r_future Annual temp of Africa.tiff", res = 300, width = 10, height = 8, units = "in")
levelplot(cutAnnual, main = "FUTURE Annual Temperature of Africa", par.settings = mapTheme2, margin = FALSE) #plot the monthly mean temperatures
dev.off()

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r1 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r1[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r1[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r1, main = "FUTURE Africa Annual bio1 RGR RCP 4.5 of A. pinnata", par.settings = mapTheme, margin = FALSE) 
Africa_future_annual_4.5_pinnata <- r1 

##annual filiculoides
annualTemp <- raster("gs45bi501.tif")
plot (annualTemp)
extentA <- extent(-20,52,-35,37.5)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutAnnual <- (crop(annualTemp, extentA)/10) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Africa Annual Temperature", margins= FALSE) #plot the monthly mean temperatures

Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495

r2 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r2[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r2[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r2, main = "FUTURE Africa Bio 1 Annual RGR RCP 4.5 of A. filiculoides", par.settings = mapTheme, margin = FALSE) 
Africa_future_annual_4.5_filiculoides <- r2 







########BIOCLIM DATA FOR Future 8.5 AFRICA ANNUAL TEMPERATURE
#BIOCLIM DATA AFRICA at 5m resolution for RCP 8.50Buisiness as usual-pessamistic #
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/future projections/GS_5min_RCP_futureproj/gs85bi50")
dir()
annualTemp <- raster("gs85bi501.tif")
plot (annualTemp)
extentA <- extent(-20,52,-35,37.5)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
mapTheme2 <- rasterTheme(region = rev(brewer.pal(9,"RdYlBu")))

cutAnnual <- (crop(annualTemp, extentA)/10) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Future Annual Temperature of Africa", par.settings = mapTheme2, margin = FALSE) #plot the monthly mean temperatures

#input Model pinnata
Pm <- 0.252 #inputting parameters
Topt <- 24.5
a <- 0.0155
b <- 0.00344

r3 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r3[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r3[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))

levelplot(r3, main = "Future Africa Annual bio 1 RGR RCP 8.5 of A. pinnata", par.settings = mapTheme, margin = FALSE) 
#writeRaster(r6, filename = "Africa_annual_8.5_pinnata.tif", format= "GTiff") #convert map to Raster file
Africa_future_annual_8.5_pinnata <- r3

##annual filiculoides
annualTemp <- raster("gs85bi501.tif")
plot (annualTemp)
extentA <- extent(-20,52,-35,37.5)
plot(zoom(annualTemp, extentA), add = TRUE)

##cutting, stack and drawing level plot
cutAnnual <- (crop(annualTemp, extentA)/10) #temperatures are already standard so don't have to divide by 10
levelplot(cutAnnual, main = "Future Africa Annual Temperature") #plot the monthly mean temperatures

Pm <- 0.184 #inputing parameters
Topt <- 24.5
a <- 0.00501
b <- 0.00495

r4 <- raster(nrow= nrow(cutAnnual), ncol= ncol(cutAnnual), xmn= extent(cutAnnual)[1], xmx=extent(cutAnnual)[2],ymn=extent(cutAnnual)[3], ymx=extent(cutAnnual)[4]) 

r4[which(cutAnnual[] <= Topt)] <- Pm * exp(-a * ((cutAnnual[which(cutAnnual[] <= Topt)] - Topt)^2))
r4[which(cutAnnual[] > Topt)] <- Pm * exp(-b * ((cutAnnual[which(cutAnnual[] > Topt)] - Topt)^2))

mapTheme <- rasterTheme(region = (brewer.pal(9,"YlGn")))
levelplot(r4, main = "Future Africa Annual bio 1 RGR RCP 8.5 of A. filiculoides", par.settings = mapTheme, margin = FALSE) 
Africa_future_annual_8.5_filiculoides <- r4   

##compiling future maps
setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/manuscript/PLOS 1/submission 2/Figures_Tiff_files")

library(gridExtra)
tiff("S2_fig2.tiff", res = 600, width = 10, height = 8, units = "in")
par(mfrow = c(2,2))
b<- levelplot(r1, main = "B                                                         ", par.settings = mapTheme, margin = FALSE) 
c <- levelplot(r2, main = "C                                                         ", par.settings = mapTheme, margin = FALSE) 
d <- levelplot(r3, main = "D                                                         ", par.settings = mapTheme, margin = FALSE) 
e <-  levelplot(r4, main = "E                                                         ", par.settings = mapTheme, margin = FALSE) 
#tiff("Fig3_r", res = 600)
grid.arrange (b,c,d,e)
dev.off()





####DIfferences###

difference_pin_4.5 <- Africa_future_annual_4.5_pinnata - Current_Africa_annual_pinnata
difference_pin_8.5 <- Africa_future_annual_8.5_pinnata - Current_Africa_annual_pinnata

difference_fil_4.5 <- Africa_future_annual_4.5_filiculoides - Current_Africa_annual_filiculoides
difference_fil_8.5 <- Africa_future_annual_8.5_filiculoides - Current_Africa_annual_filiculoides

mapTheme_diff <- rasterTheme(region = (brewer.pal(11,"RdBu"))) #make theme
my.at <- seq(-0.07, 0.07, 0.01)# make legends
my.at2 <- seq(-0.07, 0.07, 0.01)
#my.at2 <- seq(-0.184, 0.184, 0.01)
levelplot(difference_pin_4.5, main = "A. pinnata under the RCP 4.5 scenario", at=my.at, par.settings = mapTheme_diff, margin = FALSE) 
#install.packages("grid")
#library(grid)
#grid.text(expression(m^3/m^3), 0.2, 0, hjust=0.5, vjust=1)
levelplot(difference_pin_8.5, main = "A. pinnata under the RCP 8.5 scenario", at=my.at, par.settings = mapTheme_diff, margin = FALSE, cex.axis=0.5, cex.lab=0.5) 
levelplot(difference_fil_4.5, main = "A. filiculoides under the RCP 4.5 scenario", at=my.at2, par.settings = mapTheme_diff, margin = FALSE) 
levelplot(difference_fil_8.5, main = "A. filiculoides under the RCP 8.5 scenario", at=my.at2, par.settings = mapTheme_diff, margin = FALSE) 

setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/manuscript/PLOS 1/submission 2/Figures_Tiff_files")

library(gridExtra)
tiff("Fig3_r.tiff", res = 600, width = 10, height = 8, units = "in")
par(mfrow = c(2,2))
w<- levelplot(difference_pin_4.5, main = "A                                   ", at=my.at, par.settings = mapTheme_diff, margin = FALSE) 
x <- levelplot(difference_pin_8.5, main = "B                                   ", at=my.at, par.settings = mapTheme_diff, margin = FALSE) 
y <- levelplot(difference_fil_4.5, main = "C                                   ", at=my.at2, par.settings = mapTheme_diff, margin = FALSE) 
z <- levelplot(difference_fil_8.5, main = "D                                   ", at=my.at2, par.settings = mapTheme_diff, margin = FALSE) 
#tiff("Fig3_r", res = 600)
grid.arrange (w,x,y,z)
dev.off()



##summary statistics for Current Africa
mean_current_pin <- cellStats(Current_Africa_annual_pinnata, mean) #0.2204475
mean_current_fil <- cellStats(Current_Africa_annual_filiculoides, mean) #0.171971

##summary statistics for Future Africa
mean_future_pin_4.5 <- cellStats(Africa_future_annual_4.5_pinnata,mean) #0.2278979
mean_future_pin_8.5 <-cellStats(Africa_future_annual_8.5_pinnata, mean) # 0.2288697
mean_future_fil_4.5 <- cellStats(Africa_future_annual_4.5_filiculoides, mean) #0.1713162
mean_future_fil_8.5 <- cellStats(Africa_future_annual_8.5_filiculoides, mean) #0.1707884


#summary statistics for RGR change
(mean_future_pin_4.5/mean_current_pin *100) - 100 #3.37967 %
(mean_future_pin_8.5/mean_current_pin *100) - 100 #3.820511 %
(mean_future_fil_4.5/mean_current_fil *100) - 100 #-0.3807632 %
(mean_future_fil_8.5/mean_current_fil *100) - 100# -0.6876396 %




mean_future_pin_4.5<- cellStats(difference_pin_4.5, mean) # 0.007421581 #cellStats if, instead, you want to obtain a summary for all cells of a single Raster object
cellStats(difference_pin_4.5, min) # -0.07341683
cellStats(difference_pin_4.5, max) # 0.1438844

mean_future_pin_8.5<- cellStats(difference_pin_8.5, mean) # 0.008393325
cellStats(difference_pin_8.5, min) #-0.06527796
cellStats(difference_pin_8.5, max) #0.15373


mean_future_fil_4.5<-cellStats(difference_fil_4.5, mean) #-0.0006630896
cellStats(difference_fil_4.5, min) #-0.04237284
cellStats(difference_fil_4.5, max) #0.05513567


mean_future_fil_8.5<-cellStats(difference_fil_8.5, mean)#-0.001190861
cellStats(difference_fil_8.5, min) #-0.04551854
cellStats(difference_fil_8.5, max) #0.05894328

### percentage change of means from current to future


mean_future_pin_4.5/mean_current_pin *100 #3.366598
mean_future_pin_8.5/mean_current_pin *100 #3.807403

mean_future_fil_4.5/mean_current_fil *100 #-0.3855823
mean_future_fil_8.5/mean_current_fil *100 # -0.6924776


percent_difference_pin_4.5 <- Africa_future_annual_4.5_pinnata / Current_Africa_annual_pinnata
percent_difference_pin_8.5 <- Africa_future_annual_8.5_pinnata - Current_Africa_annual_pinnata

percent_difference_fil_4.5 <- Africa_future_annual_4.5_filiculoides - Current_Africa_annual_filiculoides
percent_difference_fil_8.5 <- Africa_future_annual_8.5_filiculoides - Current_Africa_annual_filiculoides



#setwd("~/Documents/EMORY_PHD/Year 3/Aim 1/Africa-wide project/Bioclim-temperature/Current_bioclim_5min/wc2.1_5m_bio")
#annualTemp <- raster("wc2.1_5m_bio_1.tif")
#plot (annualTemp)

#extentA <- getData('alt',country=c('EGY','SDN','ERI',level=0))
#plot(zoom(annualTemp, extentA), add = TRUE)


