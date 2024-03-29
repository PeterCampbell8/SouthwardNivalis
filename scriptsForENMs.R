if (!require("remotes")) {
  install.packages("remotes")
}
install_cran("knitr")
remotes::install_github("fmachados/grinnell")
library(grinnell)
library(raster)
library(terra)
library(spThin)
library(kuenm)


#-------Points Thinning-------- #General format used for all point thinning
inPoints <- read.csv("PostSixtiesLeastWeasels.csv")

occt <- thin(loc.data = inPoints, lat.col = "latitude", long.col = "longitude",   
             spec.col = "species", thin.par = 30, reps = 5, 
             locs.thinned.list.return = FALSE, write.files = TRUE,
             max.files = 1, out.dir = "Nivalis/Post60s/",
             out.base = "PostSixtiesSP30km")


#----Write TerraClim data to ASCIIs for creating the M simulations------
rastlist <- list.files(path = "./EnviroData", pattern='.tif$', 
                       all.files=TRUE, full.names=TRUE)

#import all raster files in folder using lapply
allrasters <- stack(rastlist)

#tif to ascii
raster::writeRaster(allrasters, filename = "rasterStack.asc", format = "ascii",
                    bylayer = TRUE)

#----------Run MSims------------
#Thin the cleaned M.nivalis records using the above thinning function -> thinnedPoints.csv

setwd("PATH TO YOUR WORKING DIRECTORY")
# Set your Python path
py.path <- "PATH TO YOUR PYTHON COMPILER"

# Call Python
Sys.setenv(PATH = paste(py.path, Sys.getenv()["PATH"], sep = ";"))

system("python --version")

tempdir <- file.path(getwd(), "msim")

#Define the Msim variables as a list of desired to test
kernel_variables <- c(5,3,1)
dispers_variables <- c(5,4)
events <- c(125, 250)

#Run the simulations, outputting each combination of variables as its own folder
for (x in kernel_variables) {
  for (y in dispers_variables) {
    for (z in events){
      outName = paste0(tempdir, sprintf("/eg_msim_N_%s_%s_%s",x,y,z)) 
      
      M_simulation(data = "thinnedPoints.csv", current_variables = "DIRECTORY WITH ENVIROMENTAL ASCIIS",
                   dispersal_kernel = "normal", kernel_spread = x,
                   max_dispersers = y, dispersal_events = z,
                   output_directory = outName)
    }
  }
}



#------Layering and averaging all environmental rasters, outputting an ASC-------
#Set the folder that holds the desired data (e.g., present precipitation)
rastFold <- "PPrecipAll"

#Make a SpatialRaster from all bil files in said folder
rastStack <- rast(list.files(path=rastFold, pattern=".bil$", full.names = TRUE))

#All months aggreg down to a mean at each pixel
set <- mean(rastStack)

#write raster to a file
writeRaster(set, "presentPrecipAgg.asc")

#Repeat for all environment variable groups




#------Craete PCAs and project Present data to Historic PCA loadings---------
historicFolder <- "Folder with the pre-60s time period averaged variables"
presentFolder <- "Folder with the post-60s time period averaged variables"


pca_raster(historicFolder, in_format = "ascii", project = TRUE, projection_variables = presentFolder, 
           write_to_directory = TRUE, out_format = "ascii",
           output_directory = "HistoricPCAPresentData")



#---Thin the original least weasel records within the US to 30km-----


#-------Environment Bias Layer Creation--------------
backTable <- read.csv("RodentPlusMust.csv") #Target-group background points with lat and long columns
thinTable <- backTable[,c('decimalLon','decimalLat')]

modelRast <- rast("PC1.asc")

#create a blank raster that has the same metadata as one of the environmental rasters
blankRaster <- mask(modelRast, modelRast, maskvalues=1)

#Turn the csv points into a SpatVector
counted <- extract(blankRaster, thinTable, cells=TRUE, bind = TRUE)

#Create a raster with counts of points per cell, over the model raster (environment). Background has to be non-zero for maxent, so it's given a tiny value (.001)
newRast <- rasterize(counted, modelRast, fun = function(i) {length(i)}, background = .001)

#Exclude non- m layer area. Write to file.
maskedRast <- mask(newRast, modelRast)
writeRaster(maskedRast, "RodentPlusMust.asc", overwrite = TRUE)

#Set the values as a percentage of the sum of all cells, 273291
denseRast <- maskedRast / 273291 * 1000

logRast <- log(maskedRast) / log(273291) #just for visualization

#Reassign the novalue to match the environments and write file
writeRaster(denseRast, "backgroundDensity.asc", overwrite = TRUE, NAflag = -9999)

#----------------------------------------------------------



#create the built-in markdown file to create ENM calibration models, candidate models, final models, post-60s projected models, jackknife analysis, and mmop analysis
kuenm_start('nivalis_enm_process')



#---------------Binarize rasters----------------- Threshold value = .307 (identified from sampled points in QGIS)
oldRast <- rast("calibration_med.tif")
newRast <- rast("Gset_1_med.tif")

f <- function(i) ifelse(i<.307,0,1)

newoldRast <- app(oldRast, f)
newnewRast <- app(newRast, f)

writeRaster(newoldRast, "BiasOldThreshold.tif")
writeRaster(newnewRast, "BiasNewThreshold.tif")



#------Histogram of differences in suitability between time periods
dat <- read.csv("ModernSuitMedianMap.csv") #CSV with spatially thinned modern weasels, their pre-60s and post-60s suitability, and the difference between the two.
dat <- dat[5:6]

means <- colMeans(dat)

mdif <- diff(means)

difs <- sapply(1:100000, function(x) {
  diff(colMeans(cbind(sample(unlist(dat), nrow(dat)), sample(unlist(dat), nrow(dat)))))
})

hist(difs, breaks = 50, xlim = c(-.12, .12), xlab="Difference", main = "Histogram of Differences")
abline(v = mdif, lty = "dashed")
abline(v = 2*sd(difs), col = "blue")
abline(v = 2*-sd(difs), col = "blue")


mdif1 <- mean(apply(dat, 1, diff))


difs1 <- sapply(1:100000, function(x) {
  mean(apply(cbind(sample(unlist(dat), nrow(dat)), sample(unlist(dat), nrow(dat))), 1, diff))
})

hist(difs1, breaks = 50)
abline(v = mdif1, col = "red")