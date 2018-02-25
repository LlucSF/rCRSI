#########################################################################
#     rCRMI - R package for Confocal Raman Spectroscopy Imaging data processing and visualization
#     Copyright (C) 2018 Lluc Sementé Fernàndez
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################

#' DataInput.txt
#'
#' @param PathtoFile full path to .txt file.
#'
#' Reads a .txt file containing Raman data from Reinshaw devices and convert it to an RamanR data object.
#'
#' @return a RamanR data object.
#' @export
#'

DataInput.txt <- function(numBands = 1)
{
  writeLines("Select files...")

  #Data input read
  PathToFile <- list()
  for(i in 1:(numBands))
  {
    PathToFile[[i]]<-file.choose()
  }

  rawData <- list()
  for(i in 1:(numBands))
  {
    tmp <- read.table(file = PathToFile[[i]],header = FALSE)
    if(i==1)
    {
      rawData <- tmp
    } else
     {
       rawData <- cbind(rawData,tmp[,3:4])
     }
  }
  writeLines("Starting the data structuring process...")

  #Raman Shift Axis, number of spectra & Raman band length
  RamanShift <- array(dim = c(length(union(rawData[[3]],rawData[[3]])),numBands))
  for(i in 1:numBands)
  {
    p1 <- 3+(i-1)*2
    RamanShift[,i] <- union(rawData[[p1]],rawData[[p1]])
  }
  BandLength <- (length(RamanShift)/numBands)

  #Coordenates
  Y <- union(rawData[[1]],rawData[[1]])
  X <- union(rawData[[2]],rawData[[2]])
  numPixels <- length(X)*length(Y)
  Coords <- list(X,Y)

  #Data Shaping
  RamanData <- array(dim=c(BandLength,numBands*numPixels))
  for(i in 1:numBands)
  {
    p <- 2*i+2
    for(j in 1:numPixels)
    {
      RamanData[,j+numPixels*(i-1)] <- rawData[[p]][((j-1)*BandLength+1):(j*BandLength)]
    }
  }

  writeLines("Data structuring complete...")
  return(rCRMIObj(RamanData, numPixels, numBands, RamanShift, BandLength,Coords))
}

#' rCRMIObj
#'
#' Creates a rCRMI object .
#'
#' @return a blank rCRMI data object.
#'
#'

rCRMIObj <- function(RamanData, numPixels, numBands = 1, RamanShift, BandsLenght, Coords)
{
  RmnObj <- list()

  RmnObj$numPixels    <- numPixels
  RmnObj$numBands     <- numBands
  RmnObj$numSpectr    <- numPixels*numBands
  RmnObj$BandsLength  <- BandsLenght
  RmnObj$AbsCoords    <- Coords
  RmnObj$RelCoords    <- list(trunc(Coords[[1]]-min(Coords[[1]])),trunc(Coords[[2]]-min(Coords[[2]])))
  RmnObj$Data         <- RamanData

  if (numBands == 1 & BandsLenght > 391)
    {
    RmnObj$FullSpect <- TRUE
    }
    else   RmnObj$FullSpect <- FALSE

  RmnObj$RamanShiftAxis <- RamanShift
  RmnObj$AvrgSpectr      <- RamanShift

  for (i in 1:numBands)
  {
    for (j in 1:numPixels)
    {
      RmnObj$AvrgSpectr[,i] <- 0
    }
  }
  RmnObj$AvrgSpectr  <- AverageSpectrum(numBands, numPixels, RmnObj)

  for (i in 1:numBands)
  {
  g <- ggplot2::qplot(y = RmnObj$AvrgSpectr[,i], x = RmnObj$RamanShiftAxis[,i],
                      ylim = c(min(RmnObj$AvrgSpectr[,i]), max(RmnObj$AvrgSpectr[,i])),
                      main = "Average Spectrum",geom ="line",
                      xlab = "Raman Shift (WaveNum/cm)",ylab = "Counts",color="red")
  print(g)
  }


  return (RmnObj)
}











