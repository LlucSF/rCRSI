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
#' @param path full path to .txt file.
#' @param numBands Number of .txt files to process.
#'
#' @description Reads a .txt file containing Raman data from Reinshaw devices and convert it to an RamanR data object.
#'
#' @return a rCRSIObj data object.
#' @export
#'

DataInput.txt <- function(path = NULL, numBands = 1)
{
  writeLines("Select files...")

  #Data input read
  PathToFile <- list()
  for(i in 1:(numBands))
  {
    if(is.null(path))
      {
        PathToFile[[i]] <- file.choose()
      }else
        {
          PathToFile[[i]] <- path[i]
        }
  }


  writeLines("Starting the data structuring process...")

  rawData <- list()
  for(i in 1:(numBands))
  {
    tmp <- utils::read.table(file = PathToFile[[i]],header = FALSE)
    if(i==1)
    {
      rawData <- tmp
    } else
     {
       rawData <- cbind(rawData,tmp[,3:4])
     }
  }

  #Raman Shift Axis, number of spectra & Raman band length
  RamanShift <- array(dim = c(length(union(rawData[[3]],rawData[[3]])),numBands))
  for(i in 1:numBands)
  {
    p1 <- 3+(i-1)*2
    RamanShift[,i] <- union(rawData[[p1]],rawData[[p1]])
  }
  BandLength <- (length(RamanShift)/numBands)

  #Coordenates
  X <- union(rawData[[1]],rawData[[1]])
  Y <- union(rawData[[2]],rawData[[2]])
  numPixels <- length(X)*length(Y)
  Coords <- list(X = X, Y = Y)

  #Data Shaping
  RamanData <- array(dim=c(numBands*numPixels,BandLength))
  for(i in 1:numBands)
  {
    p <- 2*i+2
    for(j in 1:numPixels)
    {
      RamanData[j+numPixels*(i-1),] <- rawData[[p]][((j-1)*BandLength+1):(j*BandLength)]
    }
  }

  writeLines("Data structuring complete...")
  return(rCRMIObj(RamanData, numPixels, numBands, RamanShift, BandLength, Coords))
}


rCRMIObj <- function(RamanData, numPixels, numBands = 1, RamanShift, BandsLenght, Coords)
{
  rCRSIObj <- list()

  rCRSIObj$numPixels    <- numPixels
  rCRSIObj$numBands     <- numBands
  rCRSIObj$numSpectr    <- numPixels*numBands
  rCRSIObj$BandsLength  <- BandsLenght
  rCRSIObj$AbsCoords    <- Coords
  rCRSIObj$RelCoords    <- list(X = 1:length(Coords$X), Y = 1:length(Coords$Y))
  rCRSIObj$PixelSize    <- c(Coords$X[2]-Coords$X[1], Coords$Y[2]-Coords$Y[1])
  rCRSIObj$Data         <- RamanData

  # ID <- gsub(x = as.character(Sys.time()),pattern = ("-"),replacement = "")
  # ID <- gsub(x = ID,pattern = (":"),replacement = "")
  # ID <- gsub(x = ID,pattern = (" "),replacement = "")
  # rCRSIObj$ID           <- ID
  #
  # if (numBands == 1 & BandsLenght > 391)
  #   {
  #   rCRSIObj$FullSpect <- TRUE
  #   }
  #    else   rCRSIObj$FullSpect <- FALSE

  rCRSIObj$RamanShiftAxis <- RamanShift
  rCRSIObj$AvrgSpectr     <- RamanShift

  for (i in 1:numBands)
  {
    for (j in 1:numPixels)
    {
      rCRSIObj$AvrgSpectr[,i] <- 0
    }
  }
  rCRSIObj$AvrgSpectr  <- AverageSpectrum(numBands, numPixels, rCRSIObj)
  return (rCRSIObj)
}

