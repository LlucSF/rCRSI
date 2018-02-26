#########################################################################
#     rCRMI - R package for Confocal Raman Spectroscopy Imaging RmnObj processing and visualization
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

#' SmoothSpectra
#'
#' @param RmnObj a RmnObj.
#' @param ZroThreshold  Boolean. Apply negative-to-zero thresholding during the wavelets preprocessing.
#' @param BackgroundSubs Boolean. Apply background substraction using a low frequency filter wavelets processing.
#'
#' @description  Smooths the data contained in a RmnObj using the discrete wavelets transform (DWT).
#'
#' @return a RmnObj object.
#' @export
#'

SmoothSpectra <- function(RmnObj, ZroThreshold = TRUE,
                             BackgroundSubs = TRUE,window = 32, ...)
{

  writeLines("Starting wavelets processing...")
  start_time <- Sys.time()
  #window <-32
  #RmnObj <- RmnObj

  for(a in 1:RmnObj$numBands)
  {
    for(b in 1:RmnObj$numPixels)
    {
      spectrum <- c(rep(x = RmnObj$Data[1,(b+((a-1)*RmnObj$BandsLength))],times = window),
                    RmnObj$Data[,(b+((a-1)*RmnObj$BandsLength))],
                    rep(x = RmnObj$Data[RmnObj$BandsLength,(b+((a-1)*RmnObj$BandsLength))],times = window))

      Wav_c <- wavelets::modwt(X = spectrum, n.levels = 4, filter = "d4",
                               boundary = "periodic", fast = TRUE)

      W<-Wav_c@W
      W$W1<-W$W1-W$W1
      W$W2<-W$W2-W$W2
      W$W3<-W$W3-W$W3
      W$W4<-W$W4-W$W4
      Wav_c@W<-W

      RmnObj$Data[,b+(a-1)*RmnObj$BandsLength]<-imodwt(Wav_c)[(window+1):(window+RmnObj$BandsLength)]

      if (BackgroundSubs == TRUE)
      {
      spectrum <- c(rep(x = RmnObj$Data[1,(b+((a-1)*RmnObj$BandsLength))],times = window),
                    RmnObj$Data[,(b+((a-1)*RmnObj$BandsLength))],
                    rep(x = RmnObj$Data[RmnObj$BandsLength,(b+((a-1)*RmnObj$BandsLength))],times = window))

      Wav_c <- wavelets::modwt(X = spectrum, n.levels=7, filter="d6",
                               boundary="periodic", fast=TRUE)

      W<-Wav_c@W
      W$W7<-W$W7-W$W7
      Wav_c@W<-W
      V<-Wav_c@V

      V$V1[,1]<-V$V1[,1]-V$V1[,1]
      V$V2[,1]<-V$V2[,1]-V$V2[,1]
      V$V3[,1]<-V$V3[,1]-V$V3[,1]
      V$V4[,1]<-V$V4[,1]-V$V4[,1]
      V$V5[,1]<-V$V5[,1]-V$V5[,1]
      V$V6[,1]<-V$V6[,1]-V$V6[,1]
      V$V7[,1]<-V$V7[,1]-V$V7[,1]
      Wav_c@V<-V

      RmnObj$Data[,b+(a-1)*RmnObj$BandsLength]<-imodwt(Wav_c)[(window+1):(window+RmnObj$BandsLength)]
      }

    }
  }

  if (ZroThreshold == TRUE)
  {
    for(i in 1:RmnObj$numBands)
    {
      for(j in 1:RmnObj$numPixels)
      {
        for(k in 1:RmnObj$BandsLength)
        {
          if((RmnObj$Data[k,j+(i-1)*RmnObj$BandsLength])<0)
          {
            RmnObj$Data[k,j+(i-1)*RmnObj$BandsLength] <- 0L
            RmnObj$Data[k,j+(i-1)*RmnObj$BandsLength]
          }
        }
      }
    }
  }

  RmnObj$AvrgSpectr <- rCRMI:::AverageSpectrum(RmnObj$numBands,RmnObj$numPixels,RmnObj)

  end_time <- Sys.time()
  execution_time <- end_time - start_time
  print(execution_time)
  writeLines("Processing complete...")

  return (RmnObj)
}

