#########################################################################
#     rCRSI - R package for Confocal Raman Spectroscopy Imaging rCRSIObj processing and visualization
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
#' @param rCRSIObj a rCRSIObj.
#' @param OddfiltLength Integer. Length of the Savitzky-Golay filter.
#' @param NtoZ  Boolean. Apply negative-to-zero thresholding during the wavelets preprocessing.
#' @param BackgroundSubs Boolean. Apply background substraction using a low frequency filter wavelets processing.
#' @param window Integer. Number of added points at the begining and ending of the signal to reduce artifacts from the wavelets transform.
#'
#' @description  Smooths the data contained in a rCRSIObj using the discrete wavelets transform (DWT).
#'
#' @return a rCRSIObj object.
#' @export
#'

SmoothSpectra <- function(rCRSIObj,method = "normal", OddfiltLength = 7,
                          NtoZ = TRUE,
                          BackgroundSubs = TRUE, window = 32
                          )
{
  if(method!="normal" & method!="wavelets")
  {
    stop("Method must be 'normal' or 'wavelets'")
  }

  writeLines("Starting processing...")
  start_time <- Sys.time()
  OldrCRSIObj <- rCRSIObj

  if(method=="wavelets")
  {
    for(a in 1:rCRSIObj$numBands)
    {
      for(b in 1:rCRSIObj$numPixels)
      {
        spectrum <- c(rep(x = rCRSIObj$Data[(b+((a-1)*rCRSIObj$numPixels)),1],times = window),
                      rCRSIObj$Data[(b+((a-1)*rCRSIObj$numPixels)),],
                      rep(x = rCRSIObj$Data[(b+((a-1)*rCRSIObj$numPixels)),rCRSIObj$BandsLength],times = window))

        Wav_c <- wavelets::modwt(X = spectrum, n.levels = 4, filter = "d4",
                                 boundary = "periodic", fast = TRUE)
        W <- Wav_c@W
        W$W1 <- W$W1-W$W1
        W$W2 <- W$W2-W$W2
        W$W3 <- W$W3-W$W3
        # W$W4 <- W$W4-W$W4
        Wav_c@W <- W

        rCRSIObj$Data[b+(a-1)*rCRSIObj$numPixels,]<-wavelets::imodwt(Wav_c)[(window+1):(window+rCRSIObj$BandsLength)]

        if (BackgroundSubs)
        {
          spectrum <- c(rep(x = rCRSIObj$Data[(b+((a-1)*rCRSIObj$numPixels)),1],times = window),
                        rCRSIObj$Data[(b+((a-1)*rCRSIObj$numPixels)),],
                        rep(x = rCRSIObj$Data[(b+((a-1)*rCRSIObj$numPixels)),rCRSIObj$BandsLength],times = window))

          Wav_c <- wavelets::modwt(X = spectrum, n.levels=9, filter="d6",
                                   boundary="periodic", fast=TRUE)
          W <- Wav_c@W
          # W$W1 <- W$W1-W$W1
          # W$W2 <- W$W2-W$W2
          # W$W3 <- W$W3-W$W3
          W$W4 <- W$W4-W$W4
          W$W5 <- W$W5-W$W5
          W$W6 <- W$W6-W$W6
          W$W7 <- W$W7-W$W7
          W$W8 <- W$W8-W$W8
          W$W9 <- W$W9-W$W9
          Wav_c@W <- W

          V <- Wav_c@V
          V$V1[,1] <- V$V1[,1]-V$V1[,1]
          V$V2[,1] <- V$V2[,1]-V$V2[,1]
          V$V3[,1] <- V$V3[,1]-V$V3[,1]
          V$V4[,1] <- V$V4[,1]-V$V4[,1]
          V$V5[,1] <- V$V5[,1]-V$V5[,1]
          V$V6[,1] <- V$V6[,1]-V$V6[,1]
          V$V7[,1] <- V$V7[,1]-V$V7[,1]
          V$V8[,1] <- V$V8[,1]-V$V8[,1]
          V$V9[,1] <- V$V9[,1]-V$V9[,1]
          Wav_c@V <- V

          rCRSIObj$Data[b+(a-1)*rCRSIObj$numPixels,] <- wavelets::imodwt(Wav_c)[(window+1):(window+rCRSIObj$BandsLength)]
        }
      }
    }

    if (NtoZ)
    {
      for(i in 1:rCRSIObj$numBands)
      {
        for(j in 1:rCRSIObj$numPixels)
        {
          for(k in 1:rCRSIObj$BandsLength)
          {
            if((rCRSIObj$Data[j+(i-1)*rCRSIObj$numPixels,k])<0)
            {
              rCRSIObj$Data[j+(i-1)*rCRSIObj$numPixels,k] <- 0L
            }
          }
        }
      }
    }
  }


  if(method=="normal")
  {
    for(i in 1:rCRSIObj$numSpectr)
    {
      rCRSIObj$Data[i,] <- signal::sgolayfilt(x = rCRSIObj$Data[i,],p = 2,n = OddfiltLength)
    }
  }

  rCRSIObj$ProAvrgSpectr <- AverageSpectrum(rCRSIObj$numBands,rCRSIObj$numPixels,rCRSIObj,F)
  Old_average <- AverageSpectrum(OldrCRSIObj$numBands, OldrCRSIObj$numPixels, OldrCRSIObj, F)

  for (i in 1:rCRSIObj$numBands)
  {
    df <- data.frame(y = c(Old_average[,i]/max(Old_average[,i]),
                           rCRSIObj$ProAvrgSpectr[,i]/max(rCRSIObj$ProAvrgSpectr[,i])),
                     x = rep(rCRSIObj$RamanShiftAxis[,i], times = 2),
                     cl = rep(c("Raw","Processed"), each = length(rCRSIObj$ProAvrgSpectr[,i])))

    g <- ggplot2::ggplot() + ggplot2::theme_bw() +
      ggplot2::geom_line(mapping = ggplot2::aes(x = df$x, y = df$y, colour = df$cl)) +
      ggplot2::labs(x ="Raman Shift", y = "Counts") +
      ggplot2::ggtitle("Data processing normalized results") +
      ggplot2::labs(colour = "Data") +
      ggplot2::scale_x_continuous(breaks = trunc(seq(from = min(df$x), to = max(df$x),length.out = 30)),minor_breaks = ggplot2::waiver())
    print(g)
  }

  end_time <- Sys.time()
  execution_time <- end_time - start_time
  print(execution_time)
  writeLines("Processing complete...")

  return (rCRSIObj)
}

