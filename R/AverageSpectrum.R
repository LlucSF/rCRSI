#########################################################################
#     rCRSI - R package for Confocal Raman Spectroscopy Imaging data processing and visualization
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

AverageSpectrum <- function(numBands, numPixels, rCRSIObj, plot = T)
{
  temp<- array(dim =(c(rCRSIObj$BandsLength,numBands)))

  for (i in 1:numBands)
  {
    for (j in 1:rCRSIObj$BandsLength)
    {
      temp[j,i] <- 0
    }
  }

  for (i in 1:numBands)
  {
    for (j in 1:numPixels)
    {
      temp[,i] <- temp[,i] + rCRSIObj$Data[(j+(numPixels*(i-1))),]
    }
    temp[,i]  <- temp[,i] / numPixels
  }

  if(plot)
  {
    for (i in 1:numBands)
    {
      g <- ggplot2::ggplot(mapping = ggplot2::aes(y = temp[,i], x = rCRSIObj$RamanShiftAxis[,i])) + ggplot2::geom_line(colour = "red") +
        ggplot2::ylim(c(min(temp[,i]), max(temp[,i]))) + ggplot2::ggtitle(label = "Average Spectrum") + ggplot2::labs(x ="Raman Shift", y = "Counts") +
        ggplot2::theme_bw() + ggplot2::scale_x_continuous(breaks = trunc(seq(from = min(rCRSIObj$RamanShiftAxis[,i]), to = max(rCRSIObj$RamanShiftAxis[,i]),length.out = 30)),minor_breaks = ggplot2::waiver())
      print(g)
    }
  }

  return(temp)
}




