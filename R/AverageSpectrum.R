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




AverageSpectrum <- function(numBands, numPixels, RmnObj)
{
  for (i in 1:numBands)
  {
    for (j in 1:numPixels)
    {
      RmnObj$AvrgSpectr[,i] = RmnObj$AvrgSpectr[,i] + RmnObj$Data[,(j+(numPixels*(i-1)))]
    }
    RmnObj$AvrgSpectr[,i]  = RmnObj$AvrgSpectr[,i] / numPixels
  }

  return(RmnObj$AvrgSpectr)
}




