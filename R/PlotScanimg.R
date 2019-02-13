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

#' PlotScanImage
#'
#' @param rCRSIObj A rCRSI object.
#' @param wavelength Vector containing pairs of indexes that define the area under one or more peaks. (Length must be numImg*2)
#' @param numImg  Number of images to be plotted.
#'
#'
#' @description Generates and plots the image of the "Are under the curve"(AUC) within the given intervals.
#'
#' @export
#'

PlotScanImage <- function(rCRSIObj, wavelength, numImg = 1, interpolate = T)
{

  if(length(wavelength) != numImg)
  {
    stop("Error: Not enough wavelength for the number of images desired")
  }

  for (i in 1:numImg)
  {
    AUCimg <- matrix(data = rCRSIObj$Data[,which.min(abs(rCRSIObj$RamanShiftAxis-wavelength))],
                     nrow = length(rCRSIObj$RelCoords$X),
                     ncol = length(rCRSIObj$RelCoords$Y))

    #Create the raster
    img <- raster::raster(nrows = nrow(AUCimg),
                          ncols = ncol(AUCimg),
                          xmn = 0,
                          xmx = ncol(AUCimg),
                          ymn = 0,
                          ymx = nrow(AUCimg)
    )

    img <- raster::setValues(img,AUCimg)

    raster::plot(x = img,
                 col = sp::bpy.colors(64),
                 useRaster = T,
                 interpolate = interpolate,colNA = "black"
    )

    graphics::title(main = paste("Wavelength: ",signif(rCRSIObj$RamanShiftAxis[which.min(abs(rCRSIObj$RamanShiftAxis-wavelength))],7)),
                    sub = paste("Pixel size: ", signif(rCRSIObj$PixelSize[1],3)," x ",signif(rCRSIObj$PixelSize[2],3)),
                    line = 2
    )
  }
}
