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

#' PlotAUCImage
#'
#' @param rCRSIObj A rCRSI object.
#' @param intervals Vector containing pairs of indexes that define the area under one or more peaks. (Length must be numImg*2)
#' @param numImg  Number of images to be plotted.
#'
#'
#' @description Generates and plots the image of the "Are under the curve"(AUC) within the given intervals.
#'
#' @export
#'

PlotAUCImage <- function(rCRSIObj, intervals, numImg = 1)
{
  for(i in 1:numImg)
  {
    if(intervals[2*i-1] > intervals[2*i])
    {
      stop("Error: In intervals imput the lower value first")
    }
  }

  if(length(intervals) != numImg*2)
  {
    stop("Error: Not enough intervals for the number images desired")
  }

  for (i in 1:numImg)
  {
    from <- which.min(abs(rCRSIObj$RamanShiftAxis-intervals[2*i-1]))
    to   <- which.min(abs(rCRSIObj$RamanShiftAxis-intervals[2*i]))

    AUCimg <- matrix(nrow = length(rCRSIObj$RelCoords$Y), ncol = length(rCRSIObj$RelCoords$X))

    for(X in 1:ncol(AUCimg))
    {
      for(Y in 1:nrow(AUCimg))
      {
        AUCimg[Y,X] <-  flux::auc(x = rCRSIObj$RamanShiftAxis[from:to],
                                  y = rCRSIObj$Data[(Y+(X-1)*length(rCRSIObj$RelCoords$Y)), from:to]
                                  )
      }
    }

    #Create the raster
    img <- raster::raster(nrows = nrow(AUCimg),
                          ncols = ncol(AUCimg),
                          xmn = 0,
                          xmx = ncol(AUCimg),
                          ymn = 0,
                          ymx = nrow(AUCimg)
                          )

    raster::values(img) <- AUCimg

    raster::plot(x = img,
                 col = sp::bpy.colors(64),
                 useRaster = T,
                 interpolate = T,
                 )

    graphics::title(main = paste("AUC: ",signif(rCRSIObj$RamanShiftAxis[from],5),"~",signif(rCRSIObj$RamanShiftAxis[to],5)),
          sub = paste("Pixel size: ", signif(rCRSIObj$PixelSize[1],3)," x ",signif(rCRSIObj$PixelSize[2],3)),
          line = 2
          )
  }
}
