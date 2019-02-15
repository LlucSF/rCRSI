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

#' PlotClusterImage
#'
#' @param path full path to .txt file.
#' @param numBands Number of .txt files to process.
#'
#' @description Reads a .txt file containing Raman data from Reinshaw devices and convert it to an RamanR data object.
#'
#' @return a rCRSIObj data object.
#' @export
#'

PlotClusterImage <- function(rCRSIObj, clusterIndex, interpolate = T)
{

  Clusimg <- matrix(data = clusterIndex,
                    nrow = length(rCRSIObj$RelCoords$X),
                    ncol = length(rCRSIObj$RelCoords$Y)
                    )


  #Create the raster
  img <- raster::raster(nrows = nrow(Clusimg),
                        ncols = ncol(Clusimg),
                        xmn = 0,
                        xmx = ncol(Clusimg),
                        ymn = 0,
                        ymx = nrow(Clusimg)
                        )

  img <- raster::setValues(img,Clusimg)

  col.palette <- rainbow(length(levels(factor(clusterIndex))))

  raster::plot(x = img,
               col = col.palette,
               useRaster = T,
               interpolate = interpolate,
               legend = FALSE,
               axes = FALSE,
               box = FALSE
               )

  graphics::title(main = "Clustering image",
                  sub = paste("Pixel size: ", signif(rCRSIObj$PixelSize[1],3)," x ",signif(rCRSIObj$PixelSize[2],3)),
                  line = 2
                  )

  graphics::legend("right",
                   legend = levels(factor(clusterIndex)),
                   col = col.palette,
                   pch = 15
                   )
}

