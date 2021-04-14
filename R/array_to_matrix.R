# Make an array a 2d matrix

#' Input is a data array where the first two dimensions are in space(x,y;lat,lon)
#' and the third dimension is time.
#'
#' Outputs a list, with the non-zero pixels being a matrix (columns = pixels,
#' rows=time steps), s
#'
#' @param data Data in a 3d array, where the first two dimensions are the
#' physical dimensions (e.g., lon and lat), and the third one is time.
#'
#' @return returns the distribution of maxT, stcs, and stcs_mvt, each a vector in a list
#' @export perm_dist
#'
array_to_matrix<- function(data){
  nrows<- dim(data)[1]
  ncols<- dim(data)[2]
  time_points<- dim(data)[3]

  out<- vector("list", length = time_points)
  for (i in 1:time_points){
    out[[i]]<- data[,,i]
  }
  names(out)<- as.character(1:time_points)

  sel<- as.numeric(!is.na(as.vector(out[[1]])))

  if(length(out) > 1){ # when imputing image time series it selects only grid cells with < 10% missing data

    for(i in 2:time_points){
      sel <- as.numeric(!is.na(as.vector(out[[i]]))) + sel
      #print(summary(sel))
    }
    sel<- sel > (time_points*.9) # keep only data points with 90% or more data
  }

  wh.sel<- which(as.vector(sel), arr.ind = TRUE)
  V = sum(sel)
  Y = matrix(NA, nrow=time_points, ncol=V)
  for (time_point in 1:time_points) {
    Y[ time_point , ] = as.vector(out[[time_point]])[sel]
  }

  # also elimiate if scores are constant (e.g., all 0)
  constant<- apply(Y, 2, function(x) diff(range(x, na.rm = TRUE)) < .Machine$double.eps ^ 0.5)
  if(sum(constant, na.rm = TRUE)>0){
    wh_constant<- which(constant)
    Y<- Y[,-wh_constant]
    sel[wh.sel][wh_constant]<- FALSE
    wh.sel<- which(as.vector(sel), arr.ind = TRUE)
  }
  return(list(Y=Y, sel = sel, wh.sel = wh.sel, nrow = nrows, ncol = ncols))
}
