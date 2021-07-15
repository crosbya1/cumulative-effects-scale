
substr(gridnum$StationKey[1], 10, 11)

as.numeric(substr(gridnum$StationKey[1], 10, 11))

statnum <- data.frame(ss = gridnum$StationKey, grid = gridnum$Gridnum, point = sapply(1:nrow(gridnum), function(x) as.numeric(substr(gridnum$StationKey[x], 9, 11))))


# Function to create a moving window over a sequentially numbered grid of points, where each row is a point and the row number reflects it's position in a matrix numbered vertically from top to bottom and then left to right. 

# Returns a list, where each list item rperesents a window, and consists of the row numbers for points within that window. 

# x.dim is the x-dimension of the matrix (number of columns)

# y.dim is the y-dimension of the matrix (number of rows)

# win.size is the dimension of the square window (number of rows x number of columns): ex. c(3, 3) for a square window or c(3, 4) for a regtangular window

win.list <- function(x.dim, y.dim, win.size){
  t.y <- y.dim - win.size[2] + 1 # The number of windows fit going down the y dimension
  t.x <- x.dim - win.size[1] + 1 # The number of windows going across the x dimension
  ts <- 1:t.y # Numbered sequence of windows in the y dimension
  tx <- 1:t.x
  
  p <- do.call(c, lapply(1:t.x, function(i){
    s.y <- ((1:win.size[1]) - 1)*y.dim + (tx[i] - 1)*y.dim
    d <- lapply(1:t.y, function(j){
      v <- ts[j]
      u <- unlist(lapply(1:win.size[1], function(x) (v:(v + win.size[2] - 1)) + s.y[x]))
      return(u)
    })
    return(d)
  }))
  return(p)
}




