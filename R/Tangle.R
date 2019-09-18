# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Calculate Tangle
#'
#' Tangle is an interative unfolding algorithm which measures complexity of a time series by forcing a 3-dimension embedding of this time series into an ellipse.
#'
#' @param x A vector of equally spaced time series data.
#' @param tau Integer value for time lag associated with Takens' embedding matrix.
#' @param eps A small value for determining if the untangling procedure should terminate.
#' @return The tangle of \code{x}
#' @examples
#' # Generate data from Lorenz attractor
#' library(nonlinearTseries)
#' lorTS <- lorenz(do.plot = FALSE)$x
#'
#' # Downsample for example speed
#' samps <- seq(1, length(lorTS), length.out = 500)
#' lorTS <- lorTS[samps]
#'
#' # Estimate tau
#' lorTau <- timeLag(lorTS)
#'
#' # Calculate Tangle
#' tangle(lorTS, tau = lorTau)


tangle <- function(x, tau = NA, eps = 5e-2){
  # Error catch

  if (!is.vector(x)){
    stop("x must be a single vector")
  }

  if (!is.numeric(x)){
    stop("x must be a numeric vector")
  }

  if (is.na(tau)){
    stop("Please select a value for tau (e.g., using the timeLag() function)")
  }

  if (tau <= 0 | tau != round(tau)){
    stop("tau must be a positive integer")
  }


  # Tangle Calculation
  xMat <- buildTakens(x, 3, tau) # Embed time series
  xMat <- scale(xMat) # Scale time series

  ## Generate "upshift" matrix
  N <- nrow(xMat)
  S <- diag(N)
  S <- S[ ,c(dim(S)[2],(1:(dim(S)[2]-1)))]
  W <- ((diag(N) + S) / 2)# Generate Updating matrix

  # Start Count
  Ku <- 1
  xMatprev <- xMat

  while(dim(convhulln(cbind(xMat[,1], xMat[,2])))[1] !=N &
        dim(convhulln(cbind(xMat[,1], xMat[,3])))[1] !=N &
        dim(convhulln(cbind(xMat[,2], xMat[,3])))[1] !=N){

    xMat <- W %*% xMat
    xMat <- scale(xMat)

    Ku <- Ku + 1

    if (norm(xMat - xMatprev[c(2:N,1),], type = "2") < eps){
      message("Tangle calculated based on average change less than eps")
      break
    }

    xMatprev <- xMat
  }

  return(log(Ku) / length(x))

}