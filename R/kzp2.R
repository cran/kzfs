# -----------------------------------------------------------------------------------
#' @title		
#' 		Check Images' Motion Scales with 2D KZ Periodogram Signals
#'
#' @description
#'    Functions used to reveal directional and scale information
#' with 2D KZ periodograms for spatial motions covered by heavy noises. 
#'
#'    One can get 2D raw periodogram with function \code{kzp2}, and smooth the
#' 2D periodogram with function \code{smooth.kzp2}. Function \code{winsize.kzp2},  
#' which is implemented in C, is used to calculate the adaptive smoothing window sizes 
#' of 2-dimensional DiRienzo-Zurbenko smoothing algorithm of 2D KZ periodogram. 
#' 
#'    Function \code{summary.kzp2} can help to summarize direction and frequency 
#'  information from smoothed 2D KZ periodogram. 
#'
#' @rdname  kzp2
#' @param     x		Data array of 2D wave field. Missing values are allowed.
#'					Limited to 2D arrays for current version.
#' @param     m		The window size for a regular Fourier transform.
#'					Default value is set to data array size. 
#' @param     k		The number of iterations for the KZFT. Default is 1.
#' @param   ... 	Arguments to be passed to methods.
#' \itemize{
#'	 \item	\code{w : } Smoothing window size. Defaults to data array size.
#'	 \item	\code{k : } The number of iteration times of KZFT
#'	 \item	\code{n : } The sampling frequency rate as a multiplication 
#'					of the Fourier frequencies
#'	 \item	\code{p : } The distance between two successive intervals as 
#'					a percentage of the total length of the data series
#' }
#' @param   rpg 	Array of raw 2D periodogram.
#' @param   spg 	Array of smoothed 2D periodogram.
#'
#' @details		KZ 2D raw spectrum is calculated based on \code{kzft::kzft}.
#'			The smoothing method is an extension of \code{kzft::smooth.kzp}.
#'			See introduction of DZ method in \code{kzft::smooth.kzp} for more 
#'			information.
#'
#' @return 		Returned value of function \code{kzp2} is a data list of 
#'			periodogram information, including data array \emph{kzp2d} for 2D 
#'			periodogram values, and two frequency vectors, \emph{freq.x} and  
#'			\emph{freq.y} for \emph{x} and \emph{y} direction, respectively.
#'
#'				\code{smooth.kzp2} only outputs the array of smoothed values.
#'			
#'				\code{kzp2.summary} returns a data list for suggested wave 
#'			paramenters, including frequecy and direction values.
#'
#' @keywords 	2D-periodogram 
#' @concept 	Kolmogorov-Zurbenko periodogram
#' @concept 	2-dimensional periodogram
#' @concept 	2D periodogram
#' @export
#' @seealso		\code{\link[kzft]{kzp}}
#'
#' @examples
#'	dx <- 100				# x range
#'	dy <- 120				# y range
#'	b <- expand.grid(x=1:dx, y=1:dy)
#'	q1 <- pi/6; f1 <- 0.2;
#'	b$v1 <- sin(f1*2*pi*(b$x*cos(q1)+b$y*sin(q1))+100*runif(1))
#'	q2 <- pi/4; f2 <- 0.08;
#'	b$v2 <- sin(f2*2*pi*(b$x*cos(q2)+b$y*sin(q2))+100*runif(1))
#'	a <- array(0,c(dx,dy))
#'	a[as.matrix(b[,1:2])] <- b$v1 + 1.5*b$v2
#'	a <- a + 10*matrix(rnorm(dx*dy,0,1),ncol=dy)
#'
#'	rp <- kzp2(a)			# raw 2D spectrum
#'
#'	fy <- rp$freq.y; fx <- rp$freq.x; rp <- rp$kzp2d
#'	
#'	# smoothing 2D spectrum 2 times
#'	sp <- smooth.kzp2(rp,0.01,k=2)	
#'	
#'	par(mfrow=c(2,1), cex=0.5)
#'	persp(x=fx, y=fy, z=rp, expand =0.5,
#'		main = "Raw 2D KZ Periodogram", ltheta=40, shade=0.75,
#'		theta=-30, phi=15, zlab="",xlab="x", ylab="y",
#'		ticktype="detailed", col="lightblue")
#'	
#'	persp(x=fx, y=fy, z=sp, expand =0.5,
#'		main = "Smoothed 2D KZ Periodogram", ltheta=40, shade=0.75,
#'		theta=-30, phi=25, zlab="",xlab="x", ylab="y",
#'		ticktype="detailed", col="lightblue")
#'	par(mfrow=c(1,1), cex=1)
#'	
#'	kzp2.summary(sp)		# direction & frequency
#'	
# ----------------------------------------------------------------------------------

kzp2 <- function(x, m = dim(x), k = 1, ...) {
   dx <- dim(x)[1]
   dy <- dim(x)[2]
   x[is.na(x)] <- 0
   if (!hasArg("n")) { n <- 1 }
   alpha <- array(0,n*c(dx,dy))
   for (i in (1:dy)) {
   	alpha[,i] <- kzft::kzft(x[,i],m=m[1],k=k,...)$fft
   }
   beta <- array(0,n*c(dx,dy))
   for (i in (1:dx)) {
   	beta[i,] <- kzft::kzft(alpha[i,],m=m[2],k=k,...)$fft
   }
   fx <- (1:(n*dx/2))/(n*dx)
   fy <- (1:(n*dy/2))/(n*dy)
   z <- (abs(beta))^2
   z <- z[1:length(fx),1:length(fy)]
   return(list(freq.x=fx,freq.y=fy,kzp2d=z))
}


# ------------------------------------------------------------------------------------
#            Smoothing Method (DZ Algorithm) for 2D Periodogram 
#' @rdname  kzp2
#' @param	   w	Smoothing window size. Defaults to data array dimension.
#' @param   dpct	A pre-specified percentage of total variation.
#'			Default value is 1\%.
#' @export
# ------------------------------------------------------------------------------------
#
#  The new simple linear 2D smoother (Dec. 23, 2015)
# 
smooth.kzp2 <- function(rpg, dpct = 0.01, w = dim(rpg), k = 1) 
{
    N1 <- dim(rpg)[1]
    N2 <- dim(rpg)[2]
    spg <- array(0, dim = c(N1, N2))
    for (h in 1:k) {
	w <- winsize.kzp2(rpg, dpct, w[1])
    	for (j in (1:N2)) {
	  for (i in (1:N1)) { 
	    dl <- diff(w[max(i-w[1], 1):i,j])
	    dr <- diff(w[i:min(i+w[1],N1),j])
	    ml <- max(w[i,j], min(which(rev(dl)>=0), length(dl)+1))
	    mr <- max(w[i,j], min(which(dr<=0), length(dr)+1))
	    lx <- max( 1, (i - ml + 1)) 
	    rx <- min(N1, (i + mr - 1))
          spg[i, j] <- mean(rpg[lx:rx, j])
	  }
    	}
    	w <- winsize.kzp2(t(spg), dpct, w[2])
    	w <- t(w)
    	for (i in (1:N1)) {
	  for (j in (1:N2)) { 
	    dl <- diff(w[i, max(j-w[1], 1):j])
	    dr <- diff(w[i, j:min(j+w[1],N2)])
	    ml <- max(w[i,j], min(which(rev(dl)>=0), length(dl)+1))
	    mr <- max(w[i,j], min(which(dr<=0), length(dr)+1))
	    ly <- max( 1, (j - ml + 1)) 
	    ry <- min(N2, (j + mr - 1))
          spg[i, j] <- mean(spg[i, ly:ry])
	  }
    	}
	rpg <- spg
    }
    return(spg)
}


# -----------------------------------------------------------------------
#      Get Adaptive Smoothing Windows Size for 2D DZ Algorithm
#
#  Window size for simple linear smoother (1 aspect of 2D)
#
#' @rdname  kzp2
#' @export
# -----------------------------------------------------------------------

winsize.kzp2 <- function(rpg, dpct, w)
{
   storage.mode(rpg) <- "double"
   w <- .Call("kzp2w", as.array(rpg), as.numeric(dpct), as.integer(w))
   w
}


# -----------------------------------------------------------------------
#       Get Direction and Frequency Information in 2D Periodogram
#' @rdname  kzp2
#' @param   rg.x	Frequency range for x direction. Defaults to c(0, 0.5). 
#' @param   rg.y 	Frequency range for y direction. 
#'			Defaults to the same value of the range for x direction.
#' @export
# -----------------------------------------------------------------------

kzp2.summary <- function (spg, rg.x, rg.y=rg.x) 
{
    d.x <- dim(spg)[1]; d.y <- dim(spg)[2];
    cat("\n")
    if (missing(rg.x)) {
	f.x <- (0.5*(1:d.x))/d.x
	cat("Assume frequency range for x is 0 to 0.5\n")
    } else {
	rg.y=rg.x
	f.x <- seq(from=rg.x[1], to=rg.x[2], length.out = d.x+1)[-1]
    }
    if (missing(rg.y)) {
	f.y <- (0.5*(1:d.y))/d.y
	cat("Assume frequency range for y is 0 to 0.5\n")
    } else {
	f.y <- seq(from=rg.y[1], to=rg.y[2], length.out = d.y+1)[-1]
    }
    px <- spikes.2d(x.fq=rep(f.x,length(f.y)),y.spm=spg)
    py <- spikes.2d(x.fq=rep(f.y,length(f.x)),y.spm=t(spg))
    if (all(py$power == px$power) & all(px$power == py$power)) {
	px <- px$freq; py <- py$freq
    } else { 
	px <- px$freq; py <- py$freq
	cat('Inconsistent results!\n')
    }
    len <- min(length(px),length(py))
    px <- px[1:len]; py <- py[1:len]
    angle <- (180/pi)*atan(py/px)
    r <- sqrt(px^2 + py^2)
    o <-(na.omit(match(f.x, px)))
    pwr <- spg[cbind(which(f.x %in% px)[o],which(f.y %in% py)[o])]
    cat("\n")
    rls <- list(x=px,y=py,z=pwr,frequency=r,direction=angle)
    return(rls)
}


