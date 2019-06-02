# -----------------------------------------------------------------------------------
#' @title		
#' 		Check Images' Motion Scales with 2D KZ Periodogram Signals
#'
#' @description
#'    Functions used to reveal directional and scale information
#' with 2D KZ periodograms for spatial motions covered by heavy noises. 
#'
#'    One can get 2D raw periodogram with function \code{kzp2}, and smooth the
#' 2D periodogram with function \code{smooth.kzp2}. 
#' 
#'    Function \code{summary.kzp2} can help to summarize direction and frequency 
#'  information from smoothed 2D KZ periodogram. The input should be a 2D KZ
#'  periodogram data with frequency range (0, 0.5] on both x- and y- axis.
#'
#' @rdname  kzp2
#' @param     x		Data array of 2D wave field. Missing values are allowed.
#'					Limited to 2D arrays for current version.
#' @param     m		The window size for a regular Fourier transform.
#'					Default value is set to data array size. 
#' @param     k		The number of iterations for the KZFT. Default is 1.
#' @param   rpg 	Array of raw 2D periodogram. Usually it is part of output of \code{kzp2}.
#' @param   spg 	Array of smoothed 2D periodogram. It could be output of \code{summary.kzp2}.
#' @param   ... 	Arguments to be passed to methods.
#' \itemize{
#'	 \item	\code{k : } The number of iteration times of KZFT
#'	 \item	\code{n : } The sampling frequency rate as a multiplication 
#'					of the Fourier frequencies
#'	 \item	\code{p : } The distance between two successive intervals as 
#'					a percentage of the total length of the data series
#' }
#'
#' @details		KZ 2D raw spectrum is calculated based on \code{kz.ft}.
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
#'			parameters, including frequency and direction values.
#'
#' @keywords 	2D periodogram 
#' @concept 	2D periodogram
#' @export
#' @seealso		\code{\link{kzpdr}}, \code{\link{kzpdr.eval}}, \code{\link{kzpdr.spikes}}
#
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

kzp2 <- function(x, k = 1, m = dim(x)/k, ...) {
   x <- x - min(x, na.rm = TRUE) + 1
   x[is.na(x)] <- 0
   dots <- list(...)
   if (hasArg("n")) { n <- dots$n } else { n <- 1 }    
   if (hasArg("rec")) { 
	dx <- length(dots$rec$x)
	dy <- length(dots$rec$x)
	 n <- 1
   } else {
	tmp <- kz.ft(x[,1],m=m[1],k=k, adpt=F, phase=F, ...)
	dx <- length(tmp$fft)
	tmp2 <- kz.ft(x[1,], m=m[2],k=k, adpt=F, phase=F, ...)
	dy <- length(tmp2$fft)
   }
   alpha <- array(0, c(dx, dim(x)[2]))
   beta  <- array(0, c(dx, dy))
   for (i in (1:dim(x)[2])) {
	if (hasArg("rec")) { 
	   tmp <- kz.ft(x[,i],m=m[1],k=k, adpt=F, phase=F, f=dots$rec$x, ...)
	} else {
   	   tmp <- kz.ft(x[,i],m=m[1],k=k, adpt=F, phase=F, ...)
	}
	alpha[,i] <- tmp$fft
   }
   for (i in (1:dx)) {
	if (hasArg("rec")) { 
	   tmp2 <- kz.ft(alpha[i,],m=m[2],k=k, adpt=F, phase=F, f=dots$rec$y, ...)
	} else {
	   tmp2 <- kz.ft(alpha[i,],m=m[2],k=k, adpt=F, phase=F, ...)
	}
	beta[i,] <- tmp2$fft
   }
   zx <- tmp$f
   zy <- tmp2$f
   if (hasArg("rec")) { 
	bx <- which(tmp$pg>0)
	by <- which(tmp2$pg>0)
	z0 <- (abs(beta))^2
	z <- array(0, c(length(zx), length(zy)))
	for (i in 1:dim(z0)[1]) {z[bx[i],by[i]] <- z0[i,i]}
   } else {
	z <- (abs(beta))^2
	z <- z[1:length(zx),1:length(zy)]
   }
   return(list(freq.x=zx,freq.y=zy,kzp2d=z,kzft2=beta))
}


# ------------------------------------------------------------------------------------
#      Smoothing Method for 2D Periodogram (2D DZ Algorithm)
#' @rdname  kzp2
#' @param	   w	Smoothing window size. 
#' @param   dpct	A pre-specified percentage of total variation.
#'			Default value is 1\%.
#' @export
# ------------------------------------------------------------------------------------
# version: Dec. 18, 2016

smooth.kzp2 <- function(rpg, dpct = 0.01, w = round(dim(rpg)/4), k=1, ...)
{
   dots <- list(...)
   if (hasArg("n"))     { n <- dots$n } else { n <- 1 }    
   if (hasArg("cut"))   { cut <- dots$cut } else { cut <- 0.1 }    
   if (hasArg("delta")) { delta <- dots$delta } else { delta <- 1 }    
    N1 <- dim(rpg)[1]
    N2 <- dim(rpg)[2]
    NAll <- which(is.na(rpg))
    spg <- rpg
    spg[is.na(spg)] <- 0
    if (n > 1) {
	vxy <- c(n,n)
	for (i in 1:vxy[1]) {
	   spg[,i] -> rule
	   vx <- which(rule > (cut))
	   if (length(vx)>0) spg[vx, i] <- NA
	}
	for (j in 1:vxy[2]) {
	   spg[j,] -> rule
	   vy <- which(rule > (cut))
	   if (length(vy)>0) spg[j, vy] <- NA
	}
     while (n > 1 & max(spg, na.rm=T)>1 & all(vxy<=2*n)) {
	vxy <- vxy + 1
	for (i in 1:vxy[1]) {
	   spg[,i] -> rule
	   vx <- which(rule > (cut))
	   if (length(vx)>0) spg[vx, i] <- NA
	}
	for (j in 1:vxy[2]) {
	   spg[j,] -> rule
	   vy <- which(rule > (cut))
	   if (length(vy)>0) spg[j, vy] <- NA
	}
     }
    }
    Nsll <- which(is.na(spg))
    rpg <- spg
    spg[is.na(spg)] <- 0
    for (h in 1:k) {
	w <- w*h
	S1 <- winsize.kzp2(spg, dpct, w[1])
	S2 <- winsize.kzp2(t(spg), dpct, w[2])
	S2 <- t(S2)
	S1 <- S1 - min(S1) + 1
	S2 <- S2 - min(S2) + 1
	ek <- 2*h 
    	for (j in (1:N2)) {
	  for (i in (1:N1)) { 
	    if (min(S2[i, j], S1[i, j]) > (delta+1)) {
		mm1 <- min(S1[max(1,i-ek):min(i+ek,N1),max(1,j-ek):min(j+ek,N2)])
		mm2 <- min(S2[max(1,i-ek):min(i+ek,N1),max(1,j-ek):min(j+ek,N2)])
		# if switch mm1 and mm2, it will be different.
		ly <- max( 1, (j - mm2 + 1)) 
		ry <- min(N2, (j + mm2 - 1))
		lx <- max( 1, (i - mm1 + 1)) 
		rx <- min(N1, (i + mm1 - 1))
		spg[i, j] <- mean(rpg[lx:rx, ly:ry], na.rm=T)
	    } else {
		spg[i, j] <- rpg[i, j]
	    }
	  }
    	}
	rpg <- spg
    }
    spg[NAll] <- NA
    spg[Nsll] <- NA
    return(spg)
}

# -----------------------------------------------------------------------
#      Get Adaptive Smoothing Windows Size for 2D DZ Algorithm
#
#  Window size for simple linear 2D smoother (1 aspect of 2D)
#  Internal function \code{winsize.kzp2} is implemented in C. 
#  It is called twice by \code{smooth.kzp2} to calculate the 
#  adaptive smoothing window sizes of 2D DiRienzo-Zurbenko 
#  smoothing algorithm.
# -----------------------------------------------------------------------
#' @useDynLib 	kzfs, .registration = TRUE

winsize.kzp2 <- function(rpg, dpct, w)
{
   storage.mode(rpg) <- "double"
   w <- .Call("kzp2w", as.array(rpg), as.numeric(dpct), as.integer(w))
   w
}


# -----------------------------------------------------------------------
#       Get Direction and Frequency Information in 2D Periodogram
#
#' @param   rg.x	Frequency range for x direction. Defaults to c(0, 0.5). 
#' @param   rg.y 	Frequency range for y direction. 
#'			Defaults to the same value of the range for x direction.
#' @param 	num	Wave numbers. Defaults to 10. 
#' @export
#' @rdname  kzp2
# -----------------------------------------------------------------------

kzp2.summary <- function (spg, rg.x, rg.y, num=10) 
{
    d.x <- dim(spg)[1]; d.y <- dim(spg)[2];
    spg[is.na(spg)] <- 0
    # cat("\n")
    if (missing(rg.x)) {
	f.x <- (0.5*(1:d.x))/d.x
	# cat("Assume frequency range for x is 0 to 0.5\n")
    } else {
	f.x <- seq(from=rg.x[1], to=rg.x[2], length.out = d.x+1)[-1]
	if (length(rg.y)>2) f.x <- rg.x
    }
    if (missing(rg.y)) {
	f.y <- (0.5*(1:d.y))/d.y
	# cat("Assume frequency range for y is 0 to 0.5\n")
    } else {
	f.y <- seq(from=rg.y[1], to=rg.y[2], length.out = d.y+1)[-1]
	if (length(rg.y)>2) f.y <- rg.y
    }
    px <- spikes.2d(x.fq=rep(f.x,length(f.y)),y.spm=spg, nm=num)
    py <- spikes.2d(x.fq=rep(f.y,length(f.x)),y.spm=t(spg), nm=num)
    if (length(px)==0 | length(py)==0) return("No spikes!")
    compwr <- intersect(px$power, py$power)
    px$freq  <- px$freq[which(px$power %in% compwr)]
    py$freq  <- py$freq[which(py$power %in% compwr)]
    px <- px$freq[1:num]; py <- py$freq[1:num]
    angle <- (180/pi)*atan(py/px)
    r <- sqrt(px^2 + py^2)
    pwr <- compwr[1:num]
    cat("\n")
    px <- px[!is.na(px)]
    py <- py[!is.na(py)]
    pwr <- pwr[!is.na(pwr)]
    r   <- r[!is.na(r)]
    angle <- angle[!is.na(angle)]
    rls <- list(x=px,y=py,z=pwr,frequency=r,direction=angle)
    return(rls)
}


