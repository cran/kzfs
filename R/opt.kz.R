# -----------------------------------------------------------------------------------
#' @title	
#'	   Improve Accuracy of KZ Periodogram Estimation with Optimization
#'
#' @description	
#'    Functions in this group are designed to improve the estimated wave parameters 
#' based on optimization of KZ directional periodograms and 2D periodograms. 
#'
#' @param	a 		Data array. Wave signals plus noise. 
#' @param	rec		Data list. For \code{optDR}, it is the outputs of function kzpdr. 
#' 				It includes the marked spectrum spike frequency values, directions.
#' @param	delta		Searching range parameter for optimization. Default is 0.005.
#' @param 	...		Other arguments. 
#' \itemize{
#'  \item	\code{k : } Integer. The iterations number of KZFT.
#'  \item	\code{n : }	The sampling frequency rate as a multiplication 
#'				of the Fourier frequencies
#' }	
#' 
#' @rdname 		opt.kz
#' 
#' @details
#' \code{optDR} optimizes estimations of directional periodograms using R
#' function \code{stats::optimize}. \code{optD2P} works for estimations of 2D 
#' periodograms; related optimization process is based on function \code{stats::optim}.
#'
#' @return
#' 	\code{optDR} will return the data frame of detailed estimation 
#'	for each direction. 
# 
#' @export
#' @seealso		\code{\link{kzpdr}}, \code{\link{kzp2}}
# 
#' @examples
# 
#'	dx <- 300				# x range
#'	dy <- 300				# y range
#'	b <- expand.grid(x=1:dx, y=1:dy)
#'	q1 <- pi/3; f1 <- 0.2;
#'	b$v1 <- sin(f1*2*pi*(b$x*cos(q1)+b$y*sin(q1))+100*runif(1))
#'	q2 <- pi/6; f2 <- 0.05;
#'	b$v2 <- sin(f2*2*pi*(b$x*cos(q2)+b$y*sin(q2))+100*runif(1))
#'	a <- array(0,c(dx,dy))
#'	a[as.matrix(b[,1:2])] <- b$v1 + 1.5*b$v2
#'	noise <- 5*matrix(rnorm(dx*dy,0,1),ncol=dy)
#'
#'	# Identifying with 2D periodogram
#'	# kzp2.demo <- kzp2(a+noise)$kzp2d
#'	QF <- kzp2.summary(kzp2.demo, num=2)
#'
#'	# Optimization of the 2D periodogram
#'	# It may take 1 to 5 minutes
#'	# kzp2.QF <- optD2P(a+noise, QF, k1=1, n1=1)
#'	kzp2.QF
#'
#'	# Optimization of directional periodogram
#'	# It may take 10 to 20 minutes
#'	# kzpdr.demo <- kzpdr(a+noise, c(0,45,15,35)*pi/180, plot=TRUE, dpct=0.05)$rec
#'	# kzpde.QF <- optDR(a+noise, kzpdr.demo)
#'
#'	QF2 <- kzpdr.eval(kzpdr.demo)
#'	opt.QF2 <- kzpdr.eval(kzpdr.QF)
#'
# -----------------------------------------------------------------------------------

optDR <- function(a, rec, delta=0.005, ...){
   ok.rec <- rec
   dots <- list(...)
   if (hasArg("k")) { k0 <- dots$k } else { k0 = 1}
   if (hasArg("n")) { n0 <- dots$n } else { n0 = 1}
   if (hasArg("off")) { off <- dots$off } else { off = 0}
   DR <- function(x) { 
	kzpdr(a, Ag, f=x, frun=T, plot=F, pair=F, k=k0, n=n0, off=off)$rec[1,"spg"]
   }
   for (i in 1:length(rec$freq)) {
	Ag   <- rec[i,]$ang
	itv  <- c(rec[i,]$freq - delta, rec[i,]$freq + delta)
	okDR <- stats::optimize(f=DR, maximum = TRUE, interval=itv)
	ok.rec[i,]$freq <- okDR[[1]]
	ok.rec[i,]$spg  <- okDR[[2]]
   }
   ok.rec
}


# ----------------------------------------------------------------------
#     Function to mark the spikes of 1D periodogram
#
#' @rdname 	opt.kz
#' @export	
# ----------------------------------------------------------------------

optD2P <- function(a, rec, ...){
   ok.rec <- rec
   dots <- list(...)
   if (hasArg("k")) { k1 <- dots$k } else { k1 = 1}
   if (hasArg("n")) { n1 <- dots$n } else { n1 = 1}
   fitD2P <- function(a, x0, y0, k0, n0) { 
	D2P <- function(xy) {
	   -sum(kzp2(a, rec=list(x=xy[1], y=xy[2]), k=k0, n=n0)$kzp2d)
	}
	stats::optim(c(x0,y0), D2P)
   }
   for (i in 1:length(rec$x)) {
	okD2P <- fitD2P(a, x0=rec$x[i], y0=rec$y[i], k0=k1, n0=n1)
	ok.rec$x[i] <- okD2P$par[1]
	ok.rec$y[i] <- okD2P$par[2]
	ok.rec$z[i] <- -1*okD2P$value
	ok.rec$frequency[i] <- sqrt(okD2P$par[1]^2 + okD2P$par[2]^2)
	ok.rec$direction[i] <- atan(okD2P$par[2]/okD2P$par[1])*180/pi
   }
   ok.rec
}

