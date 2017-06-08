# -----------------------------------------------------------------------------------
#' @title		
#'     Smooth and Plot One Dimensional Kolmogorov-Zurbenko Periodogram 
#'
#' @description	
#'    \code{kz.smpg} is designed to smooth and plot 1D KZ periodogram easily.
#'    It will calculate the raw periodogram, mark the spikes, 
#'    smooth the periodogram, and then output the plot.
#'
#' @details
#'	  The smoothing process is based on a modified DiRienzo-Zurbenko (DZ) 
#'	  method, for which the smoothing window is not symmetric around the 
#'	  value point. The smoothing algorithm is implemented in C.
#' 
#' @param      x 	The data vector for analyses. Missing values are allowed.
#' @param   dpct 	A pre-specified percentage of total variation. 
#'			Defaults to 1\%.
#' @param     rg 	The frequency range of the outputted periodogram.
#'			Default is 0 to 0.5. 
#' @param   plot	TRUE or FLASE. Flag for output periodogram plot or not. 
#'                Defaults to FLASE.
#' @param    log	TRUE or FLASE. Use log scale for output periodogram. Defaults to FLASE.
#' @param    ...	Other arguments. 
#' \itemize{
#'  \item	\code{m : } The window size for a regular Fourier transform
#'  \item	\code{k : } The number of iterations for the KZFT
#'  \item	\code{n : }	The sampling frequency rate as a multiplication 
#'				of the Fourier frequencies
#'  \item	\code{p : } The distance between two successive intervals as 
#'				a percentage of the total length of the data series
#'  \item   \code{w : }	Size of smoothing window. Default value is 20.
#'  \item \code{lvl : }														
#'		"min" or "max". Threshold strategy for marking frequency 
#'		spikes. "min" is used for cases with weak singles mixed 
#'		with dominating strong spikes. Defaults to "max".
#'  \item \code{cut : }							
#'		Set the minimum value for a marked frequency spike. Recommend
#'		to use argument \code{lvl} instead of setting this value directly.
#' }
#'
#' @return 		Data frame for outputted periodogram, including column \emph{spg} 
#'				for the periodogram values, and \emph{freq} for the frequencies.
#' @rdname 		smpg
#' @name 		smpg
#' @useDynLib 	kzfs, .registration = TRUE
#' @export		kz.smpg
#' @keywords 	KZ-periodogram
#' @concept     	Kolmogorov-Zurbenko periodogram
#' @export
#' @seealso		\code{\link[kzft]{kzp}}, \code{\link{kzp2}}, \code{\link{kz.ft}}
#' @examples 
#' ## Adapted from kzft::kzp example 2
#'	t <- 1:2000
#'	y <- 1.1*sin(2*pi*0.0339*t)+7*sin(2*pi*0.0366*t)+5*rnorm(length(t),0,1)
#'	y[sample(t,100,replace=FALSE)] <- NA
#'	
#'	\dontrun{
#'	# system.time(op <- kz.smpg(y, dpct=0.0001, rg=c(0.025,0.05),  
#'	#		plot=TRUE, log=TRUE, lvl="min", n=10, k=2))
#'	}
#'	op <- kz.smpg(y, dpct=0.0000, f=c(0.0339,0.0366), rg=c(0.025,0.05), 
#'		n=10, k=2, plot=TRUE, lvl="min", log=FALSE)
# -------------------------------------------------------------------------------------

kz.smpg <- function(x, dpct=0.01, rg=c(0,0.5), log=F, plot=F, ...) {
    if (!is.null(dim(x))) { x <- as.vector(x) }
	dots <- list(...)
	if (hasArg("k")) { k <- dots$k } else { k <- 1 } 
	if (hasArg("n")) { n <- dots$n } else { n <- 1 } 
	if (hasArg("w")) { w <- dots$w } else { w <- 20*n }
	if (hasArg("m")) { 
		m <- dots$m; 
	} else { 
		m <- floor(length(x)/(2*k))*2 
	}
	kzp.x <- kz.ft(x=x, m=m, ...)$pg
	ln <- length(kzp.x)
	omega <- seq(0,0.5,length=(ln+1))[2:(ln+1)]
	idx <- omega<=rg[2] & omega>=rg[1]	
	omega <- omega[idx]
	lmt <- min(length(omega),w)
	if (log) {
	   kzp.xp <- log((kzp.x[idx]) + ifelse(min(kzp.x[idx])==0, 0.01, 0))
	} else {
	   kzp.xp <- kzp.x[idx]
	}
	if (dpct==0) { 
	   spg.xp <- kzp.xp
	} else {
	   spg.xp <- smooth.kzp(kzp.xp, dpct=dpct, w=lmt)
	}
	if (plot) {
	   spg.xp[is.na(spg.xp)] <- 0
	   tmp <- smpg.plot(freq=omega,spg=spg.xp,Title="",dpct=dpct,...)
	}
	return(data.frame(freq=omega, spg=spg.xp, rpg=kzp.x[idx]))	
}


# -----------------------------------------------------------------------
#   1D DZ Algorithm (Improved function kzft::smooth.kzp())
#
#  @param    rpg 	Vector of raw periodogram.
#  @param 	   w  	Size of smoothing window. 
#  @rdname  smpg
#  @export
# -----------------------------------------------------------------------

smooth.kzp <- function(rpg, dpct, w = length(rpg))
{
   storage.mode(rpg) <- "double"
   sp <- .Call("kzpg", as.vector(rpg), as.numeric(dpct), as.integer(w))
   sp
}


# ------------------------------------------------------------------
#         Function to plot the periodogram plot
#
#  @rdname 	smpg
#  @param    spg	Data vector of periodogram values.
#  @param   freq	Data vector of frequency values.
#  @param  Title	String. Used to mark the periodogram plot.
#  @param  angle	Direction of periodogram. Default is missing.
#  @export
# ------------------------------------------------------------------

smpg.plot <- function(spg, freq, Title, angle, ...) {
	dots <- list(...)
	method <- "DZ"
	if (hasArg("dpct")){ dpct <- dots$dpct} else { dpct <- 0.01 }
	if (hasArg("raw")) { raw <- dots$raw } else { raw <- FALSE }
	if (raw) { dpct <- 0 }
	if (missing(Title)) Title <- "Smoothed Periodogram"
	xl <- "Frequency (cycles/unit interval)"
	tl <- ""
	if (exists("dpct")) tl <- paste(100*dpct,"%", sep="")
	tl <- paste(tl, method)
	tl <- paste(tl, length(spg),"points")
	if (!(missing(angle)))  {
	    rrg <- round(180*angle/pi,2)
	    tl <- paste(tl," : ",rrg,enc2utf8("\xB0"), sep="")
	    okag <- paste(round((180/pi)*angle,2),enc2utf8("\xB0"),sep="")
	} else {
	    okag <- ""
	}
	plot(x=freq, y=spg, main=Title, type="l", xlab=xl, ylab=" ")
	mtext(tl, col="grey40", cex=0.85, line=0.5, font=1)
	Pds <- markspikes(x.fq=c(0,freq),y.spm=c(min(spg),spg), plot=TRUE, ...)
	if (length(Pds)==0) Pds <- NA
	cat(okag," Periodogram marked at frequency",Pds,"\n")
	return(data.frame(direction=okag, freq=Pds))
}


