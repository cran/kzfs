# -----------------------------------------------------------------------------------
#' @title		
#'     Smooth and Plot One Dimensional Kolmogorov-Zurbenko Periodogram 
#'
#' @description	
#'    Functions designed to smooth and plot 1D KZ periodogram easily.
#'
#'    Function \code{smooth.kzp} is an improved version of DZ smoothing
#'    algorithm with asymmetrical window sizes, implemented in C.
#'
#'    \code{plot.smpg} outputs 1D KZ periodogram plots, marks the 
#'    frequencies for spectral spikes on smoothed periodogram.
#'
#'    \code{kz.smpg} is a user-friendly wrapper of 3 functions:
#'    \code{kzft::kzp}, \code{smooth.kzp}, and \code{plot.smpg}.
#'    It will calculate the raw periodogram, mark the spikes, 
#'    smooth the periodogram, and then output the plot.
#' 
#' @param     sd 	The data vector for analyses. Missing values are allowed.
#' @param   dpct 	A pre-specified percentage of total variation. 
#'			Defaults to 1\%.
#' @param     rg 	The frequency range of the outputted periodogram.
#'			Default is 0 to 0.5, or 2 to infinite time-steps per period
#' @param   plot	TRUE or FLASE. Flag for output periodogram plot or not. 
#'                Defaults to FLASE.
#' @param    log	TRUE or FLASE. Use log scale for output periodogram. Defaults to FLASE.
#' @rdname 		smpg
#' @name 		smpg
#' @useDynLib 	kzfs
#' @export
#' @seealso		\code{\link[kzft]{kzp}}
#' @examples 
#' ## Adapted from kzft::kzp example 2
#'	t <- 1:2000
#'	y <- 1.1*sin(2*pi*0.0355*t)+7*sin(2*pi*0.0365*t)+5*rnorm(length(t),0,1)
#'	op <- kz.smpg(y, dpct=0.0010, rg=c(0.025,0.05), plot=TRUE, log=TRUE)
#'	op <- kz.smpg(y, dpct=0.0001, rg=c(0.025,0.05), plot=TRUE, lvl="min")

# -------------------------------------------------------------------------------------

kz.smpg <- function(sd, dpct=0.01, rg=c(0,0.5), log=F, plot=F, ...) {
    if (!is.null(dim(sd))) { sd <- as.vector(sd) }
	dots <- list(...)
	if (hasArg("k")) { k <- dots$k } else { k <- 1 } 
	if (hasArg("p")) { p <- dots$p } else { p <- 1 }
	if (hasArg("n")) { n <- dots$n } else { n <- 1 }
	if (hasArg("w")) { w <- dots$w } else { w <- 20 }
	if (hasArg("m")) { m <- dots$m } else { 
		m <- floor(length(sd)/(2*k))*2 
	}
	sd[is.na(sd)] <- 0
	kzp.x <- kzft::kzp(x=sd, m, k, p, n)
	ln <- length(kzp.x)
	omega <- seq(0,0.5,length=(ln+1))[2:(ln+1)]
	idx <- omega<=rg[2] & omega>=rg[1]	
	omega <- omega[idx]
	lmt <- min(length(omega),w)
	if (log) {
		kzp.xp <- log(kzp.x[idx])
	} else {
		kzp.xp <- kzp.x[idx]
	}
	spg.xp <- smooth.kzp(kzp.xp, dpct=dpct, w=lmt)
	if (plot) {
	   spg.xp[is.na(spg.xp)] <- 0
	   tmp <- smpg.plot(freq=omega,spg=spg.xp,Title="",dpct=dpct,...)
	}
	return(data.frame(freq=omega, spg=spg.xp, rpg=kzp.x[idx]))	
}


# -----------------------------------------------------------------------
#   1D DZ Algorithm (Improved function kzft::smooth.kzp())
#
#' @param    rpg 	Vector of raw periodogram.
#' @param 	   w  	Size of smoothing window. 
#' @rdname  smpg
#' @export
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
#' @rdname 	smpg
#' @param    spg	Data vector of periodogram values.
#' @param   freq	Data vector of frequency values.
#' @param  Title	String. Used to mark the periodogram plot.
#' @param  angle	Direction of periodogram. Default is missing.
#' @param    ...	Other arguments. See details section for more information.
#' @details
#' Other arguments for function \code{kz.smpg}:
#'
#' \itemize{
#'  \item	\code{m : } The window size for a regular Fourier transform
#'  \item	\code{k : } The number of iterations for the KZFT
#'  \item	\code{n : }	The sampling frequency rate as a multiplication 
#'				of the Fourier frequencies
#'  \item	\code{p : } The distance between two successive intervals as 
#'				a percentage of the total length of the data series
#'  \item   \code{w : }	Size of smoothing window. Default value is 20.
#' }
#'
#' Other arguments for function \code{smpg.plot}:
#'
#' \itemize{
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
#'			for the periodogram values, and the frequency \emph{freq}.
#' @keywords 	KZ-periodogram
#' @concept     	Kolmogorov-Zurbenko periodogram
#' @export
# ------------------------------------------------------------------

smpg.plot <- function(spg, freq, Title, dpct, angle, ...) {
	dots <- list(...)
	method <- "DZ"
	if (missing(Title)) Title <- "Smoothed Periodogram"
	xl <- "Frequency (cycles/unit time)"
	tl <- ""
	if (!(missing(dpct)))   tl <- paste(100*dpct,"%", sep="")
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
	Pds <- markspikes(x.fq=c(0,freq),y.spm=c(min(spg),spg), ...)
	if (length(Pds)==0) Pds <- NA
	cat(okag," Periodogram marked at frequency",Pds,"\n")
	return(data.frame(direction=okag, freq=Pds))
}


