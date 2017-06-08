# -----------------------------------------------------------------------------------
#' @title		
#'     	Kolmogorov-Zurbenko Fourier Transform Function
#'
#' @description	
#'    \code{kz.ft} is improved version of Wei Yang's \code{kzft::kzft}. 
#'    It has been modified to handle missing values in the data. 
#'	Besides KZ Fourier transform, the outputs also include KZ periodogram.
#'
#'    \code{kz.ftc} is an experimental version of KZFT for signals sampled 
#' 	on continual time/space points with irregular intervals. Missing is
#' 	common in this scheme. However, you may need large window size for 
#' 	reconstruction of the signals. 
#'
#' @details
#'	If \emph{2*m*f} is not an integer, the recovered signal may have  
#'	included a phase shift. However, if the option of "phase shift 
#'	correction is enabled, the related errors caused by the phase shift 
#'	can be limited to an acceptable level. Please notice that this 
#'    method is useful for cases with low background noise and the aim 
#'    is to near perfectly recover the signal. The targeted signal also 
#'    needs to be the dominant signal in the data so that the interaction 
#'    of other signals is negligible.
#'
#'	Another way to reduce the errors caused by unmatched \emph{m} and 
#'	\emph{f} is to use the option of "adaptive window size". It will help 
#' 	you select the best window size for the given frequencies and the 
#'	data length automatically. But it only works well when it exists 
#'	\emph{m} for integer values \emph{2*m*f}.
#'
#'	These two options haven't been implemented for \code{kz.ftc}.  
#'	
#' @param      x 	The data vector. Missing values are allowed.
#' @param      m 	The window size for a regular Fourier transform 
#' @param    ...	Other arguments. 
#' \itemize{
#'  \item	\code{k : } Integer. The iterations number of KZFT.
#'  \item   \code{f : }	Vector. Selected frequencies. Default value is c(1:m)/m
#'  \item	\code{n : }	The sampling frequency rate as a multiplication 
#'				of the Fourier frequencies
#'  \item	\code{p : } The distance between two successive intervals as a
#'				percentage of the total length of the data series.
#'  \item	\code{adpt :} Logic. Flag for using adaptive window size, or not.
#'				Default is FALSE.
#'  \item	\code{phase :} Logic. Flag for correcting phase shift, or not.
#'				Default is FALSE.
#' }
#'
#' @return 		List. It includes data frame for Fourier transform matrix \code{tfmatrix},
#'				column means of Fourier transform matrix \code{fft}, vector \code{pg}
#'				and \code{f} for KZ-periodogram values and corresponding frequencies.
#' @rdname 		kzft
#' @name 		kzft
#' @export		kz.ft
#' @keywords 	KZFT
#' @concept     	Kolmogorov-Zurbenko Fourier Transform
#' @seealso		\code{\link[kzft]{kzft}}, \code{\link{kz.smpg}}, \code{\link{kzp2}}
#' @examples 
#' ## Adapted from kzft::kzp example 2
#'	t <- 1:2000
#'	y <- 1.1*sin(2*pi*0.0339*t)+7*sin(2*pi*0.0366*t)
#'	y2 <- y
#'	noise <- rnorm(length(t),0,1)
#'	y[sample(t,100,replace=FALSE)] <- NA
#'	f <- c(0.0339, 0.0366)
#'
#' ## Periodogram
#'	ft <- kz.ft(y+5*noise, f=f, k=2, m=1000, n=10)
#'	# It may take 10 ~ 20 seconds
#'	# system.time(ft <- kz.ft(y+5*noise, k=2, m=1000, n=10)) 
#'	plot(y=log(ft$pg+1), x=ft$f, type="l", xlim=c(0.025,0.045))
#'	abline(v=f, lty=21, col="red")
#'	text(x=f+0.001, y=c(2,4), f, col="red", cex=0.75)
#'
#' ## recover signal
#'	ft <- kz.ft(y+5*noise, f=f, k=3, m=500)
#'	yr <- 2*Re(rowSums(ft$tf))
#'	cor(yr, y2[1:length(yr)], use="pairwise.complete.obs")
#'	plot((y+5*noise)[1:length(yr)], type="p", cex=0.5, col="grey")
#'	points(y[1:length(yr)],type="b", col="red", cex=0.45)
#'	points(yr, type="p", cex=0.35, col="blue")
#'	mtext("Red dots: singal, Blue dots: reconstruction", cex=0.75)
#'
#'
#' ## Additional example
#'	t <- 1:2000
#'	y <- 1.1*sin(2*pi*0.011*t)+2*sin(2*pi*0.032*t)
#'	y2 <- y
#'	y[sample(t,500,replace=FALSE)] <- NA
#'	noise <- rnorm(2000,0,1)
#'	ft <- kz.ft(y + 3.0*noise, k=5, f=c(0.011,0.032), m=300, adpt=FALSE)
#'	yr <- 2*Re(rowSums(ft$tf))
#'	cor(yr, y2[1:length(yr)], use="pairwise.complete.obs")
#'	plot((y+5*noise)[1:length(yr)], type="p", col="grey")
#'	points(y2[1:length(yr)], type="l", col="red")
#'	points(y[1:length(yr)],type="p", col="red", cex=0.35)
#'	points(yr, type="p", cex=0.3, col="blue")
#'	mtext("Red: singal, Grey: singal + 5*noise, Blue: reconstruction", cex=0.75)
#'
# -------------------------------------------------------------------------------------

kz.ft <- function (x, m, ...) 
{
    data <- as.vector(x)
    dots <- list(...)
    if (hasArg("k")) { k <- dots$k } else { k <- 1 }    
    if (hasArg("n")) { n <- dots$n } else { n <- 1 }    
    if (hasArg("adpt")) { adpt <- dots$adpt } else { adpt <- FALSE }
    if (hasArg("phase")){ phase <- dots$phase } else { phase <- FALSE }
    N <- length(data)
    m <- round(m)
    if (hasArg("f")) { 
	f <- mapply(FUN=min, dots$f, abs(1-dots$f))
	delta <- (2*m*f - floor(2*m*f))/2
	if (adpt & any(delta !=0)) { m <- adapt.m(f, m=m, N=N, k=k) }
	sc <- 1/((m-1)*k+1)
    } else { 
	f <- seq(0, 1, length = n * m + 1)[-1]
	sc <- 1
	if (!hasArg("phase")) { phase <- FALSE }
    }
    if (((m - 1) * k + 1) > length(data)) 
        stop("invalid 'm' & 'k':(m-1)k+1 should be less equal length of data")
    if (hasArg("p")) { p <- dots$p } else { p <- sc }
    M <- (m - 1) * k + 1
    L <- round(M * p)
    T <- floor((N - M)/L) + 1
    kzft <- array(NA, dim = c(T, length(f)))
    omega <- 2 * pi * f 
    s <- 0:(M - 1)
    coef <- kzft::coeff.kzft(m, k)
    coefft <- coef * exp((-(0+1i)) * s %o% omega)
    data[is.na(data)] <- 0
    for (t in (1:T)) {
	tmpc <- coef 
	tmpc[is.na(x[((t - 1) * L + 1):((t - 1) * L + M)])] <- 0
	kzft[t, ] <- data[((t - 1) * L + 1):((t - 1) * L + M)] %*% (coefft / sum(tmpc))
    }
    phsf <- (2*(m - round(2*m*f)/(2*f)) * f)*pi
    delta <- sign(phsf)*(1 - cos(phsf))
    adset <- which(round(delta,10) != 0 & 2*m*f>=1 & dim(kzft)[1]>=(1/f))
    if (length(adset) > 0 & phase & dim(kzft)[1]>6) {
	value <- rep(0,length(adset))
	for (j in adset) {
	   Mok <- ifelse(round(M*f[j]) == M*f[j], M, 2*M)
	   KZFT <- array(NA, dim = c(min(Mok,T), 1))
	   DATA <- sin(2*f[j]*pi*(1:(2*Mok)))
	   for (t in (1:min(Mok,T))) {
		KZFT[t, ] <- DATA[((t - 1) * L + 1):((t - 1) * L + M)] %*% coefft[,j] 
	   }
	   a <- max(abs(KZFT))
	   b <- min(abs(KZFT))
	   ceta <- unique(Arg(KZFT[which(abs(KZFT)==b)]))
	   for (u in 1:3) {
		ellipse <- fit.CE(KZFT[1:min(Mok,T),1], a0=a, b0=b, ceta=ceta)
		a <- ellipse$par[1]
		b <- ellipse$par[2]
		ceta <- ellipse$par[3]
	   }
	   value[j] <- ellipse$value
	   d  <- (a - b)/2 
	   if (d < 0) { 
		ceta <- pi/2 + ceta 
		hash <- kzft[,j] * (cos(-ceta) + sin(-ceta)*1i)
		hash <- Re(hash)*(0.5/(0.5+abs(d))) + Im(hash)*(0.5/(0.5-abs(d)))*1i
		kzft[,j] <- hash * (cos(ceta) + sin(ceta)*1i)
	   }
	}
    }
    kzftf <- colMeans(kzft)
    kzpv <- colMeans(abs(kzft)^2) * M
    kzp  <- rep(0, n * m)
    plc <- mapply(FUN=max, 1, round(n * m * f))
    plc <- mapply(FUN=min, round(n * m / 2), plc)
    kzp[plc] <- kzpv[1:length(f)]
    kzp <- kzp[1:round(n * m / 2)]
    frule <- seq(1/(n * m),1,1/(m * n))
    frule <- frule[1:round(n * m/2)]
    if (length(f) >= length(frule)) {
	f <- ""
    } else {
	if (length(frule[which(kzp>0)])==length(f[order(f)])){
	   frule[which(kzp>0)] <- f[order(f)]
	} else {
	   if (length(f)==1) browser()		# catch error: 2017-1-15
	   vfv <- round(n * m * f[order(f)])[1:(length(f)-1)] == 
			round(n * m * f[order(f)])[2:length(f)]
	   for (j in (1:(length(f)-1))) {
		if (vfv[j]) {
		   lf <- length(frule)
		   c(frule[1:which(kzp>0)[j]], frule[which(kzp>0)[j]], 
				frule[(which(kzp>0)[j]+1):lf]) -> frule
		   frule[which(kzp>0)[j]+(0:1)] <- f[order(f)][j:(j+1)]
		   c(kzp[1:which(kzp>0)[j]], kzp[which(kzp>0)[j]], 
				kzp[(which(kzp>0)[j]+1):lf]) -> kzp
		   kzp[which(kzp>0)[j]+(0:1)] <- kzpv[order(f)][j:(j+1)]
		} else {
		   frule[which(kzp>0)][j] <- f[order(f)][j]
		}
	   }
	}
    }
    pars <- list(f=f, m=m, n=n, k=k, p=p, adpt=adpt, phase=phase, diff=delta)
    if (exists("value")) { 
	pars$value <- value; pars$a <- a; pars$b <- b; pars$ceta <- ceta 
    }
    lst <- list(tfmatrix = kzft, fft = kzftf, f=frule, pg = kzp, pars=pars)
    return(lst)
}


# -----------------------------------------------------------------------------------
#   Kolmogorov-Zurbenko Fourier Transform Function For Irregularly Sampled Data
#
#' @rdname 		kzft
#' @name 		kzft
#' @export		kz.ftc
#' @examples 
#'
#' ## Example for kz.ftc
#'
#'	t <- runif(2000)*2000
#'	f <- c(0.15, 0.1)
#'	x <- sin(2*pi*f[1]*t + pi/4)
#'	y <- sin(2*pi*f[2]*t + pi/12)
#'	y <- y[order(t)]
#'	x <- x[order(t)]
#'	tr <- t[order(t)]
#'	noise <- rnorm(length(tr),0,1)
#'	plot(y=y+x, x=tr, type="l")
#'
#'## Periodogram
#'	ft <- kz.ftc(x+y+2*noise, xt=tr, k=2, m=1000)
#'	plot(y=ft$pg, x=ft$f, type="l")
#'	abline(v=f, col="grey", lty=21)
#'	text(x=f+0.001, y=c(200,400), f, col="red", cex=0.75)
#'	mtext("Spectrum of Longitudinal Data, Selected f")
#'
#'## recover signal
#'	ft <- kz.ftc(x+y+noise, xt=tr, f=f, k=1, m=1900)
#'	yr <- rowSums(2*Re(ft$tf))
#'	iv <- 0:60
#'	plot(y=(x+y+noise)[iv], x=tr[iv], type="p", col="grey")
#'	xt <- (0:8000)/100
#'	yt <- sin(2*pi*f[1]*xt+pi/4) + sin(2*pi*f[2]*xt+pi/12)
#'	y2 <- sin(2*pi*f[1]*iv+pi/4) + sin(2*pi*f[2]*iv+pi/12)
#'	points(yt, x=xt, col="grey", cex=0.5, lwd=1, type="l")
#'	points(y2, x=iv, col="blue", cex=0.75, lwd=1, type="p")
#'	points(y=yr, x=0:(length(yr)-1), type="p", cex=0.5, lwd=1, col="red")
#'	mtext("Red: reconstruction, Grey: signal + noise", cex=0.75)
#'
# -----------------------------------------------------------------------------------

kz.ftc <- function (x, m, ...) 
{
    data <- as.vector(x)
    N <- length(data)
    dots <- list(...)
    if (hasArg("xt")) { xt <- dots$xt } else { xt <- 1:N }    
    if (hasArg("k")) { k <- dots$k } else { k <- 1 } 
    if (((m - 1) * k + 1) > length(data)) 
        stop("invalid 'm' & 'k':(m-1)k+1 should be less equal length of data")
    if (hasArg("n")) { n <- dots$n } else { n <- 1 }    
    if (hasArg("f")) { 
	f  <- dots$f 
	sc <- 1/(m*k) 
     } else { 
	f  <- seq(0, 1, length = n * m + 1)[-1]
	sc <- 1
    }
    if (hasArg("p")) { p <- dots$p } else { p <- sc }
    M <- (m - 1) * k + 1
    L <- round(M * p)
    T <- floor((N - M)/L) + 1
    kzft <- array(NA, dim = c(T, n * length(f)))
    coef <- kzft::coeff.kzft(m, k)
    lbtw <- function(ip, coef) { 
	if (floor(ip)==ceiling(ip)) return(coef[floor(ip)])
	ip <- max(min(ip, length(coef)-1),0.00001)
	coef[floor(ip)+1]*(ceiling(ip)-ip)+coef[ceiling(ip)+1]*(ip-floor(ip)) 
    }
    for (t in (1:T)) {
	sr <- which(xt > ((t - 1) * L ) & xt <= ((t - 1) * L + M ))
	st <- xt[sr] - (t - 1) * L 
	cref <- unlist(sapply(st, FUN=lbtw, coef))
	coefft <- exp(-(0+1i) * (st %o% (2 * pi * f))) * cref
	kzft[t, ] <- data[sr] %*% (coefft / sum(coef))
    }
    kzftf <- colMeans(kzft)
    kzpv <- colMeans(abs(kzft)^2) * M
    kzp  <- rep(0, n * m)
    kzp[ round(n * m * f) ]  <- kzpv[1:length(f)]
    kzp <- kzp[1:round(n * m / 2)]
    frule <- seq(1/(n * m),1,1/(m * n))
    frule <- frule[1:round(n * m/2)]
    lst <- list(tfmatrix = kzft, fft = kzftf, f=frule, pg = kzp)
    return(lst)
}



# -----------------------------------------------------------------------------------
#' @title		
#'      Reconstruct 2D Wave Signals For Given Directions with KZFT
#'
#' @description	
#'    Once you have identified the waves' directions and frequencies, 
#'    you can reconstruct the spatial wave signals with \code{kz.rc2}.
#'	Directional information is utilized to suppressive the noise.
#'
#' @details
#'	Averaging along the orthogonal direction of a wave signal will
#'	significantly reduce the noise effects and increase the accuracy
#'	of reconstruction. 
#'
#'    When the direction information is not available, the 2D signal 
#'	will be reconstructed along x-axis, but the result usually has
#'	a phase-shift even for the dominant wave pattern. 
#' 
#'    If choose to average reconstructed signal along its orthogonal 
#'    direction, \code{rlvl} should be set to control the averaging 
#'    level. If the input data has comparable size on x- or y-dimension,
#'    an experiential formula is to use the integer value around
#'    \code{sqrt(dx)*1.5}, where \code{dx} is the array size.
#'
#' @param     ds 	Matrix for data of wave field. Missing values are allowed.
#' @param      f 	Vector. Identified wave frequency. 
#' @param      m 	The window size for a regular Fourier transform 
#' @param    ...	Other arguments. 
#' \itemize{
#'  \item   \code{angle : } Vector. Identified wave direction value in degree. 
#'  \item	\code{k : } Integer, defaulting to 2. The iterations number of KZFT.
#'  \item	\code{n : }	The sampling frequency rate as a multiplication 
#'				of the Fourier frequencies
#'  \item	\code{p : } The distance between two successive intervals as 
#'				a percentage of the total length of the data series
#'  \item	\code{avg : } Logic. If average along orthogonal direction. 
#'				Default is TRUE if \code{angle} is available. 
#'  \item	\code{plot : } Logic. Flag for outputing figures. Default is FALSE. 
#'  \item	\code{compare : } Logic. Flag for drawing input image. Default is FALSE. 
#'  \item	\code{edge : } Logic. Flag for keeping the edge data in the returned
#'				   reconstructed signal. Default is FALSE. 
#'  \item	\code{rlvl : } Integer. Coefficient to control the averaging level 
#'				when \code{avg} is TRUE. Default value is 2.
#' }
#' @rdname 		kzrc2
#' @name 		kzrc2
#' @export		kz.rc2 
#' @keywords 	KZFT reconstruction
#' @concept     	Kolmogorov-Zurbenko Fourier Transform
#' @seealso		\code{\link[kzft]{kzft}}, \code{\link{kz.smpg}}, \code{\link{kzp2}}
#' @examples 
#'	dx <- 100			# The x and y scale of the wave field
#'	dy <- 100			# Enlarge them to 300 to get better result.
#'
#'	b <- expand.grid(x=1:dx, y=dy:1)
#'	q1 <- pi/6; f1 <- 0.1;
#'	b$v1 <- sin(f1*2*pi*(b$x*cos(q1)+b$y*sin(q1))+runif(1))
#'	a1 <- array(0,c(dx,dy))
#'	a1[as.matrix(b[,1:2])] <- b$v1
#'	q2 <- -pi/3; f2 <- 0.15;
#'	b$v2 <- sin(f2*2*pi*(b$x*cos(q2)+b$y*sin(q2))+runif(1))
#'	a2 <- array(0,c(dx,dy))
#'	a2[as.matrix(b[,1:2])] <- b$v2
#'	a <- array(0,c(dx,dy))
#'	a[as.matrix(b[,1:2])] <- b$v1 + 2.5*b$v2
#'	noise <- matrix(rnorm(dx*dy,0,1),ncol=dy)
#'
#'	persp(1:(dx/2), 1:(dy/2), a1[1:(dx/2), 1:(dy/2)], zlab="",
#'		main="wave #1", theta=0, phi=45, ticktype="detailed", col="lightblue")
#'	persp(1:(dx/2), 1:(dy/2), a2[1:(dx/2), 1:(dy/2)], 
#'		main="wave #2", theta=90, phi=-110, ticktype="detailed", col="lightblue")
#'	persp(1:(dx/2), 1:(dy/2), a[1:(dx/2), 1:(dy/2)], 
#'		main="wave #1 + #2 ", theta=90, phi=-110, ticktype="detailed", col="lightblue")
#'	persp(1:(dx/2), 1:(dy/2), a[1:(dx/2), 1:(dy/2)] + 5*noise[1:(dx/2), 1:(dy/2)], 
#'		main="wave #1 + #2 + 5*noise", theta=90, phi=-110, ticktype="detailed", col="lightblue")
#'
#'	image(x=1:dim(a1)[1] , y=1:dim(a1)[2], z=a1)
#'		box(); mtext("wave #1")
#'	image(x=1:dim(a2)[1] , y=1:dim(a2)[2], z=a2) 
#'		box(); mtext("wave #2")
#'	image(x=1:dim(a)[1] , y=1:dim(a)[2], z=a+0*noise)
#'		box(); mtext("wave #1 + #2 ")
#'	image(x=1:dim(a)[1] , y=1:dim(a)[2], z=a+7*noise)
#'		box(); mtext("wave #1 + #2 + 7*noise")
#'
#'	rc0 <- kz.rc2(a+0*noise, angle=c(q1,q2)*180/pi,f=c(f1,f2), m = 50, avg=FALSE)
#'	cor(as.vector(a[1:dim(rc0)[1],1:dim(rc0)[2]]), as.vector(rc0), use="pairwise.complete.obs")
#'
#'	rc0 <- kz.rc2(a+0*noise, angle=c(q1,q2)*180/pi,f=c(f1,f2), m = 50, avg=TRUE, rlvl=15)
#'	cor(as.vector(a[1:dim(rc0)[1],1:dim(rc0)[2]]), as.vector(rc0), use="pairwise.complete.obs")
#'
#'	rc <- kz.rc2(a+7*noise, angle=c(q1,q2)*180/pi,f=c(f1,f2), m = 50, avg=TRUE, rlvl=15, plot=TRUE)
#'	cor(as.vector(a[1:dim(rc)[1],1:dim(rc)[2]]), as.vector(rc), use="pairwise.complete.obs")
#'	dev.new();image(x=1:dim(rc)[1] , y=1:dim(rc)[2], z=a[1:dim(rc)[1],1:dim(rc)[2]])
#'	box();title("Signal without noise")
#'
#'	rc <- kz.rc2(a+7*noise, angle=q2*180/pi, f=f2, m = 50, avg=TRUE, rlvl=21, plot=TRUE)
#'	cor(as.vector(a2[1:dim(rc)[1],1:dim(rc)[2]]), as.vector(rc), use="pairwise.complete.obs")
#'	dev.new();image(x=1:dim(rc)[1] , y=1:dim(rc)[2], z=a2[1:dim(rc)[1],1:dim(rc)[2]])
#'	box();title("Signal without noise")
#'
# -----------------------------------------------------------------------------------

kz.rc2 <- function(ds, f=0.25, m=round(min(dim(ds)/2)), ...) {
   dots <- list(...)
   if (hasArg("k"))   { k <- dots$k } else { k <- 1  }
   if (hasArg("edge"))   { edge <- dots$edge } else { edge <- FALSE  }
   if (hasArg("plot"))   { plot <- dots$plot } else { plot <- FALSE  }
   if (hasArg("angle"))  {angle <- dots$angle} else { angle <- 0     }
   if (hasArg("compare"))  {compare <- dots$compare} else { compare <- FALSE  }
   if (hasArg("rlvl"))  { rlvl <- dots$rlvl} else { rlvl <- 2  }
   if (hasArg("avg") )   {  
	avg <- dots$avg  
   } else {  
	if (hasArg("angle")) { avg <- TRUE } else { avg <- FALSE }
   }
   agls <- angle%%180 
   agls <- ifelse(agls >= 180, agls%%180, agls)
   agls <- ifelse(agls <= -90, agls%%180, agls)
   agls <- ifelse(agls > 90, agls-180, agls)
   df <- a2d(ds)
   Ang <- abs(90*round(angle/90))
   for (i in 1:length(agls)) {
	F1 <- efg(f[i], agls[i], Ang[i])
	if (F1 > 0.5) { F1 <- (1 - F1) }
  	xy <- getwave(df, Ang[i]*pi/180)[,1:4]
	wv <- split(xy$obs,xy$e)
	xy <- xy[order(xy$e),]
	rv <- list()
  	for (j in 1:length(wv)) {
	   rv[[j]] <- as.vector(rep(NA, length(wv[[j]])))
	   if (length(wv[[j]]) < ((m - 1) * k + 1)) next
  	   xr <- 2*Re(kz.ft(wv[[j]], m=m, f=F1, k=k, ...)$tf) 
	   rv[[j]][1:length(xr)] <- as.vector(xr)	
 	}
	drv <- unlist(rv)
	xy$xr <- drv
	if (sum(diff(xy$x)==-1) > sum(diff(xy$x)==1)) { scx = -1 } else { scx = 1}
	if (sum(diff(xy$y)==-1) > sum(diff(xy$y)==1)) { scy = -1 } else { scy = 1}
	mtxr <- df2mt(xy[,c(1,2,5)], scale=c(1,1))
	if (avg) {
	   xyz <- getwavf(xy[,c(1,2,5)], (agls[i]-90)*pi/180, f[i], rlvl)
	   drv <- aggregate(xyz$xr, by=list(xyz$e), FUN=mean, na.rm=T)
	   names(drv) <- c("e","avg")
	   xyz <- merge(xyz, drv, all=TRUE)
	   mtxr2 <- df2mt(xyz[,c(2,3,5)], scale=c(1,1))
	}
	if (i==1) { 
	   if (avg) { ds2 <- mtxr2 } else { ds2 <- mtxr } 
	} else {
	   if (avg) { ds2 <- ds2 + mtxr2 } else { ds2 <- ds2 + mtxr } 
	}
   }
   if (!edge) {
	M <- (m-1)*k 
   	dx <- dim(ds2)[1]-M
   	dy <- dim(ds2)[2]-M
	ds2 <- ds2[(1:dx)+(1-scx)*M/2, (1:dy)+(1-scy)*M/2]
   }
   dx <- dim(ds2)[1]
   dy <- dim(ds2)[2]
   if (plot) {
	txt1 <- toString(paste(round(angle[i],2), enc2utf8("\xB0"), " ",sep=""))
	txt2 <- paste("f=", toString(paste(round(f,3))),sep="")
	if (compare) {
	   dev.new(); graphics::image(x=1:dx, y=1:dy, z=ds[1:dx,1:dy])
	   graphics::box(); graphics::title("Signal")
	   mtext(paste(txt2, ", d=", txt1, sep=""), cex=0.75, line=0.2)
	}
	dev.new(); graphics::image(x=1:dx, y=1:dy, z=ds2)
	graphics::box(); graphics::title("Reconstruction")
	mtext(paste(txt2, ", d=", txt1, sep=""), cex=0.75, line=0.2)
   }
   return(ds2)
}


