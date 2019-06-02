# -----------------------------------------------------------------------------------
#' @title		
#'	    Average Periodogram for Spatial Data in Given Directions 
#'
#' @description	
#'    Functions in this group are designed to check periodogram for data series 
#' in a given direction or a list of directions. 
#'
#'    \code{kzpdr} samples the data of wave field, and outputs the average
#' pattern of periodogram for series in a given direction. A collection of 
#' these pattern records will be sent to \code{kzpdr.eval} or \code{kzpdr.estimate} 
#' to estimate the wave frequencies and directions.   
#'
#' @rdname 	kzpdr
#' @param     ds 	Data array. Only 2 dimensional arrays are allowed for current version.
#' @param   plot	TRUE or FLASE. Flag for outputting designed periodogram plot or not.
#'                Defaults to FLASE. In \code{kzpdr}, the plot is the mean periodogram 
#'			for data series in a given direction. 
#' @param  angle	Vector or single numeric value in radians.
#' @param   pair	Logic. Defaults to TRUE, i.e., check the given directions and their
#'			orthogonal opposition at the same time.
#' @param 	...	Other arguments. 
#' \itemize{
#'  \item 	For function \code{kzpdr}, it could be the following arguments
#'		(right of equals signs are their default setting):
#'	\itemize{
#'	 \item  \code{w = 20 : } smoothing window size.
#'	 \item  \code{dpct = 0.01 : } a percentage of total variation of periodogram;
#'		  smoothing window is extended until variation within the window 
#'		  reaches this number. See DZ method in \code{kzft::kzp} for details.
#'	 \item  \code{rg = c(0,0.5) : } the frequency range for the periodogram.
#'	 \item  \code{raw = FALSE : } if use the raw periodogram directly. 
#'	 \item  \code{log = FALSE : } if use log scale for periodogram.
#'	 \item  \code{frun = FALSE :} if force to run the sampling on given directions.
#'		  Defaults to check records and not sample on duplicate directions
#'	 \item  \code{min.ln = 0.6 :} the minimum ratio of sampling data length vs. 
#'		  original data length to product a periodogram for a direction.
#'	}	
#'
#'  \item	In \code{kzpdr.3d} function, it could be arguments of 
#'		the perspective plot, like \code{theta}, \code{phi}, etc.,
#'		please refer function \code{graphics::persp} for more 
#'		information.
#'
#'  \item	For \code{kzpdr.valid}, \code{level} control the cross-validation
#'		process: integer number \code{k} means to run cross-validation
#'		by excluding \code{k} pairs of directional samples each time. 
#'		Default value is 1.
#' }
#'
#' @details		
#'		    \code{kzpdr} is used to sample the spatial data and generates 
#'		periodograms in orthogonal direction pairs; the frequencies of spikes 
#'		for each directional periodogram are identified and recorded as the 
#'		function output. The spike pattern of average periodograms for spatial
#'		directions can help to identify wave frequencies and directions.
#'
#'		    Function \code{kzpdr.3d} will provide 3D perspective plot as the   
#'		global view for periodograms of data series in a given direction.
#'
#' @return	    
#'		    The returned data list of function \code{kzpdr} includes the 
#'		data frame for frequencies of spikes on mean periodograms of each
#'		checked direction. It also includes a vector recording the md5sum
#'		value of the spatial wave data array for internal control.
#'
#'		    Function \code{kzpdr} will output the periodogram plots when 
#'		option \code{plot} is set as TRUE. The frequencies of marked spikes 
#'		will also be print out for each sampling direction.
#'			
#'		    \code{kzpdr.3d} returns back the data frame for re-gridded mean 
#' 		periodogram for data series in given direction, as showed in the 
#'		perspective plot.
#'
#' @keywords 	directional-periodogram
#' @concept 	directional periodogram
#' @export
#' @seealso		\code{\link{kzp2}}, \code{\link{kzpdr.tol}}, \code{\link{kzpdr.eval}}
#'   			\code{\link{kzpdr.valid}}, \code{\link{kzpdr.spikes}}
#' @examples
#'	dx <- 300			 
#'	dy <- 300			 
#'
#'	b <- expand.grid(x=1:dx, y=1:dy)
#'	q1 <- pi/3; f1 <- 0.2;
#'	b$v1 <- sin(f1*2*pi*(b$x*cos(q1)+b$y*sin(q1))+100*runif(1))
#'	q2 <- pi/6; f2 <- 0.05;
#'	b$v2 <- sin(f2*2*pi*(b$x*cos(q2)+b$y*sin(q2))+100*runif(1))
#'
#'	a <- array(0,c(dx,dy))
#'	a[as.matrix(b[,1:2])] <- b$v1 + 1.5*b$v2
#'	persp(1:dx, 1:dy, a, theta=90, phi=-110, 
#'		ticktype="detailed", col="lightblue")
#'	a <- a + 5*matrix(rnorm(dx*dy,0,1),ncol=dy)
#'	persp(1:dx, 1:dy, a, theta=90, phi=-110, 
#'		ticktype="detailed", col="lightblue")
#'
#'	# It may take about 30 seconds
#'	# o <- kzpdr.3d(a, -pi/6) 
#'
#'	# Load pre-saved data to save running-time
#'	data(kzpdr.demo);
#'
#'	# sampling, it may take a few minutes 
#'	# system.time(kzpdr.demo <- kzpdr(a, pi/12, pair=FALSE, plot=TRUE)) 
#'	# system.time(kzpdr.demo <- kzpdr(a, pi/12, plot=TRUE))
#'	# system.time(kzpdr.demo <- kzpdr(a, c(0, pi/6, pi/4, pi/3), plot=TRUE))
#'  
#'	kzpdr.spikes(kzpdr.demo)
#'
#'	# For identification of the wave parameters, see kzpdr.estimate
# -----------------------------------------------------------------------------------

kzpdr <- function(ds, angle, plot=F, pair=T, ...) {
   dots <- list(...)
   if (hasArg("dpct")){ dpct <- dots$dpct} else { dpct <- 0.01 }
   if (hasArg("log")) {  log <- dots$log } else { log <- FALSE }
   if (hasArg("raw")) {  raw <- dots$raw } else { raw <- FALSE }
   if (hasArg("frun")) { frun <- dots$frun } else { frun <- FALSE }
   if (hasArg("min.ln")) { min.ln <- dots$min.ln } else { min.ln <- 0.6 }
        if (hasArg("k")) { k <- dots$k } else { k <- 1 } 
        if (hasArg("n")) { n <- dots$n } else { n <- 1 } 
        if (hasArg("w")) { w <- dots$w } else { w <- 20*n }
        if (hasArg("cp")) { cp <- dots$cp } else { cp <- 0 } 
   md5 <- digest::digest(ds, algo="md5")
   rec0 <- kzpdr.rec(ls(1),md5)
   OL <- data.frame(freq = 0, spg = 0, id = 0)[0,]
   OL0 <- OL
   angs <- unique(angle)
   angs <- ifelse(angs >= pi, angs%%pi, angs)
   angs <- ifelse(angs <= (-pi/2), angs%%pi, angs)
   angs <- ifelse(angs > pi/2, angs-pi, angs)
   if (pair) {
	angle <- angs[1]
	for (i in 1:length(angs)) {
	   addg <- ifelse(angs[i] <= 0, (angs[i] - pi/2)%%pi, angs[i] - pi/2)
	   angle <- unique(round(c(angle, angs[i], addg),10))
	}
   }
   agls <- ckGV(as.vector(angle), rec0)
   okag <- setdiff(angle, agls)
   if (frun) {
	if (length(rec0)>0) {
	   rec0 <- rec0[!(round(rec0$ang,6) %in% round(okag,6)),]
	}
	agls <- angle
   } else {
	okstr <- " Periodogram marked at frequency"
	if (length(okag)>0) {
	   okag <- paste(round((180/pi)*okag,2),enc2utf8("\xB0"),sep="")
	   for (i in 1:length(okag)) {
		cat(okag[i],okstr,rec0[rec0[,1]==okag[i],"freq"],"\n")
	   }
	}
   }
   if (length(agls)==0) { return(list(MD5=md5,rec=rec0)) }
   df <- a2d(ds)
   for (i in 1:length(agls)) {
  	wv <- getwave(df, agls[i], cp=cp)
	wv <- split(wv$obs,wv$e)
	lm <- ifelse(abs(tan(agls[i])) >= 1, dim(ds)[2], dim(ds)[1])  
	OL <- OL0
  	for (j in 1:length(wv)) {
	   if (length(wv[[j]]) < (min.ln * lm))  next
	   m <- floor(length(wv[[j]])/(2*k))*2
	   kzp <- kz.ft(wv[[j]], m=m, adpt=F, phase=F, ...)[c(3,4)]
	   OL <- rbind(OL,kzp)
  	}
  	y <- agrid(OL, scale=min(OL[,1]))
	names(y) <- c("freq", "rpg")
	if (hasArg("f") & all(dots$f %in% unique(OL$f))) {
  	   z <- aggregate(OL[OL[,2]!=0,]$pg,by=list(OL[OL[,2]!=0,]$f), FUN=mean)
	   y[y$rpg>0,c(1,2)] <- c(z[,1],z[,2])
	   raw <- TRUE
	} else {
	   y$rpg[dim(y)[1]] <- mean(y$rpg[-dim(y)[1]])
	   if (hasArg("n") & n>1) { y[1:round(1.5*n),2] <- NA }
	}
	if (hasArg("cut")) {
	   cut0 <- dots$cut
	} else {
	   cut0 <- mean(y$rpg, na.rm=TRUE) + 2*sd(y$rpg, na.rm=TRUE)
	}
	if (raw) {
	   y$spg <- y$rpg; dpct <- 0;
	} else {
	   y$spg <- smooth.kzp(y$rpg, dpct=dpct, w=w)
	}
	if (plot) {
	   # dev.new(); 
 	   rec <- smpg.plot(spg=y$spg,freq=y$freq,cut=cut0,
			Title="Mean Periodogram",angle=agls[i], ...)
	} else {
	   rec <- markspikes(x.fq=c(0,y[,1]),y.spm=c(min(y[,2]),y[,2]), cut=cut0,plot=FALSE, ...)
	   okag <- paste(round((180/pi)*agls[i],2),enc2utf8("\xB0"),sep="")
	   if (length(rec)==0) rec=NA
	   rec <- data.frame(direction=okag, freq=rec)
	}
	pgrt <- y[y$freq %in% rec$freq,c("freq","spg")]
	rec <- merge(rec, pgrt)[,c(2,1,3)]
	rec$max <- order(rec$spg, decreasing = TRUE)
	rec$ang <- as.numeric(gsub(enc2utf8("\xB0"), "", rec$direction))*pi/180
	rec$ang <- ifelse(rec$ang<=(-pi/2), rec$ang%%pi, rec$ang)
   	if (!exists("rec_svg")) {
	    rec_svg <- rbind(rec,rec0)
   	} else {
	    rec_svg <- rbind(rec,rec_svg)
   	}
   }
   rec_svg <- unique(rec_svg)
   rec_svg$direction <- rec_svg$ang*180/pi
   return(list(MD5 = md5, rec = rec_svg, cut=cut0, dpct=dpct,
	max = max(y$rpg, na.rm=TRUE), 
	mean=mean(y$rpg, na.rm=TRUE), sd=sd(y$rpg, na.rm=TRUE)))
}



# --------------------------------------------------------------------------
#  	Draw 3D perspective plots for periodogram of a given direction
#
#  View all periodogram for series in a given direction simultaneously.
#
#' @rdname	kzpdr
#' @export
# --------------------------------------------------------------------------

kzpdr.3d <- function(ds, angle, ...) {
   dots <- list(...)
   if (hasArg("min.ln")){ min.ln <- dots$min.ln} else { min.ln <- 0.6 }
   SP <- data.frame(freq = 0, spg = 0, id = 0)[0,]
   lm <- ifelse(abs(tan(angle[1])) >= 1, dim(ds)[2], dim(ds)[1])  
   df <- a2d(ds)
   wv <- getwave(df, angle)
   for (j in 1:length(wv)) {
	if (length(wv[[j]]) < (min.ln * lm))  next
	sp <- kz.smpg(wv[[j]], ...)
	sp$id <- j
	SP <- rbind(SP,sp)
   }
   SP <- SP[,c("freq","id","spg")]
   id.one <- min(SP$id)
   SP$id <- SP$id - id.one + 1
   z <- agrid(SP, scale=c(min(SP$freq), 1), math="min")[,1:2]
   z$freq <- z$freq / min(z$freq)
   a.sp <- array(0, c(max(z$freq), max(z$id)))
   a.sp[as.matrix(z)] <- SP[,3]
   freq <- (1:round(dim(a.sp)[1]/2))/(round(dim(a.sp)[1]/2)*2)
   series <- (1:round(dim(a.sp)[2]/2)) + id.one - 1
   z <- a.sp[(1:round(dim(a.sp)[1]/2)), 1:round(dim(a.sp)[2]/2)]
   par(cex=0.90)
   persp(freq, series, z, 
	main = paste("Periodograms on ", 
	round(180*angle/pi,2), enc2utf8("\xB0"), sep=""),
	theta=-40, phi=25, ticktype="detailed", col="lightblue")
   return(z)
}


