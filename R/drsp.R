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
#' to estimate the wave frequecies and directions.   
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
#'			    The average periodograms for a few pairs of orthogonal spatial
#'			directions can help to identify frequencies and directions of waves.
#'
#'			    First, function \code{kzpdr} samples the spatial data and generates 
#'			periodograms in orthogonal direction pairs, and the frequencies of spikes 
#'			for each directional periodogram are identified and recorded as the output. 
#'
#'				Then, \code{kzpdr.spikes} can be used to summarized the outputs of  
#'			\code{kzpdr}. Function \code{kzpdr.eval} or \code{kzpdr.estimate} all
#'			can be used to estimate the wave parameters (frequencies and directions).
#'			\code{kzpdr.estimate} can provide visualization of the results, and it is
#'			more convenient to use.
#'
#'			    Usually, if noise level is low, periodograms of a few direction pairs 
#'			may provide satisfied results. But when the noise is high, you may need to 
#'			intensively sample on different directions over the spatial data array with 
#'			\code{kzpdr}. Generally speaking, when the number of samples increases,
#'			the estimation will become more stable and reliable. 
#'
#'				Some other functions could be helpful in this procedure. \code{kzpdr.search}
#'			will search the feasible tolerance settings for wave parameter estimation; 
#'			Function \code{kzpdr.valid} can provide some kind of cross-validation information
#'			for the estimation results of \code{kzpdr.eval} and \code{kzpdr.estimate};
#'			Function \code{kzpdr.3d} will provide 3D perspective plot as the global view for  
#'			periodograms of data series in a given direction.
#'
#' @return		    The returned data list of function \code{kzpdr} includes the 
#'			data frame for frequencies of spikes on mean periodograms 
#'			of each checked direction. 
#'
#'			    Both \code{kzpdr.eval} and \code{kzpdr.estimate} will provide suggested 
#' 			wave frequency and direction values. The data frame of detailed estimation 
#'			for each direction are also include in their outputs. Beside these,   
#'			\code{kzpdr.estimate} can generate 3D or 2D plots for the supports of each 
#'			suggested wave on direction-frequency parameter plane. 
#'			
#'			    \code{kzpdr.3d} returns back the data frame for re-gridded mean 
#' 			periodogram for data series in given direction, as showed in the 
#'			perspective plot.
#'
#' @keywords 	directional-periodogram
#' @concept 	Kolmogorov-Zurbenko periodogram
#' @concept 	directional periodogram
#' @concept 	average periodogram
#' @export
#' @seealso		\code{\link[kzft]{kzp}}, \code{\link{kzp2}}
#'
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
#'	# persp(1:dx, 1:dy, a, theta=90, phi=-110, 
#'	#	ticktype="detailed", col="lightblue")
#'	a <- a + 5*matrix(rnorm(dx*dy,0,1),ncol=dy)
#'	# persp(1:dx, 1:dy, a, theta=90, phi=-110, 
#'	#	ticktype="detailed", col="lightblue")
#'
#'	# It may take a few minutes
#'	# o <- kzpdr.3d(a, -pi/6)
#'
#'	# sampling, it may take a few minutes 
#'	# e <- kzpdr(a, c(0, pi/4, pi/3, -pi/3, pi/18), plot=TRUE)
#'
#'	# load pre-saved data to save running-time
#'	data(kzpdr.demo); e <- kzpdr.demo; rm(kzpdr.demo)
#'
#'	# counting spikes
#'	kzpdr.spikes(e)
#'
#'	# search for tolerance 
#'	kzpdr.search(e, t.D = c(1,2,3), t.F = c(0.005))
#'
#'	# estimate the wave parameters
#'	kzpdr.eval(e, t.D = 3, t.F = 0.01)
#'
#'	# visualization
#'	kzpdr.estimate(e)
#
# -----------------------------------------------------------------------------------

kzpdr <- function(ds, angle, plot=F, pair=T, ...) {
   dots <- list(...)
   if (hasArg("w"))   {   w <- dots$m } else { w <- 20 }
   if (hasArg("dpct")){ dpct <- dots$dpct} else { dpct <- 0.01 }
   if (hasArg("log")) {  log <- dots$log } else { log <- FALSE }
   if (hasArg("raw")) {  raw <- dots$raw } else { raw <- FALSE }
   if (hasArg("frun")) { frun <- dots$frun } else { frun <- FALSE }
   if (hasArg("min.ln")) { min.ln <- dots$min.ln } else { min.ln <- 0.6 }
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
	rec0 <- rec0[!(round(rec0$ang,6) %in% round(okag,6)),]
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
	angle <- agls[i]
	atc <- abs(tan(angle))
  	wv <- getwave(df, angle)
	OL <- OL0
	lm <- ifelse(atc >= 1, dim(ds)[2], dim(ds)[1])  
  	for (j in 1:length(wv)) {
	   if (length(wv[[j]]) < (min.ln * lm))  next
  	   sp <- kz.smpg(wv[[j]], ...)
  	   sp$id <- j
  	   sp$dr <- angle*(180/pi)
	   if (log) { sp$spg <- log(sp$spg) }
         OL <- rbind(OL,sp)
  	}
	sc <- min(OL[,1])
  	y <- agrid(OL[,c(1,2,4)], scale=sc)
	names(y) <- c("freq", "rpg", "id")
	y$direction <- angle
	y$rpg[dim(y)[1]] <- mean(y$rpg[-dim(y)[1]])
	if (raw) {
	   y$spg <- y$rpg
	} else {
	   y$spg <- smooth.kzp(y$rpg, dpct=dpct, w=w)
	}
	ok.n <- (min(1/dim(ds)[1], 1/dim(ds)[2]))/sc
	if (hasArg("n") & ok.n>1) { y[1:ok.n,2] <- y[ok.n,2] }
   	cut0 <- mean(y$rpg, na.rm=TRUE) + 2*sd(y$rpg, na.rm=TRUE)
	if (plot) {
	   dev.new(); 
 	   rec <- smpg.plot(spg=y$spg,freq=y$freq,cut=cut0,dpct=dpct,
			Title="Mean Periodogram",angle=angle)
	} else {
	   rec <- markspikes(x.fq=c(0,y[,1]),y.spm=c(min(y[,2]),y[,2]),plot=FALSE)
	   okag <- paste(round((180/pi)*angle,2),enc2utf8("\xB0"),sep="")
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
   return(list(MD5 = md5, rec = rec_svg))
}


# -----------------------------------------------------------------------------
#     Check spikes on directional periodograms for waves' directions
#
#' @rdname kzpdr
#' @param  rec	  Data frame or list of the outputs of function \code{kzpdr}. 
#'			  Includes the marked frequency values and corresponding directions.
#'			  Defaults is searching for the available records in the environment. 
#' @export
# -----------------------------------------------------------------------------

kzpdr.spikes <- function(rec = ls(1)) {
   fok <- kzpdr.pairs(rec)
   rc <- fok$rc; fc <- fok$fc; 
   ag1 <- fok$ag1; ag2 <- fok$ag2;
   cat("\nChecked", dim(fc)[1], "directions in", length(rc),"orthogonal pairs:\n\n")
   print(noquote(fc$direction)); cat("\n")
   tmpv <- data.frame(cbind(ag1,ag2))
   tmpv$v <- paste("(",ag1,enc2utf8("\xB0"), ", ",ag2,enc2utf8("\xB0"),")", sep="")
   print(noquote(tmpv$v)); 
   fg <- split(fc[,1], fc[,2])
   wvnmbr <- round(fivenum(fc$count)[2:4])
   if (length(sys.parents())>1) { cat("\n"); return(wvnmbr) }
   if (all(wvnmbr == wvnmbr[2])) {
	wvn <- paste(wvnmbr[2])
   } else {
	wvn <- paste(range(wvnmbr)[1],range(wvnmbr)[2],sep=" to ")
   }
   if (length(fg) > 1) {
	angles <- as.numeric(gsub(enc2utf8("\xB0"), "", fg[[1]])) - 90
	angles <- ifelse(angles < -90 | angles > 180, angles%%180, angles)
	cat("\nSuggested wave direction: ", format(paste(angles,enc2utf8("\xB0"),sep=""),width=4), "\n")
	cat("\nSuggested wave numbers: ", format(wvn, width=6), "\n\n")
      mypst <- function(x){ y=x[1];for (i in 1:length(x)) y=paste(y,x[i],sep=", "); 
			ifelse(length(x)<=1,x,y)}
	cat("Note: The spike counts for direction", mypst(fg[[1]]))
	cat(" is only", names(fg)[1], "\n")
	cat("      Other directions have", paste(names(fg)[-1],"spikes"))
	if (length(fg) > 2) cat(" or more")
	cat("\n\n")
   } else {
	cat("\nNo conclusion for wave directions.\n")
	cat("\nSuggested wave numbers: ", wvn, "\n\n")
	cat("Note: All directions have", names(fg)[1])
	cat(" spike(s) on their periodograms. \n\n")
   }
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


# -----------------------------------------------------------------------------
#     Evaluate directional spectrum data for waves' frequencies and directions
#
#' @param	t.D	Tolerance of direction (in degree).
#' @param	t.F	Tolerance of frequency.
#' @rdname kzpdr
#' @export
# -----------------------------------------------------------------------------

kzpdr.eval <- function(rec = ls(1),	t.D = 2.5, t.F = 0.01, ...) {
   dots <- list(...)
   if (hasArg("itr")){ itr <- dots$itr} else { itr <- 1 }
   if (length(sys.parents())>1) { mute = TRUE } else { mute = FALSE }
   if (is.data.frame(rec)) {
	tmpv <- c("dir","g1","g2","sf","agp.2")
	if (all(names(rec)[1:5] %in% tmpv)) {
	   cmbf <- rec; rm(rec);
	}
   }
   if (!exists("cmbf")) {
	if (!mute) wvnmbr <- kzpdr.spikes(rec)
	cmbf <- kzpdr.proj(rec)
	rm(rec)
	cat("\n")
   }
   loopct <- 1; loopct2 <- 1;
   allpairs <- unique(cmbf$dir)
   if (length(allpairs)<3) {
	cat("No enough data!\n\n")
	return(list(df=cmbf[0,]))
   }
   sv.cmbf <- cmbf
   cmbf <- cmbf[order(cmbf$dir, cmbf$agp.2, cmbf$sf),]
   while (TRUE) {
   	cmbf$sd <- shade(cmbf$apha, cmbf$g1)
	t.f <- max(tight(cmbf$sf, cmbf$g1, t.F), t.F, na.rm = TRUE)
	coref <- commcore(split(cmbf$sf,cmbf$dir), STEP = t.f, DIFF=cmbf$sd, digit=2)
	if (length(coref$suggest)==0) return(list(df=cmbf[0,]));
	cmbf$grp.f <- coref$grp
	cmbf <- cmbf[cmbf$grp.f %in% coref$suggest,]
	if (dim(cmbf)[1]<=1) return(list(df=cmbf[0,]))
   	cmbf$sd <- shade(cmbf$apha, cmbf$g1)
  	t.d <- max(tight(cmbf$agp.2*180/pi, cmbf$g1, t.D), t.D, na.rm = TRUE)
   	cored <- commcore(split(cmbf$agp.2*180/pi,cmbf$dir), STEP = t.d, DIFF=cmbf$sd)
	cmbf$grp.d <- cored$grp
	while (length(cored[[1]])==0) {
	   if (!mute) cat("Cannot find expected spike in at least one direction. \n")
	   return(list(df=cmbf[0,]))
	}
	cmbf <- cmbf[cmbf$grp.d %in% cored$suggest,]
	if (dim(cmbf)[1]<=1) return(list(df=cmbf[0,]))
	cmbf$dif.f <- abs(cmbf$sf - cmbf$grp.f)
	cmbf$dif.d <- abs(cmbf$agp.2*180/pi - cmbf$grp.d)
	cmbf <- cmbf[order(cmbf$dir, cmbf$grp.f, cmbf$grp.d, cmbf$dif.f, cmbf$dif.d),]
	cmbf$dif.gf <- c(999, diff(cmbf$grp.f))
	cmbf$dif.gd <- c(999, diff(cmbf$grp.d))
	cmbf$dif.g1 <- c(999, diff(cmbf$g1))
	cmbf$dif.g2 <- c(999, diff(cmbf$g2))
	cmbf$c1 <- cmbf$dif.gd==0 & cmbf$dif.gf==0 & cmbf$dif.g1==0 & cmbf$dif.g2==0
	cmbf <- cmbf[,-c(12:15)]
	tmpw <- unique(cmbf[cmbf$c1, c("dir", "grp.f", "grp.d")])
	cmbf$tv <- paste(cmbf$dir, cmbf$grp.f, cmbf$grp.d)
	tmpw$tv <- paste(tmpw$dir, tmpw$grp.f, tmpw$grp.d)
	cmbf$tw <- cmbf$tv %in% tmpw$tv
	twf <- split(cmbf$sf, f=cmbf$tv)
	twf <- sapply(twf, FUN="max")-sapply(twf, FUN="min")
	twd <- split(cmbf$agp.2, f=cmbf$tv)
	twd <- sapply(twd, FUN="max")-sapply(twd, FUN="min")
	cmbf$c2 <- ifelse(cmbf$c1, twf < t.F & twd < t.D, NA)
	cmbf <- cmbf[!(cmbf$c1 & cmbf$c2), ]
	cmbf <- cmbf[,c("dir","g1","g2","sf","agp.2","apha","grp.d", "grp.f")]
	if (dim(cmbf)[1]<=1) return(list(df=cmbf[0,]))
	loopct <- loopct + 1
	loopct2 <- loopct2 + 1
      if (loopct > 3) break
	if (any(is.na(cored$grp))| any(is.na(coref$grp))) { loopct <- loopct - 1 }
   }
   if (dim(cmbf)[1]>1) {
	cmbf$ok <- sign(cmbf$grp.d)*cmbf$grp.f + 1000*cmbf$grp.d
   	cmbf$sd <- shade(cmbf$apha, cmbf$g1)
	core_ok <- commcore(split(cmbf$ok,cmbf$dir), STEP = 0.0002, DIFF=cmbf$sd, digit=2)
	if (length(core_ok[[1]])==0) {
	   if (itr<=3) return(kzpdr.eval(cmbf[,1:9], itr=itr+1))
	   return(list(df=cmbf[0,1:9]))
	} else { 
	   cmbf <- cmbf[cmbf$ok %in% core_ok$suggest,1:9] 
	}
	if (dim(cmbf)[1]<2) { 
	   if (!mute) cat("Cannot find expected spike in at least one direction. \n")
	   return(list(df=cmbf[0,]))
	}
	cmbf$sd <- shade(cmbf$apha, cmbf$g1)
	t.f <- min(tight(cmbf$sf, cmbf$g1, t.F)*2, t.F, na.rm = TRUE)*1.0005
  	t.d <- min(tight(cmbf$agp.2*180/pi, cmbf$g1, t.D)*2, t.D, na.rm = TRUE)*1.0005
	cored <- commcore(split(cmbf$agp.2*180/pi,cmbf$dir), STEP = t.d, digit=8, DIFF=cmbf$sd)
	coref <- commcore(split(cmbf$sf,cmbf$dir), STEP = t.f, digit=8, DIFF=cmbf$sd)
	if (!(length(coref$suggest)==0 | length(cored$suggest)==0)) {
	   cmbf$grp.f2 <- coref$grp
	   cmbf$grp.d2 <- cored$grp
	   cmbf <- cmbf[cmbf$grp.f2 %in% coref$suggest,]
	   cmbf <- cmbf[cmbf$grp.d2 %in% cored$suggest,]
	} else if (length(cored$suggest) > 0) {
	   cmbf$grp.d2 <- cored$grp
	   cmbf <- cmbf[cmbf$grp.d2 %in% cored$suggest,]
	   tmpv <- as.vector(sapply(split(cmbf$sf,cmbf$grp.d2),FUN=mean))
	   tmpw <- lapply(lapply(split(cmbf$sf,cmbf$grp.d2),FUN=is.na), "!")
	   cmbf$grp.f2 <- rep(tmpv, sapply(tmpw, sum))
	} else if (length(coref$suggest) > 0) {
	   cmbf$grp.f2 <- coref$grp
	   cmbf <- cmbf[cmbf$grp.f2 %in% coref$suggest,]
	   tmpv <- as.vector(sapply(split(cmbf$agp.2,cmbf$grp.f2),FUN=mean))
	   tmpw <- lapply(lapply(split(cmbf$agp.2,cmbf$grp.f2),FUN=is.na), "!")
	   cmbf$grp.d2 <- rep(tmpv, sapply(tmpw, sum))
	   cmbf$grp.d2 <- cmbf$grp.d2*180/pi
	} else {  return(list(df=cmbf[0,])) }
	cmbf$grp.d <- round(cmbf$grp.d2)
	cmbf$grp.f <- round(cmbf$grp.f2,2)
	t.d <- cored$tolerance
	t.f <- coref$tolerance
   }
   nms <- c(names(cmbf)[c(1,4:6)],"grp.f", "grp.d")
   cmb <- data.frame(cmbf[order(cmbf$grp.d, cmbf$grp.f),nms],row.names = NULL)
   names(cmb)[1:4] <- c("Pair","Frequency","Direction","degree")
   cmb$degree <- round(cmb$Direction*180/pi,2)
   suggestion <- unique(cmbf[,c("grp.d","grp.f")])
   leftover <- setdiff(allpairs, unique(cmb$Pair))
   left.df <- sv.cmbf[sv.cmbf$dir %in% leftover, -c(2:3,7)]
   names(left.df)[1:4] <- c("Pair","Frequency","Direction","degree")
   if (length(unique(cmb$Pair))>2) {
	if (length(unique(cmb$grp.d))<18) {
	   cat("Detected direction :",format(paste(suggestion$grp.d,enc2utf8("\xB0"), sep=""),width=6),"\n")
	   cat("Detected frequency :",format(paste(suggestion$grp.f), width=6),"\n\n")
	}
      suggestion <- unique(cmbf[,c("grp.d2","grp.f2")])
   }
   if (!mute) cat("  t.d =", round(t.d,4), "   t.f =", round(t.f,6), "\n\n")
   rng.dif <- function(x) { max(x) - min(x) }
   tof <- aggregate(x=cmb$Frequency, by=list(cmb$grp.f), FUN=rng.dif)
   tod <- aggregate(x=cmb$degree, by=list(cmb$grp.d), FUN=rng.dif)
   ok <- list(direction=unlist(suggestion[,1]), 
		   frequency=unlist(suggestion[,2]), df=cmb,
		   support=sort(unique(cmb$Pair)), 
		   all.pairs=sort(unique(sv.cmbf$dir)), unused=NA,
		   tolerance=data.frame(freq=max(tof$x), 
		   direction=max(tod$x), row.names = "") )
   if (length(leftover)>0) {
	ok$unused <- leftover
	tmp <- data.frame(left.df, row.names = NULL)
	# if (dim(tmp)[1]<5*length(leftover)) ok$exclude.df <- tmp
   } else {
	ok <- ok[-which(names(ok)=="unused")]
   }
   if (!mute) {
   if (!max(length(ok$direction),length(ok$frequency)) %in% wvnmbr){
	cat("  Note: the expected wave number is", unique(wvnmbr), "\n\n\n")
   } }
   return(ok)
}


# -----------------------------------------------------------------------------
#   Search for appropriate tolerances setting
#
#' @rdname kzpdr
#' @export
# -----------------------------------------------------------------------------

kzpdr.search <- function(rec = ls(1), t.D = seq(1,10,1), t.F = 0.01) {
   rec <- kzpdr.rec(rec)
   wvnmbr <- kzpdr.spikes(rec)
   fok <- kzpdr.pairs(rec)
   ag1 <- fok$ag1
   ag2 <- fok$ag2
   cmbf <- kzpdr.proj(rec)
   if (length(ag1)<3) {
	cat("No enough data!\n\n")
	return("")
   }
 for (t_F in t.F) {
   for (t_D in t.D) {
   ort <- kzpdr.eval(cmbf, t.D=t_D, t.F=t_F)
   if (length(unique((ort$df)$Pair))>2) {
	if (length(unique((ort$df)$grp.d))>12) {
	   cat(">12 motion directions mixed together. No clear suggestion.\n\n")
	}
	t.d <- ort$tolerance[1,2]
	t.f <- ort$tolerance[1,1]
	alldf <- ort$df
   } else { 
	t.d <- t.D
	t.f <- t.F
	alldf <- ort$df[0,]
   }
   if (length(ort$frequency)==0) {
	rof=NA 
   } else {
   	rof <- round(ort$frequency,2)
   }
   if (length(ort$direction)==0) {
	rod=NA
   } else {
	rod <- round(ort$direction,0)
   }
   vnm <- max(length(rof), length(rod))
   sumrec0 <- data.frame(td=t_D, tf=t_F, rf=rof, rd = rod, nm = vnm)
   if (length(ort$tolerance$freq)==0) { 
	sumrec0$tlf <- NA
   } else {
	sumrec0$tlf <- ort$tolerance$freq
   }
   if (length(ort$tolerance$direction)==0) {
	sumrec0$tld <- NA
   } else {
	sumrec0$tld <- ort$tolerance$direction
   }
   sumrec0$order <- order(sumrec0$rf, sumrec0$rd)
   if (exists("sumrec")) {
	sumrec <- rbind(sumrec0, sumrec)
   } else {
	sumrec <- sumrec0
   }
   if (is.na(rof[1])) {
	cat("Search stopped at t.D =",t_D,", t.F =",t_F,"\n\n")
	break
   }
   cat("  t.D =", format(t_D, width=6) , "    t.F =", t_F, "\n")
   cat("  t.d =", format(sumrec0$tld[1], width=6), "    t.f =", 
		round(sumrec0$tlf[1],6), "\n\n")
 } }
   srch <- sumrec[sumrec$nm %in% wvnmbr,]
   srch <- srch[order(srch$order, srch$rd),]
   tmp0 <- split(srch$rd, srch$order)
   mydif <- function(x) { return(x-median(x)) }
   tmp1 <- unlist(lapply(tmp0, FUN=mydif))
   tmp1[abs(tmp1) > min(srch$td)] <- 0
   srch$grd <- unlist(tmp0) - tmp1
   srch$grd <- paste(srch$grd,enc2utf8("\xB0"),sep="")
   # ----------------------------------------
   srch <- srch[order(srch$order, srch$rf),]
   tmp0 <- split(srch$rf, srch$order)
   tmp1 <- unlist(lapply(tmp0, FUN=mydif))
   tmp1[abs(tmp1) > min(srch$tf)] <- 0
   srch$grf <- unlist(tmp0) - tmp1
   # ---------------------------------------- 
   mypst <- function(x){ y="";for (i in 1:length(x)) y=paste(y,x[i],sep=" "); y}
   sref <- aggregate(srch[,c("grd","grf")], by=list(srch$tf, srch$td), FUN=mypst)
   final.grp = split(x=sref[,1:2], f=list(sref$grf,sref$grd))
	rule = row.names(sref)
	mark = rule
	for (i in 1:length(final.grp)) {
	   mark[rule %in% row.names(final.grp[[i]])] <- i
	}
   sref <- sref[order(mark),]
   names(sref) <- c("t.F", "t.D", "Directions","Frequencies")
   mycat <- function(x, mk) {
      myline <- " ----------------------------------------------------------------"
      myline <- paste(myline, "---------------------------------------------------", sep="")
	y <- diff(as.numeric(mk))
	y <- c(0, y)
      for (i in 1:dim(x)[1]){
	  if (i == 1) { 
	    w = rep(1,dim(x)[2])
	    for (j in 1:dim(x)[2]){ 
		w[j] <- max(nchar(names(x)[j]), nchar(x[1,j]))+2
		cat(format(names(x)[j], width=w[j], zero.print = "", justify="right"),"  ") 
	    }
	    cat("\n")	
	    cat(substr(myline,1,4+sum(w+2))," \n")
	  }
	  if (abs(y[i])>0) cat(substr(myline,1,4+sum(w+2))," \n")
	  for (j in 1:dim(x)[2]){ 
		cat(format(x[i,j], width=w[j], zero.print = "", justify="right", trim=FALSE),"  ") 
	  }
	  cat("\n")
	}
   }
   cat("\nFeasible tolerance settings in searching range:\n\n")
   cat("\n"); mycat(sref,mark); cat("\n");
}


# -----------------------------------------------------------------------------
#    Cross-validation procedure to test noise effects  
#
#' @rdname kzpdr
#' @export
# -----------------------------------------------------------------------------

kzpdr.valid <- function(rec = ls(1), t.D = 3, t.F = 0.01, ...) {
   dots <- list(...)
   if (hasArg("level")){ level <- dots$level} else { level <- 1 }
   rec <- kzpdr.rec(rec)
   wvnmbr <- kzpdr.spikes(rec)
   fok <- kzpdr.pairs(rec)
   ag1 <- fok$ag1
   ag2 <- fok$ag2
   cmbf <- kzpdr.proj(rec)
   if (length(ag1)<4) {
	cat("No enough data for cross-validation!\n\n")
	return("")
   }
   ort <- kzpdr.eval(cmbf, t.D=t.D, t.F=t.F)
   if (length(unique((ort$df)$Pair))>2) {
	if (length(unique((ort$df)$grp.d))>12) {
	   cat(">12 motion directions mxied together. No clear suggestion.\n\n")
	}
	t.d <- ort$tolerance[1,2]
	t.f <- ort$tolerance[1,1]
	alldf <- ort$df
   } else { 
	t.d <- t.D
	t.f <- t.F
	alldf <- ort$df[0,]
   }
   cat("  t.D =", format(t.D, width=4, justify="right"), "    t.F =", t.F, "\n")
   cat("  t.d =", format(t.d, width=4, justify="right"), "    t.f =", round(t.f,6), "\n\n")
   rof <- unique(ort$df[,c("grp.f","grp.d")])$grp.f
   rod <- unique(ort$df[,c("grp.f","grp.d")])$grp.d
   vnm <- max(length(rof), length(rod))
   if (length(rof)==0)  rof=NA 
   if (length(rod)==0)	rod=NA
   sumrec0 <- data.frame(td=t.D, tf=t.F, rf=rof, rd = rod, nm = vnm, k=0, i=0)
   if (length(ort$tolerance$freq)==0) { 
	sumrec0$tlf <- NA
   } else {
	sumrec0$tlf <- round(ort$tolerance$freq, 6)
   }
   if (length(ort$tolerance$direction)==0) {
	sumrec0$tld <- NA
   } else {
	sumrec0$tld <- round(ort$tolerance$direction, 6)
   }
   sumrec0$order <- order(sumrec0$rf, sumrec0$rd)
   sumrec0$exclude <- NA
   if (exists("sumrec")) {
	sumrec <- rbind(sumrec0, sumrec)
   } else {
	sumrec <- sumrec0
   }
   lvl <- max(min((length(ag1)-3), level),1)
   mypst <- function(x){ y="";for (i in 1:length(x)) y=paste(y,x[i],sep=""); y}
   for (k in lvl:lvl) {
	cat("\n------------------ Cross-validation ----------------- \n")
	notin <- combn(1:length(ag1), k)
	for (i in 1:dim(notin)[2]) {
	   selected <- paste("(",ag1[-notin[,i]],", ",ag2[-notin[,i]],")",sep="")
	   unselected <- paste("(",ag1[notin[,i]],enc2utf8("\xB0, "),ag2[notin[,i]],enc2utf8("\xB0)"),sep="")
	   cat("\n", k, ".", i," Exclude", unselected, "\n\n")
	   ort <- kzpdr.eval(cmbf[cmbf$dir %in% selected,], t.D=t.D, t.F=t.F)
	   rof <- unique(ort$df[,c("grp.f","grp.d")])$grp.f
	   rod <- unique(ort$df[,c("grp.f","grp.d")])$grp.d
	   vnm <- max(length(rof), length(rod))
	   sumrec0 <- data.frame(td=t.D, tf=t.F, rf=rof, rd = rod, nm = vnm, k=k, i=i)
	   sumrec0$tlf <- round(ort$tolerance$freq, 6)
	   sumrec0$tld <- round(ort$tolerance$direction, 6)
   	   sumrec0$order <- order(sumrec0$rf, sumrec0$rd)
	   sumrec0$exclude <- mypst(unselected)
	   sumrec  <- rbind(sumrec, sumrec0)
	   if (length(ort$tolerance)>0) {
	   	cat("  t.d =", format(sumrec0$tld[1], width=4, justify="left"), 
			"   t.f =", round(sumrec0$tlf[1],6), "\n")
	   }
	   if (length(ort$direction)>8) {
		cat("Couldn't give clear suggestion.\n\n")
		next
	   }
	   if (dim(ort$df)[1]==0) { 
		# cat("No suggestion.\n\n")
		next 
	   }
	   if (length(unique((ort$df)$Pair))>2) {
		ort$unused <- c(ort$unused, unselected)
		alldf <- rbind(alldf, ort$df)
	   }
	}
   }
   cat("\n\n-------------------- Summary -------------------- \n\n")
   srch <- sumrec[sumrec$nm %in% wvnmbr,]
   support <- unique(srch[,c("i","k")])
   consistent = FALSE
   if (dim(support)[1]<2) {
      srch <- sumrec
	cat("\n  No consistent result for cross-validation. ")
	cat("\n  Please try different tolerance settings, ")
	cat("\n  or collect more sampling. \n\n")
	t.d <- max(t.D, t.d) + 1
   } else {
   	srch <- srch[order(srch$order, srch$rd),]
   	tmp0 <- split(srch$rd, srch$order)
   	mydif <- function(x) { return(x-median(x)) }
   	tmp1 <- unlist(lapply(tmp0, FUN=mydif))
   	tmp1[abs(tmp1) > min(srch$td)] <- 0
   	srch$grd <- unlist(tmp0) - tmp1
   	srch$grd <- paste(srch$grd,enc2utf8("\xB0"),sep="")
   	# ----------------------------------------
   	srch <- srch[order(srch$order, srch$rf),]
   	tmp0 <- split(srch$rf, srch$order)
   	tmp1 <- unlist(lapply(tmp0, FUN=mydif))
   	tmp1[abs(tmp1) > min(srch$tf)] <- 0
   	srch$grf <- unlist(tmp0) - tmp1
   	# ---------------------------------------- 
   	mypst <- function(x){ y="";for (i in 1:length(x)) y=paste(y,x[i],sep=" "); y}
   	srch <- srch[order(srch$exclude),]
	srch$id <- rownames(srch)
	srch[is.na(srch$exclude),"exclude"] <- ""
   	tmp0 <- aggregate(srch[,c("grd","grf","id")],by=list(srch$exclude), FUN=mypst)
   	tmp1 <- aggregate(srch[,c("tld","tlf")],by=list(srch$exclude), FUN="min")
	tmp0 <- merge(tmp0, tmp1)
   	sref <- aggregate(tmp0[,c("tld","tlf")],by=tmp0[,c("grd","grf")], FUN="min")
   	scnt <- aggregate(tmp0[,c("tlf")],by=tmp0[,c("grd","grf")], FUN="length")
   	names(scnt)[3] <- "len"
	sref <- merge(sref, scnt)
   	scnt <- aggregate(tmp0$id,by=tmp0[,c("grd","grf")], FUN="paste")
   	names(scnt)[3] <- "ID"
	sref <- merge(sref, scnt)
	# ~~~~~~~~~~~~~~~~~~~
	tmp0 <- srch[which(srch$tlf %in% sref$tlf),c("grd","grf","id","exclude")]
	tmp1 <- srch[which(srch$tld %in% sref$tld),c("grd","grf","id","exclude")]
	scnt <- aggregate(tmp0[,c("id","grd","grf")], by=list(tmp0[,c("exclude")]), FUN=mypst)
	tmp0 <- merge(sref,scnt)
	scnt <- aggregate(tmp1[,c("id","grd","grf")], by=list(tmp1[,c("exclude")]), FUN=mypst)
	tmp1 <- merge(sref,scnt)
	lenf <- tmp0[tmp0$len>1, c("len","id")]
	lend <- tmp1[tmp1$len>1, c("len","id")]
	lenf$id <- substr(lenf$id,2,nchar(lenf$id))
	lend$id <- substr(lend$id,2,nchar(lend$id))
	lenf$id <- strsplit(lenf$id," ")
	lend$id <- strsplit(lend$id," ")
	SGNf <- srch[srch$id %in% as.numeric(unlist(lenf$id)),-c(3:6)]
	SGNd <- srch[srch$id %in% as.numeric(unlist(lend$id)),-c(3:6)]
	SGNx <- unique(rbind(SGNf, SGNd))
	lenx <- unique(rbind(lenf, lend))
	if (dim(lenx)[1]==0 | dim(SGNx)[1]==0) {
	   cat("\n\n  Results are inconsist on cross-validation level", lvl,".\n\n")
	} else {
	mylen <- function(id) { 
	   for (j in 1:dim(lenx)[1]) { if (id %in% lenx$id[[j]]) return(lenx$len[j]) }
	}
	SGNx$len <- lapply(SGNx$id, FUN=mylen)
	SGN0 <- srch[srch$k==0,-c(3:6)]
	if (dim(SGN0)[1]==0) {
	   consistent = FALSE
	} else {
	   consistent = TRUE
	   SGN0$len <- 1
	   SGNx <- rbind(SGNx, SGN0)
	}
	sref <- SGNx[,c(9,8,7,4,5,11,6)]
	names(sref)[1:6] <- c("Frequency","Direction","Excluded Pair","t.f","t.d","Times")
	sref[sref$order>1,c(4:6)] <- 0
	sref[sref$order>1,c(3)] <- ""
      mycat <- function(x, mk) {
         myline <- " ----------------------------------------------------------------"
         myline <- paste(myline, "---------------------------------------------------", sep="")
	   y <- diff(as.numeric(mk));  y <- c(0, y)
         for (i in 1:dim(x)[1]){
		if (i == 1) { 
    	    	   w = rep(1,dim(x)[2])
    	    	   for (j in 1:dim(x)[2]){ 
    		   	w[j] <- max(nchar(names(x)[j]), nchar(x[1,j]))+2
    		   	cat(format(names(x)[j], width=w[j], zero.print = "", justify="right"),"  ") 
		   }
		   cat("\n")	
		   cat(substr(myline,1,4+sum(w+2))," \n")
    	  	}
    	  	if (abs(y[i])>0) cat(substr(myline,1,4+sum(w+2))," \n")
    	  	for (j in 1:dim(x)[2]){ 
		   cat(format(x[i,j], width=w[j], zero.print = "", justify="right", trim=FALSE),"  ") 
    	  	}
    	  	cat("\n")
	   }
      }
      cat("\n"); mycat(sref[,-7],SGNx$i); cat("\n");
	cat("  Note: t.f = tolerance of frequency,  t.d = tolerance of direction (in degree)\n")
	cat("\n  Note: the expected wave number is ", unique(wvnmbr), "\n\n")
	# ~~~~~~~~~~~~~~~~~~~
	scnt <- split(sref[,1], f=SGNx$i)
	c1 <- all(unlist(lapply(scnt, FUN="==", scnt[[1]])))
	scnt <- split(sref[,2], f=SGNx$i)
	c2 <- all(unlist(lapply(scnt, FUN="==", scnt[[1]])))
   	if (!(c1 & c2)) {
	   cat("    It may need to be re-grouped based on new tolerance setting,\n")
         cat("    or collect more sampling. \n\n")
   	} else {
	   if (consistent) {
		cat("    Results are consistent on cross-validation level", lvl, ".\n\n")
	   } else {
		cat("    Error may exist in the excluded directions. Or it may\n")
		cat("    need to be re-grouped based on new tolerance setting.\n\n")
	   }
	}
   }}
}


# -----------------------------------------------------------------------------
#     Evaluate directional spectrum data for waves' frequencies and directions
#
#' @param  scale	  The scale of gridding data. 
#' @rdname kzpdr
#' @export
# -----------------------------------------------------------------------------

kzpdr.estimate <- function(rec = ls(1), scale=c(0.005, 1), ...) {
   if (is.data.frame(rec)) {
	tmpv <- c("dir","g1","g2","sf","agp.2")
	if (all(names(rec)[1:5] %in% tmpv)) {
	   cmbf <- rec
	}
   }
   if (!exists("cmbf")) {
	cmbf <- kzpdr.proj(rec)
	cat("\n")
   }
   dots <- list(...)
   if (hasArg("sd"))  { sd  <- dots$sd  } else { sd  <- 0 } 
   if (hasArg("tol")) { tol <- dots$tol } else { tol <- 1 }
   if (hasArg("D3"))  { D3  <- dots$D3  } else { D3  <- FALSE }
   if (hasArg("raw")) { raw <- dots$raw } else { raw <- TRUE  }
   cmbf$sd <- shade(cmbf$apha, cmbf$g1, cc=6) 
   cmbf <- cmbf[order(cmbf$dir, cmbf$agp.2, cmbf$sf),]
   allpairs <- sort(unique(cmbf$dir))
   if (length(allpairs)<3) {
	cat("No enough data!!\n\n")
	return("")
   }
   dc <- max(length(allpairs)-sd, 1)
   mt <- cmbf[order(cmbf$agp.2, cmbf$sf),c("sf","agp.2","sd","dir")]
   mt$id <- as.numeric(row(mt)[,1])
   mt$agp.2 <- mt$agp.2*180/pi
   names(mt)[1:4] <- c("freq","direction","sd","pair")
   mt <- closure(mt, 3, 2*scale)
   	mt2 <- mt[mt$closest,-which(names(mt)=="grpd")]
   mt <- closure(mt, 4, 2*scale)
   	mt3 <- mt[mt$closest,-which(names(mt)=="grpd")]
	mt2 <- unique(rbind(mt2, mt3))
   mt <- closure(mt, 5, 2*scale)
   	mt3 <- mt[mt$closest,-which(names(mt)=="grpd")]
	mt2 <- unique(rbind(mt2, mt3))
   mt <- closure(mt, 10, 3*scale)
   	mt3 <- mt[mt$closest,-which(names(mt)=="grpd")]
	mt2 <- unique(rbind(mt2, mt3))
   mt <- closure(mt, 20, 3*scale)
   	mt3 <- mt[mt$closest,-which(names(mt)=="grpd")]
	mt2 <- unique(rbind(mt2, mt3))
   mt <- closure(mt, 40, 4*scale)
   	mt3 <- mt[mt$closest,-which(names(mt)=="grpd")]
	mt2 <- unique(rbind(mt2, mt3))
   mt <- closure(mt, dc, scale)
   	mt3 <- mt[mt$closest,-which(names(mt)=="grpd")]
	mt2 <- unique(rbind(mt2, mt3))
   mt2 <- distill(mt2, dc)
   if (dim(mt2)[1]>0) {
	tmpv <- lapply(split(mt2, f=mt2$grp), FUN=tolerance, scale, 1)
	rule <- sapply(lapply(tmpv, FUN=">",  scale*3), any)
	if (sum(rule)>0) {
	   mt2 <- mt2[!(mt2$grp %in% names(rule)[as.vector(rule)]),]
	}
	mt  <- mt[!(mt$id %in% mt2$id),1:5]
   }
   neighbors <- function(id, V, cc) { 
	nm <- names(V)
	nm <- nm[nm!="id"]
	cpt <- V[V$id==id,nm]
	vx <- abs(V[,nm[1]] - cpt[,nm[1]]) <= cc[1]
	vy <- abs(V[,nm[2]] - cpt[,nm[2]]) <= cc[2]
	vid <- V[vy & vx,]
	freq   <- median(vid[,nm[1]])
	degree <- median(vid[,nm[2]])
	group  <- min(vid[,"id"])
	gpls  <- paste(sort(vid[,"id"]))
	return(list(len=dim(vid)[1], freq=freq, degree=degree, gpls=gpls, group=group))
   }
   nbr <- sapply(mt[,c("id")], mt, FUN=neighbors, cc=2*scale*tol)
   mt$gpls  <- nbr[5*(1:dim(mt)[1])-1]
   mt$nbr  <- as.numeric(unlist(nbr["len",]))
   mt$grp  <- unlist(nbr["group",])
   checkgroup <- function(id, V) {
	nm <- which(V$id==id)
	cpt <- unlist(V[V$id==id,"gpls"])
	V$mark <- sapply(lapply(V$gpls, "%in%", cpt),any)
	mygrp <- min(as.numeric(V[V$mark, "id"]))
	mygrp <- min(as.numeric(unlist(V[V$mark, "gpls"])), mygrp)
	return(mygrp)
   }
   mt$grp <- sapply(mt$id, mt, FUN=checkgroup)
   mt <- mt[order(mt$grp),]
	tmpv <- sapply(split(mt$grp>0, mt$grp), sum)
   mt$nbr <- rep(as.vector(tmpv), as.vector(tmpv))
   mt3 <- distill(mt[,-6], 2)
   mt2 <- rbind(mt2, mt3)
   mt2 <- mt2[order(mt2$grp),]
   mt3 <- mt2						# all meaningful groups
   mt  <- mt[!(mt$id %in% mt2$id),-6]
   mt$nbr <- 1
   mt$grp <- 0
   mtz <- rbind(mt2[,1:7], mt[,1:7])
   mtz <- mt[order(mt$grp),]
	tmpv <- tolerance(mt2, scale, tol)
   t.d2 <- tmpv[2]
   t.f2 <- tmpv[1]
   mt2 <- mt2[(mt2$nbr >= dc - mt2$sd & mt2$nbr >= 2) | mt2$nbr >= 6, ]	
	# selected wave parameters
	tmpv <- tolerance(mt2, scale, tol)
   t.d <- tmpv[2]
   t.f <- tmpv[1]
	support <- sort(unique(mt2$pair))
	mxt <- paste("Based on", length(allpairs), "pairs of directions")
	newrt.f <- sapply(split(mt2$freq, f=mt2$grp), FUN="median")
	newrt.d <- sapply(split(mt2$direction, f=mt2$grp), FUN="median")
   newrt <- unique(mt2[,c("medf","medd","nbr")])
   newrt <- newrt[!is.na(newrt$nbr),]
   names(newrt) <- c("freq","direction","count")
	row.names(newrt) <- NULL
	newrt$freq <- newrt.f
	newrt$direction <- newrt.d
   row.names(mt2) <- NULL
   mt4 <- split(mt2[, c(1,2,4)], f=mt2$grp)
   if (!D3) {
	mtz  <- mtz[mtz$freq < 0.5,]
	mt3z <- mt3[mt3$freq < 0.5,]
	if (!raw) mtz <- mt2
	mtz$df <- as.numeric(factor(mtz$pair))
	mcol <- c(gray(0.5), rainbow(dc-1))
	par(mfrow=c(1,1), cex=1)
	plot(x=mtz[,1], y=mtz[,2], xlab = "Frequency", ylim = c(-95, 95),  
		ylab = enc2utf8("Direction (\xB0)"), xlim = c(0,0.5), yaxp = c(-100, 120, 11),
		main="Supports of possible wave parameters", cex.axis = 0.85,
		col=gray(0:(dc/1.25) /(dc))[mtz$df], pch=mtz$nbr, lwd=1, cex=0.2+0.2*mtz$df)
	points(x=mt3z[,1], y=mt3z[,2], col=mcol[mt3z$nbr], pch=mt3z$nbr, 
			lwd=mt3z$nbr/2+0.5, cex=1.1)
	points(x=mt2[,1], y=mt2[,2], col=mcol[mt2$nbr], pch=mt2$nbr, 
			lwd=mt2$nbr/2+0.5, cex=1.25)
	if (dim(newrt)[1]>0) {
	   for (i in 1:dim(newrt)[1]) {
		pz <- -0.1
	   	px <- seq(from=pz,to=newrt$freq[i],by=scale[1])
	   	py <- rep(newrt$direction[i], length(px))
	   	points(x=px, y=py, col="gray41", type="l", lty=21, lwd=0.5)
	   	pz <- -110 
	   	py <- seq(from=pz,to=newrt$direction[i],by=scale[2])
	   	px <- rep(newrt$freq[i], length(py))
	   	points(x=px, y=py, col="gray41", type="l", lty=21, lwd=0.5)
	   }
	}
	nbls <- sort(unique(mt3$nbr))
	nbl2 <- sapply(nbls,min,3)
	mtext(mxt, line=0.25, col="gray21", cex=0.75)
	legend("topright",legend=nbls, bty="o", cex=0.75, 
		pch=nbls, horiz = TRUE, col=mcol[nbls], pt.lwd=nbl2)
	return(list(suggestion=newrt, detail=mt4,
			grouping_tolerance=data.frame(freq=t.f2, direction=t.d2,row.names = ""),
			tolerance=data.frame(freq=t.f, direction=t.d,row.names = ""),
			pairs=allpairs, support=support))
   }
   # ------------------------------------------- grids by group and takes maximum values
   mtz <- mt3[, 1:7]
   mtz <- mtz[mtz$freq < 0.5,]
   mt.z <- aggregate(mtz[,c("freq","direction")], by=list(mtz$grp), FUN=median)
   tmpv <- aggregate(mtz[,c("nbr","sd")], by=list(mtz$grp), FUN=max)
   mt.z$nbr <- tmpv$nbr
   mt.z$sd  <- tmpv$sd
   names(mt.z)[1] <- "grp"
   smt <- agrid(mt.z[,c("freq","direction","nbr","sd","grp")], scale=scale, math="max")
   smt$g2p <- agrid(mt.z[,c("freq","direction","sd","grp")], scale=scale, math="sum")$grp
   names(smt) <- c("freq", "direction", "count", "sd", "group", "grp2")
   smt$x <- smt$freq/scale[1]
   smt$y <- (smt$direction + 91)/scale[2]
   smt$id <- as.numeric(row.names(smt))
   prep3drsp <- function(scale, SmT) {
	   freq <- seq(from=0.0, to=0.5, by=scale[1])
	   direction <- seq(from=-90, to=90, by=scale[2])
	   asmt <- array(0, dim=c(length(freq), length(direction)))
	   asmt[as.matrix(SmT[,c("x","y")])] <- SmT$count
	   persp(x=freq, y=direction, z=asmt, theta=35, phi=20,
		main="Supports of possible wave parameters",
		zlab="count", ticktype="detailed", col="lightblue")
   }
   if (raw) { 
	prep3drsp(scale, smt[,c("freq","direction","count","x","y")])
	smt2 <- smt[smt$count > 1,]
	smt2 <- smt2[,c("freq","direction","count")]
	row.names(smt2) <- NULL
	myrt <- list(suggestion=newrt[,1:3], detail=mt4, grid=smt2, 
		grouping_tolerance=data.frame(freq=t.f2, direction=t.d2,row.names = ""),
		tolerance=data.frame(freq=t.f, direction=t.d,row.names = ""),
		pairs=allpairs, support=support)
	names(myrt)[5] <- "grids w. counts > 1"
	row.names(myrt) <- NULL
   } else {
	smt2 <- smt[smt$count >= (dc - smt$sd) & smt$count >= 2, ]
	smt2 <- smt2[,c("freq","direction","count","x","y")]
	if (dim(smt2)[1]>0) prep3drsp(scale, smt2)
	smt2 <- smt2[,c("freq","direction","count")]
	row.names(smt2) <- NULL
	myrt <- list(suggestion=newrt[,1:3], detail=mt4, grid=smt2, 
		grouping_tolerance=data.frame(freq=t.f2, direction=t.d2,row.names = ""),
		tolerance=data.frame(freq=t.f, direction=t.d,row.names = ""),
		pairs=allpairs, support=support)
   } 
   mtext(mxt, line=0.25, col="gray21", cex=0.75)
   return(myrt)
}


