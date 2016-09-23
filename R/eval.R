# -----------------------------------------------------------------------------------
#' @title	
#'	Evaluate Directional Spectrum Data for Wave Frequencies and Directions
#'
#' @description	
#'    Functions in this group are designed to estimate wave parameters based on 
#' directional periodogram records. 
#'
#'    \code{kzpdr} samples the data of wave field, and outputs the average
#' pattern of periodogram for series in a given direction. A collection of 
#' these pattern records will be sent to \code{kzpdr.eval} or \code{kzpdr.estimate} 
#' to estimate the wave frequecies and directions.   
#'
#' @param	t.D	Tolerance of direction in degree. Default is 2.
#' @param	t.F	Tolerance of frequency. Default value is 0.01.
#' @param 	...	Other arguments. 
#' \itemize{
#'  \item 	 D3	Logic. Deafult is FALSE. If TRUE, output 3D perspective plot; 
#'			otherwise, 2D plot on frequency-direction surface. 
#'  \item scale	A two element vector for grid on frequency-direction plant. 
#'			The first element is for frequency. The second is for degree
#'			of direction. Default is c(0.005,1).
#'  \item   ...	... 
#' }	
#' @inheritParams kzpdr.spikes
# 
#' @rdname 		eval
# 
#' @details
#'	    The average periodograms for a few pairs of orthogonal spatial
#'	directions can be used to identify frequencies and directions of waves.
#'
#'	    First, function \code{kzpdr} samples the spatial data and generates 
#'	periodograms in orthogonal direction pairs, and the frequencies of spikes 
#'	for each directional periodogram are identified and recorded as the output. 
#'
#'	    Then, \code{kzpdr.spikes} can be used to summarized the outputs of  
#'	\code{kzpdr}. Function \code{kzpdr.eval} or \code{kzpdr.estimate} all
#'	can be used to estimate the wave parameters (frequencies and directions).
#'	\code{kzpdr.estimate} is based on clustering-closure and the tolerances
#'	could be decided automatically. It also provides visualization of the results,
#'	thus this function is more convenient to use.
#'
#'	    Usually, if noise level is low, periodograms of a few direction pairs 
#'	may provide satisfied results. But when the noise is high, you may need to 
#'	intensively sample on different directions over the spatial data array with 
#'	\code{kzpdr}. Generally speaking, when the number of samples increases,
#'	the estimation will become more stable and reliable. 
# 
#' @return
#' 	    Both \code{kzpdr.eval} and \code{kzpdr.estimate} will return suggested 
#' 	wave frequency and direction values. The data frame of detailed estimation 
#'	for each direction are also include in their returned data list. Beside these,   
#'	\code{kzpdr.estimate} can generate 3D or 2D plots for the supports of each 
#'	suggested wave on direction-frequency parameter plane. 
# 
#' @keywords 	directional-periodogram
#' @concept 	Kolmogorov-Zurbenko periodogram
#' @concept 	directional periodogram
#' @concept 	average periodogram
#' @export
#' @seealso		\code{\link{kzpdr}}, \code{\link{kzpdr.valid}}, \code{\link{kzp2}}
#'   			\code{\link{kzpdr.tol}}, \code{\link{kzpdr.spikes}}
# 
#' @examples
# 
#'	# load pre-saved data to save running-time
#'	data(kzpdr.demo);  
#'
#'	# estimate the wave parameters
#'	kzpdr.eval(kzpdr.demo, t.D = 3, t.F = 0.01)
#'
#'	# estimation & visualization
#'	kzpdr.estimate(kzpdr.demo)
#'
#'	# For validation of the estimation, see \code{kzpdr.valid}
#'	# For reconstruction of the signals, see \code{kzrc}
# -----------------------------------------------------------------------------------

kzpdr.eval <- function(rec = ls(1),	t.D = 2, t.F = 0.01, ...) {
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
   if (length(allpairs)<3) stop("No enough data!\n\n")
	# return(list(df=cmbf[0,]))
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
	   cat("Detected direction :",
		format(paste(suggestion$grp.d,enc2utf8("\xB0"), sep=""),width=6),"\n")
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
##   Wave Parameter Estimation And Visualization Based on Clustering Closure
#
#' @inheritParams kzpdr.spikes
#' @rdname	eval  
#' @export
# -----------------------------------------------------------------------------

kzpdr.estimate <- function(rec = ls(1), ...) {
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
   if (hasArg("scale")) { scale <- dots$scale } else { scale <- c(0.005, 1) }
   cmbf$sd <- shade(cmbf$apha, cmbf$g1, cc=6) 
   cmbf <- cmbf[order(cmbf$dir, cmbf$agp.2, cmbf$sf),]
   allpairs <- sort(unique(cmbf$dir))
   if (length(allpairs)<3) stop("No enough data!!\n\n")
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


