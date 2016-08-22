# --------------------------------------------------------------------------
#' @title		
#'		Internal Function Used In Wave Parameter Estimation
#'
#' @description
#'    Function \code{commcore} is designed to find the common elements 
#' in a given list, which is the pool of all potential parameters. The
#' common elements should appear in most data groups of paired directional
#' periodograms.
#'
#' 	\code{kzpdr.rec} will search all variables in R global environment  
#' for available records of directional periodograms, and return the 
#' specified part as output. 
#'
#'	\code{kzpdr.pairs} returns a list of all available pairs of sampling 
#' directions and the counts number for the spikes on each direction.
#'
#'	\code{kzpdr.prj} is used for calculating the projected parameters
#' based on directional periodograms.
#'
#'	Function \code{ckGV} refers to a given record set of directional
#' periodogram and checks duplicated directions. The returned value is a 
#' list of sampling directions without those duplications.
#'
#'	\code{shade} decides if a spectral signals would sometime disappear 
#' in the given sampling directions. The return value is a vector of 0/1. 
#' 0 means the signal should be seen in all given directions. 0 means not. 
#' 2 and larger may happen when there are very closed sampling angles.
#'
#'	\code{tight} is used to get the minimum difference in a group of
#' potential directions or frequencies. 
#'
#'	\code{tolerance} is used to calculate the tolerance of the estimations.
#'
#'	Function \code{closure} returns the closure of the nearest neighbors. 
#'
#'	\code{distill} excludes the duplicate supports of a cluster. Usually,
#' duplications are supports for the same parameters but from the same pairs,
#' grouping under the given tolerance condition.
#'
#' @param	pool		List for the potential parameters.
#' @param	STEP		The tolerance interval for a data element.
#'				Defaults to 1.
#' @param	DIFF		The allowed absents times for a suggested 
#'				common element. Default value is 1.
#' @param   digit		Significant digits. Defaults to 0.
#' @rdname	kzprj
#' @export
#' @keywords   internal
# ---------------------------------------------------------------------------

commcore <- function(pool, STEP=1, DIFF=1, digit=0) {
   u <- as.vector(unlist(pool))
   flg <- names(unlist(pool))
   flg <- gsub(")[0-9]+",")",flg)
   a <- round(u, digit)
   o <- order(u)
   v <- data.frame(v=a,u=u,d=DIFF,flg=flg)
   v <- v[o,]
   while (TRUE) {
	v$upp <- c(999,diff(v[,1]))
	v$dwn <- abs(c(rev(diff(rev(v[,1]))),999))
	v$x <- ifelse((v$upp == 0 | v$dwn == 0), 1, 0)
	v[v$upp > STEP,]$upp <- 0
	v[v$dwn > STEP,]$dwn <- 0
	v$delta <- round(v$dwn - v$upp,digit)
	v$cc <- c(0, v[1:(dim(v)[1]-1), "delta"])
	v$cc <- round(sign(v$cc)*(abs(v$delta)),digit)
	v$cc <- ifelse(v$delta<0,v$cc,0)
	v$grp <- v[,1] + v$delta + v$cc
	v$v <- v$grp
	v <- v[order(v$v),]
	if (all((v$delta + v$cc)==0)) break
   }
   v$one <- 1
   t <- aggregate(v[,c("one","u")], by=list(v$flg,v$grp), FUN=sum)
   names(t)[1:2] <- c("flg", "grp")
   t$u <- t$u / t$one
   t$one <- 1; 
   w <- aggregate(t[,c("one","u")], by=list(t$grp), FUN=median)
   w$nm <- round(w$u/w$one,digit)
   w2 <- aggregate(t[,c("one","u")], by=list(t$grp), FUN=sum)
   DIF <- aggregate(v$d, by=list(v$grp), FUN=max)
   condition.1 <- (w2$one >= length(pool)-max(DIF[,2])) | (w2$one > 6)
   condition.2 <- (w2$one >= length(pool)-DIF[,2]) | (w2$one > 5)
   sgstn  <- unique(w[which(condition.1),"nm"])
   sgstn2 <- unique(w[which(condition.2),"nm"])
   cc <- STEP; tol <- STEP
   itv <- setdiff(sgstn,sgstn2)
   if (length(itv)>0) {
   	igv <- w[w$nm %in% itv,1]
   	tmp <- abs(sapply(igv, FUN="-", v$v))
      tmp[round(tmp,8)==0] <- NA
   	cc <- max(apply(tmp, FUN="min", MARGIN=2, na.rm=TRUE))
	flag <- matrix(rep(cc,2),ncol=2)
	for (i in 1:length(cc)) {
	   flag[i,1] <- which(tmp[,i]==cc[i], arr.ind =TRUE)[1]
	   flag[i,2] <- i
	}
	cc <- ifelse(is.na(tmp[flag]),cc,tmp[flag])
	if (cc >= 3*STEP) { tol <- cc ; cc <- STEP ; sgstn <- sgstn2 }
   }
   nearest <- function(x, Y, cc=1) { 
	na.omit(ifelse(abs(x - Y) > cc, NA, Y))[1]
   }
   if (length(sgstn)>0) {
	v$ok <- sapply(v$u, FUN=nearest, sgstn, cc)
	a[o] <- unlist(v$ok)
	er <- aggregate(na.omit(abs(v$u - v$ok)), 
		by=list(v[!is.na(v$ok),]$ok), FUN=max)
	er <- round(er$x, digits=nchar(v$ok[1])-1)
   } else {
	er <- NULL
   }	
   return(list(suggest=sgstn, grp=unlist(a), tolerance=tol, error=er))
}


# -----------------------------------------------------------------------------
#    Calculating of the projected parameters based on directional periodogram
#
#' @rdname kzprj
#' @export
# -----------------------------------------------------------------------------

kzpdr.proj <- function(rec = ls(1)) {
   rec <- kzpdr.rec(rec)[,-c(3,4)]
   fok <- kzpdr.pairs(rec)
   ag1 <- fok$ag1
   ag2 <- fok$ag2
	s1 <- rec[round(rec$ang*180/pi,2) %in% ag1,]
	s2 <- rec[round(rec$ang*180/pi,2) %in% ag2,]
	s1$ag1 <- round(s1$ang*180/pi,2)
	s1$ag2 <- s1$ag1 + 90
	s2$ag2 <- round(s2$ang*180/pi,2)
	s2$ag1 <- s2$ag2 - 90
	s1$direction <- paste("(",s1$ag1,", ",s1$ag2,")",sep="")
	s2$direction <- paste("(",s2$ag1,", ",s2$ag2,")",sep="")
	names(s1) <- c("dir","f1","agp.1","g1","g2")
	names(s2) <- c("dir","f2","agp.2","g2","g1")
	s1$fp <- sapply(s1$f1/cos(s1$agp.1),min, 100)
	s2$fp <- sapply(s2$f2/cos(s2$agp.2),min, 100)
	s3 <- s1[s1$f1>0.29 & s1$fp>0.5 & s1$fp<=1,]
	s4 <- s2[s2$f2>0.29 & s2$fp>0.5 & s2$fp<=1,]
	s3 <- s3[s3$g2>=45 & s3$g2<=90,]
	s4 <- s4[s4$g2>=45 & s4$g2<=90,]
	s3$f1 <- 1 - s3$f1
	s4$f2 <- 1 - s4$f2
	s1 <- rbind(s1, s3)
	s2 <- rbind(s2, s4)
	cmb  <- merge(s1[,1:5], s2[,1:5])
	cmb$f0 <- sqrt(cmb$f1^2 + cmb$f2^2)
	cmb$b1 <- cmb$g2>45 & cmb$g2<=90
	cmb$sf <- ifelse(cmb$b1, sin(cmb$agp.2), -sin(cmb$agp.1))*cmb$f0
	cmb$beta <- atan(cmb$f2/cmb$f1) # + cmb$agp.1
	cmb$apha <- atan(cmb$f1/cmb$f2) # + cmb$agp.2 
	cmb$apha -> cmb$agp.2
	cmb$beta -> cmb$agp.1 
	cmb$apha <- round(cmb$apha*180/pi, 0)
	cmb$beta <- round(cmb$beta*180/pi, 0)
	cmb <- cmb[,-c(4,6,8,9)]
	s1 <- cmb[,c("dir", "g1", "g2", "sf", "agp.2", "apha")]
	s2 <- cmb[,c("dir", "g1", "g2", "sf", "agp.1", "beta")]
	names(s2) <- c("dir", "g1", "g2", "sf", "agp.2", "apha")
	s1$formula <- "f2/f1"
	s2$formula <- "f1/f2"
	cmbf <- rbind(s1, s2)
 	cmbf$b2 <- ifelse(cmbf$formula=="f1/f2", cmbf$g1, cmbf$g2)
	cmbf$apha  <- cmbf$b2 + cmbf$apha
	cmbf$agp.2 <- (cmbf$b2*pi/180) + cmbf$agp.2
	 s3 <- cmbf
	 s3$apha  <- (s3$b2 - s3$apha)
	 s3$agp.2 <- (s3$b2*pi/180 - s3$agp.2)
	 cmbf <- rbind(cmbf, s3)
	cmbf$apha  <- ifelse(cmbf$apha > 90, cmbf$apha-180, cmbf$apha)
	cmbf$apha  <- ifelse(cmbf$apha < -90, (cmbf$apha)%%180, cmbf$apha)
	cmbf$agp.2 <- ifelse(cmbf$agp.2 > pi/2, cmbf$agp.2-pi, cmbf$agp.2)
	cmbf$agp.2 <- ifelse(cmbf$agp.2 < -pi/2, (cmbf$agp.2)%%pi, cmbf$agp.2)
	cmbf <- cmbf[order(cmbf$dir,cmbf$sf,cmbf$agp.2),]
	cmbf <- unique(cmbf[,1:6])
	cmbf$sd <- shade(cmbf$apha, cmbf$g1)
   return(cmbf)
}


# -----------------------------------------------------------------------------
# 	Find available records of directional periodograms in R environment
#
#' @param  rec  Data frame or list of the outputs from function \code{kzpdr}. 
#'		    It includes the marked spike frequencies and sampling directions.
#'		    Default is for searching from all variables in the R environment. 
#' @param  md5  MD5 code for the specified data array. The returned data will be
#'		    the records for the array with the same MD5 code. The default is
#'		    to return the first record variable in alphabet order.
#' @rdname kzprj
#' @export
# -----------------------------------------------------------------------------

kzpdr.rec <- function(rec = ls(1), md5="") {
   if (missing(rec)) {
	nm <- ls(1)
   } else if (is.character(rec) & is.vector(rec)) {
	nm <- rec
   } 
   if (exists("nm")) {
	oknm <- c("direction","freq","spg","max")
	nc <- rep(FALSE,length(nm))
	for (i in (1:length(nm))) {
	   rec <- get0(nm[i])
	   nc[i] <- is.list(rec) & !is.data.frame(rec)
	   nc[i] <- nc[i] & ("MD5" %in% names(rec))
	   nc[i] <- nc[i] & ("rec" %in% names(rec))
	}
	for (i in which(nc)) {
	   rec <- get0(nm[i])
	   nc[i] <- nc[i] & all(oknm %in% names(rec$rec))
	   if (nchar(md5)==32) {
	      nc[i] <- nc[i] & (rec$MD5 == md5)
	   }
	}
	if (length(which(nc)) > 1) {
	   cat("\nMultiple kzpdr objects. Use the first one: ")
	   cat(nm[which(nc)[1]],"\n\n")
	} else if (length(nc)==0) {
	   return(NULL)
	}
	nc <- which(nc)
	rec <- get0(nm[nc[1]])
	md5 <- as.vector(rec$MD5)
	if (length(md5)>0) cat("\n",md5,"\n\n")
   }
   if (length(rec)==0) return(rec)
   if ("MD5" %in% names(rec)) rec <- rec$rec
   rec$ang <- as.numeric(gsub(enc2utf8("\xB0"), "", rec$direction))*pi/180
   rec$ang <- ifelse(rec$ang<=(-pi/2), rec$ang%%pi, rec$ang)
   rec[,1] <- paste(round(rec$ang*180/pi,2), enc2utf8("\xB0"), sep="")
   rec <- unique(rec)
   return(rec)
}


# -----------------------------------------------------------------------------
#    Get all the pairs of directional periodograms
#    and their counts for the spikes on each direction
# 
#  Returned value is a list of counts and direction pairs.
#
#' @rdname kzprj
#' @export
# -----------------------------------------------------------------------------

kzpdr.pairs <- function(rec = ls(1)) {
   rec <- kzpdr.rec(rec)
   fc <- aggregate(rec$freq, by=list(rec$direction),FUN=length)
   names(fc) <- c("direction","count")
   fc <- fc[order(as.numeric(gsub(enc2utf8("\xB0"), "", fc[,1]))),]
   fg <- split(fc[,1], fc[,2])
   rc <- as.numeric(gsub(enc2utf8("\xB0"), "", as.vector(fc[,1])))
   rc -> tmpv
   rc <- rev(rc[rc > 0])
   rc <- rc[(rc - 90) %in% tmpv]		# evaluate in pairs
   ag2 <- rc; ag1 <- ag2 - 90;
   return(list(fc=fc, fg=fg, rc=rc, ag1=ag1, ag2=ag2))
}


# -----------------------------------------------------------------
#		Check and delete duplicated sampling directions
# 
#' @param 	angle	 	Vector of sampling directions in radian.
#' @rdname	kzprj
#' @export
# -----------------------------------------------------------------

ckGV <- function(angle, rec) {
   if (length(rec)>0) {
	relc <- rec
	aglc <- unique(relc$direction)
	okag <- paste(round((180/pi)*angle,2),enc2utf8("\xB0"),sep="")
	okat <- which(okag %in% aglc)
	if (length(okat)>0) angle <- angle[-okat]
   }
   return(angle)
}


# ------------------------------------------------------------------------------
#   Check if spectral signals would sometimes disappear in given directions
# 
#' @param	alpha		Vector of checked directions in degree.
#' @param	   gs    	Vector of directions for all available periodograms.
#' @param	   cc	   	Tolerance of wave parameters. Defaults to 3.
#' @rdname	kzprj
#' @export
# -----------------------------------------------------------------------------

shade <- function(alpha, gs, cc=3){
   neighbors <- function(x, Y) { 
	if (sum(abs(Y-x)<= cc) == 0) return(0) 
	length(unique(Y[abs(Y-x)<= cc]))
   }
   if (length(alpha)==0) return(NA)
   df <- data.frame(alpha=alpha, g1=gs, g2=90+gs)
   df$c1 <- sapply(df$alpha, FUN=neighbors, df$g1)
   df$c2 <- sapply(df$alpha, FUN=neighbors, df$g2)
   df$c1 <- df$c1 + df$c2
   df$c2 <- (abs(df$alpha - df$g1)<=cc) + (abs(df$alpha - df$g2)<=cc)
   df$d1 <- df$c1 - df$c2
   return(df$d1)
}

# ------------------------------------------------------------------------------
#  		Get the minimum gap of potential directions or frequencies
#
#  Used to find the tolerance of the estimations. 
#  The argument \code{alpha} and \code{gs} could be frequencies for this function.
#' @rdname	kzprj
#' @export
# -----------------------------------------------------------------------------

tight <- function(alpha, gs, cc=3){
   sv <- split(alpha, gs)
   sv <- lapply(sv,sort)
   svd <- lapply(sv,diff)
   ok <- unlist(svd)
   ok <- ok[round(ok,8)!=0]
   ok <- ok[ok<cc*2]
   if (length(ok)>0) { 
	mstp <- max(ok, na.rm = TRUE)
      return(mstp*1.1)
   } else {
	svd <- diff(sort(alpha))
	ok <- svd[round(svd,8)!=0]
	ok <- ok[ok<cc*2]
	if (length(ok)==0) return(ok)
	Mstp <- max(ok, na.rm = TRUE)
	return(Mstp*1.1)
   }
}

# ------------------------------------------------------------------------------
#  Calculate the tolerance of the estimations. 
#' @param	 tm2		Data frame. Internal data table for potential parameters.
#' @param  scale		The scale of gridding data. 
#' @rdname	kzprj
#' @export
# -----------------------------------------------------------------------------

tolerance <- function(tm2, scale, tol=1) {
   rg.f <- split(tm2$freq, f=tm2$grp)
	if (length(rg.f)>0) {
	   t.f <- max(sapply(rg.f, FUN="max") - sapply(rg.f, FUN="min"))
	} else {
	   t.f <- scale[1]*tol
	}
   rg.d <- split(tm2$direction, f=tm2$grp)
	if (length(rg.d)>0) {
	   t.d <- max(sapply(rg.d, FUN="max") - sapply(rg.d, FUN="min"))
	} else {
	   t.d <- scale[2]*tol
	}
   return(c(t.f, t.d))
}

# ------------------------------------------------------------------------------
#  Find the closure of nearest neighbors. 
#' @param	mdc		Integer. Specified approximate size of the closure.
#' @param	 sc		Numeric. Specified tolerance level of the cluster.
#' @rdname	kzprj
#' @export
# -----------------------------------------------------------------------------

closure <- function(tm2, mdc, sc) {
   distance <- function(id, V, mdc) { 
	nm <- names(V)
	nm <- nm[nm!="id"]
	cpt <- V[V$id==id,nm]
	sd <- cpt$sd
	vx <- abs(V[,nm[1]] - cpt[,nm[1]])/0.5
	vy <- abs(V[,nm[2]] - cpt[,nm[2]])/180
	vz <- sqrt(vx^2 + vy^2)
	vid <- V[order(vz)[1:max((mdc-sd),2)],"id"]
	return(list(mdist=max(vz[vid]), mid=paste(sort(vid))))
   }
   nearest.neighbors <- function(id, V, mdc) {
	nm <- names(V)
	nm <- nm[nm!="id"]
	vid <- as.numeric(unlist(V[V$id==id,"grpd"]))
	v <- V[V$id %in% vid, ]
	sd <- v[v$id==id, "sd"]
	vx <- (paste(v$grpd) == paste(v[v$id==id,"grpd"]))
	len  <- sum(vx)
	if (sum(vx) >= (mdc - sd)) {
	   group  <- max(vid)
	   freq   <- median(v$freq)
	   degree <- median(v$direction)
	} else {
	   group  <- 0
	   freq   <- NA
	   degree <- NA
	}
	return(list(len=len, freq=freq, degree=degree, group=group))
   }
   dist <- sapply(tm2[,c("id")], tm2, FUN=distance, mdc)
   tm2$grpd <- dist[2*(1:dim(tm2)[1])]
   nbr <- sapply(tm2$id, tm2, FUN=nearest.neighbors, mdc)
   tm2$nbr  <- unlist(nbr["len",])
   tm2$grp  <- unlist(nbr["group",])
   tm22 <- tm2[tm2$grp>0,]
   rg.f <- split(tm22$freq, f=tm22$grp)
   t.f <- lapply(lapply(rg.f, FUN="range"),FUN="diff")
   rg.d <- split(tm22$direction, f=tm22$grp)
   t.d <- lapply(lapply(rg.d, FUN="range"), FUN="diff")
   rule <- (t.f < 2*sc[1]) & (t.d < 2*sc[2])
   tm2$closure <- tm2$grp > 0
   tm2$closest <- tm2$grp %in% names(rule[rule])
   return(tm2)
}

# ------------------------------------------------------------------------------
#  			Exclude the duplicate supports in a cluster
#
#  Duplication is support for one potential parameter but from the same pair of
#  directional periodogram, grouping under the given tolerance condition.
#' @rdname	kzprj
#' @export
# -----------------------------------------------------------------------------

distill <- function(tm2, mdc) {
   tm2$nbr <- as.numeric(tm2$nbr)
   tm2$grp <- as.numeric(tm2$grp)
   tm2 <- tm2[order(tm2$nbr, tm2$grp, decreasing = TRUE),]
	tmpv <- split(rownames(tm2),f=tm2$id)
	tmpv <- as.vector(sapply(tmpv, FUN=function(m) m[1]))
   tm2 <- tm2[rownames(tm2) %in% tmpv, 1:7]
   tm2 <- tm2[order(tm2$grp),]
	tmpv <- split(tm2$pair,f=tm2$grp)
	tmpw <- sapply(sapply(tmpv,unique),FUN=function(m) !is.na(m))
	tmpw <- as.vector(sapply(tmpw, sum))
	tmpv <- as.vector(sapply(sapply(tmpv,FUN=function(m) !is.na(m)), sum))
   tm2$nbr <- rep(tmpw, tmpv)
   tm2 <- tm2[tm2$nbr >= mdc - tm2$sd, ]
   if (dim(tm2)[1]==0) return(tm2)
   tm2 <- tm2[order(tm2$grp, tm2$pair),]
	tmpv <- split(tm2$pair, f=paste(tm2$pair, tm2$grp))
	tmpv <- sapply(tmpv, FUN=function(m) !is.na(m))
	tmpv <- lapply(lapply(tmpv, sum), ">", 1)
	tmpv <- names(tmpv[unlist(tmpv)])
   tm2$tw <- paste(tm2$pair, tm2$grp) %in% tmpv
   tm2 <- tm2[!(tm2$tw & tm2$nbr==1), ]
	tmpv <- split(tm2[!tm2$tw,]$freq, f=tm2[!tm2$tw,]$grp)
	tmpw <- as.vector(sapply(tmpv, median))
	tmpv <- split(tm2$freq, f=tm2$grp)
	tmpv <- as.vector(sapply(lapply(tmpv, '>', 0),sum))
   tm2$medf <- rep(tmpw, tmpv)
	tmpw <- split(tm2[!tm2$tw,]$direction, f=tm2[!tm2$tw,]$grp)
	tmpw <- as.vector(sapply(tmpw, median))
   tm2$medd <- rep(tmpw, tmpv)
	tm2$dif.f <- abs(tm2$freq - tm2$medf)
	tm2$dif.d <- abs(tm2$direction - tm2$medd)
   tm2 <- tm2[order(tm2$grp, tm2$pair, tm2$dif.f, tm2$dif.d),]
	tm2$g1 <- substring(gsub(", [0-9]+)","",tm2$pair),2)
	tm2$g2 <- gsub("[(+-]+[0-9]+, ","",tm2$pair)
	tm2$g2 <- gsub("[)]","",tm2$g2)
	tm2$g1 <- as.numeric(tm2$g1)
	tm2$g2 <- as.numeric(tm2$g2)
	tm2$dif.gf <- c(999, diff(tm2$grp))
	tm2$dif.g1 <- c(999, diff(tm2$g1))
	tm2$dif.g2 <- c(999, diff(tm2$g2))
   tm2$c1 <- tm2$dif.gf==0 & tm2$dif.g1==0 & tm2$dif.g2==0
   tm2 <- tm2[!(tm2$c1),1:10]
   tm2 <- tm2[,-8]
   return(tm2)
}

