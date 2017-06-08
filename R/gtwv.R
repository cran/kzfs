# -------------------------------------------------------------------
#' @title		
#'	Internal Functions For \code{kzpdr}, \code{kzp2}, \code{kz.ft}, and \code{kzrc2} 
#'
#' @description	
#'   A group of internal functions used by \code{kzpdr}, \code{kzp2}, \code{kz.ft}, 
#'   and \code{kzrc2}.
#' \itemize{
#'
#'   \item \code{a2d} transfers 2D array to data frame. If input is data 
#' 	frame, it will return the original data frame.
#'
#    \item \code{adapt.m} for KZFT window size that matches the frequency.
#'
#'   \item \code{agrid} aggregates data based on given grid scale.
#'
#'   \item \code{best.cor} returns the largest correlation and related 
#'   lags on x- or y-direction for two image matrices.
#'
#'   \item \code{df2mt} transfers data frame to matrix.
#'
#'   \item \code{efg} gives projected wave frequency on given direction.
#'
#    \item \code{fit.CE} fits ellipse for a set of points on complex plane.
#'
#'   \item \code{getwave} is designed to extract data series along a given
#' 	direction in a 2D field.
#'
#'   \item \code{getwavf} extracts data series along a given direction in
#' 	a 2D field with phase arranged according to the wave frequency.
#'
#'   \item \code{markspikes} and \code{spikes.2d} are functions to mark the 
#'  	spikes of the 1D periodogram and 2D periodogram, respectively.
#'  
#    \item \code{mycat}	print the table with separation-lines between groups
#'
#' }
#'
#' @param	  df	  Data frame of signal values and positions.
#' @param      a	  2D array for position and signal values.
#' @param  angle	  Direction or vector of directions in radians.
#' @param  scale	  Vector for scale of each dimension. For example,
#'			  for a \emph{x * y} grid, the scale is \emph{c(1/x, 1/y)}.
#' @param 	math	  Function to aggregate the data. Defaults to "mean".
#'
#' @param	x.fq	  Vector of frequency values for x axis.
#' @param  y.spm	  Vector of spectrum values for y axis.
#' @param   plot	  If need to add marks on the periodogram. Defaults to TRUE.
#' @param 	 ...	  Other arguments, i.e. the "cut" threshold, etc.
#' \itemize{
#'	\item	\code{ cut : }  
#'		Set the minimum value for a marked frequency spike. Recommend to
#'		use argument \code{lvl} instead of setting this value directly.
#'	\item	\code{ lvl : }
#'		"min" or "max". Threshold strategy for marking frequency spikes.
#'		Essentially it will set the "cut" threshold value as different
#'		level. "min" is used for cases of weak singles dominating by 
#'		some strong singles. Defaults to "max". 
#' }
#' @export	   getwave
#' @keywords   internal
# ----------------------------------------------------------------------------------

# -----------------------------------------------------------------
#   	Transfer array to data frame 
#
#  Internal use only. If input is data frame, it will do nothing.
#' @rdname	getwave
#  @return	Data frame of signal values and grid positions.
#' @export	 
# -----------------------------------------------------------------

a2d <- function(a) {
  if (is.array(a)) {
      a <- na.omit(data.frame(expand.grid(x=1:dim(a)[1], 
		y=1:dim(a)[2]), obs=as.vector(a)))
  	return(a) 
  } else if (is.data.frame(a)) {
	return(a)
  } else {
	return(data.frame(a))
  }
}

# --------------------------------------------------------------------------
#    Get Adaptive Window Size for KZFT
#
#  @param	m	   Window size, used as reference. 
#  @param	N	   Length of data series.
#  @param	k	   Integer. Iteration times, as in \code{kz.ft}.
#  @param	f	   Vector. Selected frequency, as in \code{kz.ft}.
#  @param	tol	   Tolerance level for \emph{2*m2*f} difference from
# 			   the closest interger number. The default strategy, 
# 			   tol = 0, pefers minimum tolerance. For tol = 0.5*f,
# 			   the tolerance will be less than half data step,
# 			   and the result will be more close to m2.
#  @rdname	getwave
#  @export
# -------------------------------------------------------------------------

adapt.m <- function(f, m, N, k, tol=0, na.rm=FALSE) {
   H <- (floor((N-1)/k) + 1):1
   H <- H[ H > max(1/f)]
   if (length(H)==0) { if (na.rm) return(NA) else return(m) }
   nfc <- function(f){ 
	dh1 <- abs(2*H*f - floor(2*H*f)) 
	dh2 <- abs(2*H*f - ceiling(2*H*f))
	dh  <- data.frame(dh1, dh2)
	dh$rt <- ifelse(dh$dh1<=dh$dh2, dh$dh1, dh$dh2)
	return(dh$rt)
   }
   tol.m <- sapply(f, FUN= nfc)
   if (length(tol.m)==1) {
	rule <- sum(tol.m)
   	tol.M <- tol.m
   } else {
   	rule <- apply(tol.m, MARGIN=1, function(x){sum(x==0)})
   	tol.M <- rowSums(tol.m)
   }
   nm <- which(rule==max(rule))
   h <- which(tol.M==min(tol.M[nm]))
   g <- which(abs(H[nm]-m)==min(abs(H[nm]-m)))
   if (length(g)==0) { 
	if (length(h)==0) return(m)
	return(H[h[1]])
   }
   if (length(h)==0) { 
	return(H[g[1]])
   }
   if	(min(tol.M[g]) > max(1000*tol.M[h])) { g <- h }
   return(H[g[1]])
}


# -----------------------------------------------------------------
#    Aggregate data based on given grid points
#
#  v2: 2013-Nov-13
#' @rdname	getwave
#' @export
# -----------------------------------------------------------------

agrid <- function(df, scale, math="mean")
{	
	if (!is.data.frame(df)) df <- data.frame(df)
	vars <- length(scale)
	if (vars>1) {
	   maxx <- apply(df[,1:vars], 2, max, na.rm = TRUE)
	   minx <- apply(df[,1:vars], 2, min, na.rm = TRUE)
	} else {
	   maxx <- max(df[,1]); minx <- min(df[,1])
	}
	maxg <- ceiling(maxx/scale)
	ming <- floor(minx/scale)
	gs <- maxg - ming + 1
	z1 <- data.frame(df[,1:vars])
	for (i in 1:vars) { 
	   z1[,i] <- (as.integer(as.character(round(df[,i]/scale[i])))) 
	}
	if (vars>1)  { 
	    ls.z1 <- as.list(z1[,1:vars]) 
	} else {
	    ls.z1 <- list(z1[,1:vars]) 
	}
	z2 <- aggregate( df[,-c(1:vars)], by = ls.z1, FUN=math[1], na.rm=T)
	if (length(math)>1) {
	   for (i in 2:(length(math))) { 
		z2[,(vars+i)] <- aggregate( df[,-c(1:vars)], 
		by = ls.z1, FUN=math[i], na.rm=T)[,(vars+1)]
	   }
	}
	if (vars==1)  xs <- z2[,1:vars]*scale  
	if (vars >1)  xs <- sweep(z2[,1:vars],2,scale,"*")
	zs <- data.frame(cbind(xs,z2[,-c(1:vars)]))
      return(zs);
}


# -------------------------------------------------------------------------
#    Find the best correlation coefficient around "zero"
#
#  Suppose we have two image matrices a1 and a2, their are very similar 
#  but there are some phase shift between them. Therefore their largest 
#  correlation is not cor(a1, a2). {best.cor} searches for the largest
#  correlation coefficient and their spatial lags on x- or y-direction.
#
#  v2.1: 2017-Feb-22	comment out "break"s
#' @param	rc		Array. Reconstructed signal.
#' @param	sig		Array. Original signal (without noise).
#' @rdname	getwave
#' @export
# -----------------------------------------------------------------

best.cor <- function(rc, sig, ...){
	dots <- list(...)
	if (hasArg("f"))  { f <- dots$f  } else { f <- 0.5 }
	if (hasArg("q"))  { q <- dots$q  } else { q <- 0   }
	if (hasArg("level")) { 
	   level <- dots$level 
	} else { 
	   level <- max(ceiling(1/(2*f))) 
	}
	m.x <- ceiling(dim(rc)/2)[1]
	m.y <- ceiling(dim(rc)/2)[2]
	mycor <- rep(0, 2*(m.x+m.y)+6)
	x.cor <- rep(0, 2*(m.x+m.y)+6)
	y.cor <- rep(0, 2*(m.x+m.y)+6)
	mycor[1] <- stats::cor(as.vector(sig[(1:dim(rc)[1]),(1:dim(rc)[2])]), 
			 as.vector(rc), use="pairwise.complete.obs")
	if (is.na(mycor[1])) {
	   rt <- list(max.cor=NA, cor.ls=NA, shift=c(0, 0),
			x=1, y=1, ip=1, level=level, m.x=m.x, m.y=m.y)
	   return(rt)
	}
	j <- 2
	for (ix in (1:m.x)) {
	   j <- j + 1
	   x.cor[j] <- min(max(ix,0),dim(sig)[1])
	   if (dim(rc)[1]<=x.cor[j]) { x.cor[j]=0 ; break }
	   mycor[j] <- stats::cor(as.vector(sig[1:(dim(rc)[1]-x.cor[j]),1:dim(rc)[2]]), 
				 as.vector(rc[(1+x.cor[j]):dim(rc)[1],1:dim(rc)[2]]), 
				 use="pairwise.complete.obs")
	}
	j <- j + 1; j01 <- j
	for (ix in (1:m.x)) {
	   j <- j + 1
	   x.cor[j] <- min(max(ix,0),dim(sig)[1])
	   if (dim(rc)[1]<=x.cor[j]) { x.cor[j]=0 ; break }
	   mycor[j] <- stats::cor(as.vector(sig[(1+x.cor[j]):dim(rc)[1],1:dim(rc)[2]]), 
				 as.vector(rc[1:(dim(rc)[1]-x.cor[j]),1:dim(rc)[2]]), 
				 use="pairwise.complete.obs")
	   x.cor[j] <- -x.cor[j]
	}
	j <- j + 2; j02 <- j
	for (iy in (1:m.y)) {
	   j <- j + 1
	   y.cor[j] <- min(max(iy,0),dim(sig)[2])
	   if (dim(rc)[2]<=y.cor[j]) { y.cor[j]=0 ; break }
	   mycor[j] <- stats::cor(as.vector(sig[1:dim(rc)[1],1:(dim(rc)[2]-y.cor[j])]), 
				 as.vector(rc[1:dim(rc)[1],(1+y.cor[j]):dim(rc)[2]]), 
				 use="pairwise.complete.obs")
	}
	j <- j + 1; j03 <- j
	for (iy in (1:m.y)) {
	   j <- j + 1
	   y.cor[j] <- min(max(iy,0),dim(sig)[2])
	   if (dim(rc)[2]<=y.cor[j]) { y.cor[j]=0 ; break }
	   mycor[j] <- stats::cor(as.vector(sig[1:dim(rc)[1],(1+y.cor[j]):dim(rc)[2]]), 
				 as.vector(rc[1:dim(rc)[1], 1:(dim(rc)[2]-y.cor[j])]), 
				 use="pairwise.complete.obs")
	   
	   y.cor[j] <- -y.cor[j]
	}
	mx <- max(mycor[1:j], na.rm=T)
	ip <- which(mycor==mx)[1]
	mx0 <- max(mycor[c(1,3,j01+(1:3),j02+(1:3),j03+(1:3))], na.rm=T)
	ip0 <- which(mycor==mx0)[1]
	if (mx0 < (mycor[1]+0.1)) { ip0 <- 1; mx0 <- mycor[1] }
	if (mx  < (mx0 + 0.1)) { ip <- ip0; mx <- mx0 }
	rt <- list(max.cor=mx, cor.ls=mycor[1:j], 
			shift=c(x.cor[ip], y.cor[ip]),
			x=x.cor[ip], y=y.cor[ip], ip=ip, m.x=m.x, m.y=m.y)
	if (!hasArg("level")) return(rt)
	mycor <- c(mycor, rep(0,(level*2+1)^2))
	j <- j + 1
	f.x <- abs(round((1/f)*cos(q)))
	f.y <- abs(round((1/f)*sin(q)))
	for (ix in (-level:level)) for (iy in (-level:level)) {
	   j <- j + 1
	   x.cor[j] <- min(max(f.x+ix,1),dim(sig)[1])
	   y.cor[j] <- min(max(f.y+iy,1),dim(sig)[2])
	   mycor[j] <- stats::cor(as.vector(sig[(1:dim(rc)[1])+max(f.x+ix,1),
				 (1:dim(rc)[2])+max(f.y+iy,1)]),
				  as.vector(rc), use="pairwise.complete.obs")
	}
	mx <- max(mycor[1:j],na.rm=T)
	ip <- which(mycor==mx)[1]
	if (mx < (mycor[1]+0.1)) { ip <- 1; mx <- mycor[1] }
	rt <- list(max.cor=mx, cor.ls=mycor[1:j], 
			shift=c(x.cor[ip], y.cor[ip]),
			x=x.cor[ip], y=y.cor[ip], ip=ip, level=level, m.x=m.x, m.y=m.y)
	return(rt)
}


# -----------------------------------------------------------------
#   	Transfer data frame to matrix
#
#' @rdname	getwave
#  @return	Matrix
#' @export	 
# -----------------------------------------------------------------

df2mt <- function(df, scale) 
{     
	vars <- length(scale)
	mi <- apply(df[,1:vars],2,min)
	ma <- apply(df[,1:vars],2,max)
	y <- sweep(sweep(df[,1:vars],2,mi,"-"),2,abs(scale),"/")+1
	y <- as.matrix(y)
	ln <- abs((ma - mi)/scale) + 1
	x <- array(NaN,ln)
	x[y] <- df[,(vars+1)]
   	r <- list(vars); for (v in 1:vars) { r[[v]] <- mi[v]:ma[v] }
	dimnames(x) <- r
	return(x)
}


# -----------------------------------------------------------------
#   	Projected frequency (spatial frequecy on x- or y-direction)
#
#  Internal use only. 
#' @param	f	Wave frequency.
#' @param	d	Wave direction.
#' @param	ag	Sampling direction.
#  @return	Numerical. Projected frequency.
#' @rdname	getwave
#' @export
# -----------------------------------------------------------------

efg <- function(f, d, ag){ 
   delta <- d - ag
   ag <- ifelse(ag >  45,  90-ag, ag)
   ag <- ifelse(ag < -45, -90-ag, ag)
   rt <- abs(cos(delta*pi/180)*(eval(f)/cos(ag*pi/180)))
   rt <- ifelse(rt>0.5, 1-rt, rt)
   return(rt)
}


# -----------------------------------------------------------------
# 	Extract data series on given direction
#
#' @param	cp	  Sampling scheme. Values 0, 2, and others : normal;
#'			  1: interpolation; 3: errors
#' @export	 
#' @rdname	getwave
# ------------------------------------------------------------------
 
getwave <- function(df, angle, cp=0) {
   dx <- max(df[,1]);
   dy <- max(df[,2]);
   ceta <- angle%%pi;
   if (ceta>(pi/2)) df <- df[order(df$x,-df$y),]
   atca <- abs(tan(ceta))
   if (round(1/atca) >= dy) {	
	df$e <- df$y
	if (cp==1) cp <- 0
   } else if (round(atca) >= dx) {
	df$e <- df$x
	if (cp==1) cp <- 0
   } else if (atca > 1) {
	ry <- round(dy/2) + 1
	df$e  <- df$x - (df$y - ry)/tan(ceta)
  } else {
	rx <- round(dx/2) + 1
	df$e <- df$y - tan(ceta)*(df$x - rx)
  }
  df$er <- round(df$e)- df$e 
  if (cp==1) {
	df$e1 <- floor(df$e)
	df$e2 <- ceiling(df$e)
	df$w1 <- 1-abs(abs(df$e)-abs(df$e1))
	df$w2 <- abs(abs(df$e)-abs(df$e1))
	if (atca > 1) { okxy <- "y" } else { okxy <- "x" }
	wv1 <- df[,c("obs",okxy,"w1","e1")]
	names(wv1)[c(1,4)] <- c("obs1","e")
	wv2 <- df[,c("obs",okxy,"w2","e2")]
	names(wv2)[c(1,4)] <- c("obs2","e")
	v12 <- merge(wv1, wv2, by=c("e", okxy), all = T)
	v12$obs <- v12$w1*v12$obs1 - v12$w2*v12$obs2
	v12[is.na(v12$obs2) & v12$w1>0.5,]$obs <- 
		v12[is.na(v12$obs2) & v12$w1>0.5,]$obs1
	v12[is.na(v12$obs1) & v12$w2>0.5,]$obs <- 
		v12[is.na(v12$obs1) & v12$w2>0.5,]$obs2
	v12 <- v12[!is.na(v12$obs),]
	if (atca > 1) { 
	   v12 <- v12[order(v12$e, v12$y),]
	} else { 
	   v12 <- v12[order(v12$e, v12$x),]
	}
	rt <- v12[,c(okxy, "e", "obs")]
	return(rt)
  } else {
	df$e	<-  round(df$e)
	df$o <- df$er
	if (length(df[df$er == -0.5, ]$e) > 1) {
	   df[df$er == -0.5, ]$e  <- df[df$er == -0.5, ]$e - 1
	   df[df$er == -0.5, ]$er <- 0.5
	}
	if (atca > 1) { 
	   df <- df[order(df$e, df$y),]
	} else { 
	   df <- df[order(df$e, df$x),]
	}
	if (cp==3) df$obs <- df$er
  }
  return(df)
}


# -----------------------------------------------------------------
#      Fit Ellipse for A Set of Points on the Complex Plane
#
#  \emph{a} and \emph{b} are the axes of the ellpise. Both could be
#  the long axis or the short axis. \emph{ceta} is the angle between
#  the x axis and the \emph{a} axis. The ellipse central is on 0+0i.
#
#  @param	a0	   The initial value of \emph{a}.
#  @param	b0	   The initial value of \emph{b}. 
#  @param	ceta	   Angle in radian. 
#  @rdname	getwave
#  @export
# -----------------------------------------------------------------

fit.CE <- function(xy, a0 = max(Mod(xy)), b0 = min(Mod(xy)), 
   ceta = Arg(xy[which(Mod(xy)==a0)]), ...){
   CE <- function(abc){
	p1 <- abc[2]*(Re(xy)*cos(abc[3]) + Im(xy)*sin(abc[3]))
	p2 <- abc[1]*(Im(xy)*cos(abc[3]) - Re(xy)*sin(abc[3]))
	sum((abc[1]*abc[2] - sqrt(p1^2 + p2^2))^2)
   }
   stats::optim(c(a0, b0, ceta), CE, ...)
}


# ------------------------------------------------------------------
# 	Sampling scheme for reconstruction
#' @param	  f1	  	Wave frequency
#' @param	  angle	Wave direction
#' @param	  rlvl	Coefficient to control the averaging level.
#' @rdname	getwave
# ------------------------------------------------------------------

getwavf <- function(df, angle, f1, rlvl=1) {
   dx <- max(df[,1]);
   dy <- max(df[,2]);
   atca <- abs(tan(angle))
   emark <- function(x, y, hy, hx=1, r=0) { 
	round(rlvl*(x*hx + sign(cos(angle))*0.4999 - y*hy + r*hy))
   }
   pseek <- function(hy, hx) {
	p = 1; q = 1; e=1; r0=0
	while (p==q & e==p) {
	   e <- emark(x=0, y=0, hy=hy, r=r0, hx=hx)
	   p <- emark(x=1, y=0, hy=hy, r=r0, hx=hx)
	   q <- emark(x=2, y=0, hy=hy, r=r0, hx=hx)
	   r0 <- r0 - 1
	   if (abs(r0) > 1000) { break; browser() }
	}
	r0 <- r0 + 2
	return(r0)
   }
   if (round(1/atca) >= dy) {	
	df$e <- df$y 
   } else if (round(atca) >= dx) {
	df$e <- df$x 
   } else {
	hy <- f1*cos(angle)
	hx <- f1*sin(angle)
	if (atca > 1) { 
	   r0 <- pseek(hy, hx)
	   df$e <- emark(x=df$x, y=df$y, hx=hx, hy=hy, r=r0)
	} else {
	   r0 <- pseek(hx, hy)
	   df$e <- emark(x=df$y, y=df$x, hx=hy, hy=hx, r=r0)
	}
   }
   return(df)
}



# -------------------------------------------------------------------------------
#     Print the result table with separation-lines between groups
#
#  @param	x	   Data frame for result table.
#  @param	mk	   Vector. The marked values for different groups.
#  @rdname	getwave
#  @export
# -------------------------------------------------------------------------------

mycat <- function(x, mk) {
   myline <- " ----------------------------------------------------------------"
   myline <- paste(myline, "-------------------------------------------", sep="")
   y <- diff(as.numeric(mk))
   y <- c(0, y)
   for (i in 1:dim(x)[1]){
	if (i == 1) { 
	   w = rep(1,dim(x)[2])
	   for (j in 1:dim(x)[2]){ 
		w[j] <- max(nchar(names(x)[j])+1, max(nchar(x[,j]),na.rm = T))+1
		cat(format(names(x)[j], width=w[j], 
		zero.print = "", justify="right")," ") 
	   }
	   cat("\n")
	   cat(substr(myline,1,4+sum(w)+dim(x)[2]*1)," \n")
	}
	if (abs(y[i])>0) cat(substr(myline,1,4+sum(w)+dim(x)[2]*1)," \n")
	for (j in 1:dim(x)[2]){ 
	   cat(format(x[i,j], width=w[j], 
		zero.print = "", justify="right")," ") 
	}
	cat("\n")
   }
}


# ----------------------------------------------------------------
#     Function to locate the spikes of the 2D periodogram
#
#' @rdname	getwave
#' @export
# ----------------------------------------------------------------
spikes.2d <- function(x.fq, y.spm, nm=10) {
 	rg <- c(0,max(x.fq))
	len <- length(x.fq)
	step <- min(abs(x.fq[1:(len-1)]-x.fq[2:len]))
	cut <- mean(y.spm)+2*sd(y.spm)
	vat <- x.fq[y.spm>cut]
	ss <- unique(sort(vat))
	slen <- length(which(diff(ss)>1.001*step))
	while (slen>nm) {
	   cut <- cut + 0.5*sd(y.spm[y.spm>cut])
	   vat <- x.fq[y.spm>cut]
	   ss <- unique(sort(vat))
	   slen <- length(which(diff(ss)>1.001*step))
	}
	idx <- c(1:len)[y.spm>cut]
	cut0 <- min(cut,y.spm[1])
	end <- min(c(1:len)[y.spm<=cut0])-1
 	if (end>0) { 
		exd <- c(1:end)
		vat <- vat[-exd]
		idx <- idx[-exd]
		y <- y.spm[-exd]
	} else if (end==0) { 
		y <- y.spm 
	}
	if (length(vat)==0) { return(vat) }
	Ruler <- function(y.spm, cut) {
	   vat <- x.fq[na.omit(y.spm) >= cut]
	   if (length(vat)==0) { return(vat) }
	   aty <- y.spm[y.spm >= cut];
	   delta <- abs(diff(vat))
	   eq2 <- c((round(delta,10) == round(step,10)),FALSE)
	   seg <- 1; eq1 <- eq2+0;
	   for (j in 1:length(vat)) { 
		eq1[j] <- seg
		if (!eq2[j]) { seg <- seg + 1 }
	   }
	   aty <- as.vector(sapply(split(aty, f=eq1),FUN=max))
	   return(aty)
	}
	bigy <- y[y >= cut]
	mark.all <- as.vector(max(y.spm))
	for (cut0 in bigy[order(bigy)]) {
	   mark <- Ruler(y.spm, cut0)
	   mark.all <- unique(c(mark, mark.all))
	}
	idx <- c(1:len)[which(y.spm %in% mark.all)]
	vat <- x.fq[idx]
	aty <- y.spm[idx]
	exd <- c(0) 
	l <- dim(y.spm)[1]
	idx.win <- c(1,l,l+1,l-1,l+2,l-2,2*l,2*l+1,2*l-1,2*l+2,2*l-2,2)
	# idx.win <- c(1,l,l+1,l-1)
	while (length(exd)>0) {
		delta <- abs(vat[1:(length(vat)-1)]-vat[2:length(vat)])
		eq2  <- c((round(delta,10) <= round(step,10)),FALSE)
		idlt <- diff(idx); 
		eid  <- c(idlt %in% idx.win, FALSE)
		cmp2 <- c(((aty[1:(length(aty)-1)]-aty[2:length(aty)])<=0),FALSE)
		exd <-  c(1:(length(delta)+1))[(eq2 & cmp2 & eid)]
		exd <- exd[!is.na(exd)]
		if (length(exd)>0) { 
			vat <- vat[-exd]; idx <- idx[-exd]; aty <- aty[-exd]
		}
		delta <- abs(vat[1:(length(vat)-1)]-vat[2:length(vat)])
		eq1  <- c(FALSE,(round(delta,10) <= round(step,10)))
		idlt <- diff(idx); 
		eid  <- c(FALSE, idlt %in% idx.win)
		exd <-  c(1:(length(delta)+1))[eq1 & eid]
		exd <- exd[!is.na(exd)]
		if (length(exd)>0) {
			vat <- vat[-exd]; idx <- idx[-exd]; aty <- aty[-exd]
		}
	}
	return(list(freq=vat[order(aty,decreasing = TRUE)],
			power=aty[order(aty,decreasing = TRUE)]))
}

# ----------------------------------------------------------------------
#     Function to mark the spikes of 1D periodogram
#
#' @rdname 	getwave
#' @export	markspikes
# ----------------------------------------------------------------------

markspikes <- function(x.fq, y.spm, plot=TRUE, ...) {
	dots <- list(...)
	if (hasArg("off")) { off <- dots$off } else { off <- 0 }
	if (hasArg("rg"))  {  rg <- dots$rg  } else {  rg <- c(0,max(x.fq)) }
	if (hasArg("cut")) { 
	   cut <- dots$cut 
	} else { 
	   cut0 <- max(y.spm, na.rm=TRUE) + 1
	   myc <- 1
	   while (TRUE) {
		y <- y.spm[ y.spm < cut0 ]
	   	cut <- mean(y, na.rm=TRUE) + 1.65*sd(y, na.rm=TRUE)
		myc <- myc + 1
		# abline(h=cut,col=myc, lty=21)
		if (is.na(cut)) cut <- cut0
		if (cut0 == cut) break
		cut0 <- cut
		if (myc >= 2) break
	   }
	}
	if (hasArg("lvl"))  { lvl <- dots$lvl } else { lvl <- "max" }
	if (hasArg("n"))    { n <- dots$n } else { n <- 1 }
	cut0 <- quantile(y.spm, probs = 0.925, na.rm = T)
	# abline(h=cut0, col="red", lty=2, lwd=0.1)
	if (lvl == "max") { 
		cut <- max(cut0, cut)
	} else if (lvl == "mean") {
		cut <- mean(cut0, cut)
	} else {
		cut <- min(cut0, cut)*0.95 + 0.05*max(cut0, cut)
	}
	len <- length(x.fq)
	y.max <- max(y.spm, na.rm=TRUE)*0.99
	y.min <- min(y.spm, na.rm=TRUE)
	# abline(h=cut,col="grey2", lty=21, lwd=0.1)
	step <- min(abs(x.fq[1:(len-1)]-x.fq[2:len]))
	vat <- x.fq[na.omit(y.spm) > cut]
	if (length(vat)==0) { 
	   x.fq <- x.fq[-c((len-5):len)]
	   y.spm <- y.spm[-c((len-5):len)]
	   len <- length(x.fq)
	   cut <- y.max
	   vat <- x.fq[na.omit(y.spm) > cut] 
	}
	cut0 <- min(cut,y.spm[1],na.rm=TRUE)
	y.spm[is.na(y.spm)] <- 0
	end <- min(c(1:len)[y.spm[1:len]<=cut0],na.rm=TRUE)-1
 	if (end>0) { exd <- c(1:end); vat <- vat[-exd];y <- y.spm[-exd] }
      if (end==0){ y <- y.spm } 
	Ruler <- function(y.spm, cut) {
	   vat <- x.fq[na.omit(y.spm) >= cut]
	   if (length(vat)==0) { return(vat) }
	   aty <- y.spm[y.spm >= cut];
	   delta <- abs(diff(vat))
	   eq2 <- c((round(delta,10) == round(step,10)),FALSE)
	   seg <- 1; eq1 <- eq2+0;
	   for (j in 1:length(vat)) { 
		eq1[j] <- seg
		if (!eq2[j]) { seg <- seg + 1 }
	   }
	   aty <- as.vector(sapply(split(aty, f=eq1),FUN=max))
	   return(aty)
	}
	bigy <- y[y >= cut]
	mark.all <- as.vector(max(y.spm))
	for (cut0 in bigy[order(bigy)]) {
	   mark <- Ruler(y.spm, cut0)
	   mark.all <- unique(c(mark, mark.all))
	}
	idx <- c(1:len)[which(y.spm %in% mark.all)]
	vat <- x.fq[idx]
	aty <- y.spm[idx]
	yrs <- round(vat,6)
	step <- 0.5*n*step
	if (length(vat)==0) { return(vat) }
	   delta <- abs(diff(vat))
	   eq2 <- c((round(delta,10) <= round(step,10)),FALSE)
	   cmp2 <- c(diff(aty)>=0,FALSE)
	   exd <- c(1:length(aty))[(eq2 & cmp2)]
	   exd <- exd[!is.na(exd)]
	if (length(exd)>0) {vat <- vat[-exd]; aty <- aty[-exd]; yrs <- yrs[-exd]}
	   delta <- abs(vat[1:(length(vat)-1)]-vat[2:length(vat)])
	   eq1 <- c(FALSE,(round(delta,10) <= round(step,10)))
	   exd <- c(1:length(aty))[eq1]
	   exd <- exd[!is.na(exd)]
	if (length(exd)>0) {vat <- vat[-exd]; aty <- aty[-exd]; yrs <- yrs[-exd]}
	   base <- max((y.max*0.02),0.25)
	   hat <- aty + base*rep(c(2,8),1+length(aty)/2)[1:length(aty)]
	   hat <- hat + (y.max-hat)/5
	   hat[hat>=y.max] <- y.max - (y.max-y.min)*0.009*
	   	sample.int(100,length(hat[hat>=y.max]),replace = TRUE)
	   idx <- vat>=(rg[1]-step) & vat< (rg[2]-(off+1)*step)
	   if (any(idx) & plot) {
		abline(v = vat[idx], lty = 21, lwd = .1, col = "grey40")
		text(vat[idx],hat[idx],yrs[idx],col=2)
	   }
	return(vat[idx])
}

