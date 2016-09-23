# -------------------------------------------------------------------
#' @title		
#'	Internal Function For Directional KZ Periodogram \code{kzpdr}
#'
#' @description	
#'   A group of internal functions used by \code{kzpdr} function.
#' \itemize{
#'   \item \code{getaway} is designed to extract data series along 
#' a given line in a 2D field.
#'
#'   \item \code{agrid} is to aggregate data based on given grid points.
#'
#'   \item \code{a2d} transfers 2D array to data frame. If input is data frame, 
#' it will return the original data frame.
#'
#'   \item \code{df2mt} transfers data frame to matrix.
#
#	 \item \code{mycat}	print the table with separation-lines between groups
#  
#    \item \code{markspikes} and \code{spikes.2d} are functions to mark the 
#  spikes of the 1D periodogram and 2D periodogram, respectively.
#' }
#'
#' @param	  df	  Data frame of signal values and positions.
#' @param      a	  2D array for position and signal values.
#' @param  angle	  Direction or vector of directions in radians.
#' @param  scale	  Vector for scale of each dimension. For example,
#'			  for a \emph{x * y} grid, the scale is \emph{c(1/x, 1/y)}.
#' @param 	math	  Function to aggregate the data. Defaults to "mean".
#'
#  @param	x.fq	  Vector of frequency values for x axis.
#  @param  y.spm	  Vector of spectrum values for y axis.
#  @param   plot	  If need to add marks on the periodogram. Defaults to TRUE.
#  @param 	 ...	  Other arguments, i.e. the "cut" threshold, etc.
#  \itemize{
# 	\item	\code{ cut : }  
# 			  Set the minimum value for a marked frequency spike. Recommend to
# 			  use argument \code{lvl} instead of setting this value directly.
# 	\item	\code{ lvl : }
# 			  "min" or "max". Threshold strategy for marking frequency spikes.
# 			  Essentially it will set the "cut" threshold value as different
# 			  level. "min" is used for cases of weak singles dominating by 
# 			  some strong singles. Defaults to "max". 
#  }
#' @export	   getwave
#' @keywords   internal
# ----------------------------------------------------------------------------------

getwave <- function(df, angle) {
   dx <- max(df[,1]);
   dy <- max(df[,2]);
   ceta <- angle%%pi;
   if (ceta>(pi/2)) df <- df[order(df$x,-df$y),]
   atca <- abs(tan(ceta))
   if (round(1/atca) >= dy) {	
	df$e <- df$y 
   } else if (round(atca) >= dx) {
	df$e <- df$x 
   } else if (atca > 1) {
	ry <- round(dy/2)+1; 
	df$e <- round(df$x - (df$y - ry)/tan(ceta))
   } else {
	rx <- round(dx/2)+1; 
	df$e <- round(df$y - tan(ceta)*(df$x - rx))
   }
   return(df)
}


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
	y <- sweep(sweep(df[,1:vars],2,mi,"-"),2,scale,"/")+1
	y <- as.matrix(y)
	ln <- ((ma - mi)/scale) + 1
	x <- array(NaN,ln)
	x[y] <- df[,(vars+1)]
   	r <- list(vars); for (v in 1:vars) { r[[v]] <- mi[v]:ma[v] }
	dimnames(x) <- r
	return(x)
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



# -----------------------------------------------------------------
#    Get Adaptive Window Size for KZFT
#
#  @param	m	   Window size, used as reference. 
#  @param	N	   Length of data series.
#  @param	k	   Integer. Iteration times, as in \code{kz.ft}.
#  @param	f	   Vector. Selected frequency, as in \code{kz.ft}.
#  @param	tol	   Tolerance level for \emph{2*m2*f} difference from
#			   the closest interger number. The default strategy, 
#			   tol = 0, pefers minimum tolerance. For tol = 0.5*f,
#			   the tolerance will be less than half data step,
#			   and the result will be more close to m2.
#  @rdname	getwave
#  @export
# -----------------------------------------------------------------

adapt.m <- function(f, m, N, k, tol=0) {
   H <- (floor((N-1)/k) + 1):1
   H <- H[ H > max(1/f)]
   nfc <- function(f){ 
	dh1 <- abs(2*H*f - floor(2*H*f)) 
	dh2 <- abs(2*H*f - ceiling(2*H*f))
	dh  <- data.frame(dh1, dh2)
	dh$rt <- ifelse(dh$dh1<=dh$dh2, dh$dh1, dh$dh2)
	return(dh$rt)
   }
   tol.m <- sapply(f, FUN= nfc)
   rule <- apply(tol.m, MARGIN=1, function(x){sum(x==0)})
   nm <- which(rule==max(rule))
   tol.M <- rowSums(tol.m)
   h <- which(tol.M==min(tol.M[nm]))
   g <- which(abs(H[nm]-m)==min(abs(H[nm]-m)))
   if (length(g)==0) { 
	if (length(h)==0) return(m)
	return(H[h[1]])
   }
   if (length(h)==0) { 
	return(H[g[1]])
   }
   if	(tol.M[g] > 1000*tol.M[h]) { g <- h }
   return(H[g[1]])
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
#  @rdname	getwave
#  @export
# ----------------------------------------------------------------
spikes.2d <- function(x.fq, y.spm) {
 	rg <- c(0,max(x.fq))
	cut <- mean(y.spm)+2*sd(y.spm) 
	len <- length(x.fq)
	step <- min(abs(x.fq[1:(len-1)]-x.fq[2:len]))
	vat <- x.fq[y.spm>cut]
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
	   aty <- y[y >= cut];
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
#  @rdname 	getwave
#  @export	markspikes
# ----------------------------------------------------------------------

markspikes <- function(x.fq, y.spm, plot=TRUE, ...) {
	dots <- list(...)
	if (hasArg("off")) { off <- dots$off } else { off <- 0 }
	if (hasArg("rg"))  { rg <- dots$rg } else { rg <- c(0,max(x.fq)) }
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
	step <- n*min(abs(x.fq[1:(len-1)]-x.fq[2:len]))
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
	   aty <- y[y >= cut];
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
	   idx <- vat>=(rg[1]-step) & vat< (rg[2]-(off+3)*step)
	   if (any(idx) & plot) {
		abline(v = vat[idx], lty = 21, lwd = .1, col = "grey40")
		text(vat[idx],hat[idx],yrs[idx],col=2, cex=0.75)
	   }
	return(vat[idx])
}

