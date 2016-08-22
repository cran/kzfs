# -----------------------------------------------------------------------------------
#' @title		
#'			Yet Another Multi-dimensional Kolmogorov-Zurbenko Filter 
#'
#' @description	
#'    This implement of spatial KZ-filter works for any dimensions. It is designed 
#' for cases with sparse data in large time-space. 
#' 
#' @param     ss 		Data frame with value column after time/space columns.
#' @param     window	Vector for window size of each dimension.
#' @param	  scale	Vector for scale of each dimension.
#' @param     k		Iteration times of KZ filter. Defaults to 1.
#' @param     edges	Logic. Defaults to TRUE. FLASE means clear the data that  
#'				are located outside the time-space range of input data.
#' @return		Data framework with the similar structure as the input data.
#  @details		See introduction of KZ filters in \code{kza::kz}.
#' @concept       Kolmogorov-Zurbenko filter
#' @keywords 	KZ-filter
#' @export
#' @seealso		\code{\link[kza]{kz}}
#' @examples
#'	zs <- rbind(c(0,5,1,40),c(12,6,1,10),c(6,7,1,20),c(15,15,4,80))
#'	colnames(zs) <- c("x","y","z","v")
#'	zs <- kzmd(data.frame(zs), scale=c(1,1,1), window=c(3,5,3), k=4)
#'	u <- zs[zs$z==1, -3]
#'	x = sort(unique(u$x))
#'	y = sort(unique(u$y))
#'	z=df2mt(u, scale=c(1,1))	# Transfer from data frame to matrix.
#'	image(x=x, y=y, z=z)
#'	
# ----------------------------------------------------------------------------------


kzmd <- function(ss, window, scale, k = 1, edges=TRUE)
{	
	Tstart <- Sys.time()
	idx = FALSE
	if (!is.data.frame(ss)) ss <- data.frame(ss)
	vars <- length(window)
	x <- ss[1:vars]
	y <- ss[1+vars]

      cat(paste("Checking the data ...            ", date()),"\n")

	if (length(dim(x))==0) { dim(x) <- c(length(x),1); colnames(x) <- "x" }    
	if (length(dim(y))==0) { dim(y) <- c(length(y),1); colnames(y) <- "y" }    
	if (length(dim(window))==0) { dim(window) <- c(length(window),1) }
	if (length(dim(scale))==0) { dim(scale) <- c(length(scale),1) }
	xname <- colnames(x); yname <- colnames(ss)[1+vars];
	maxx <- apply(x, 2, max, na.rm = TRUE)
	minx <- apply(x, 2, min, na.rm = TRUE)
	cat("   Data size : ",paste(c(nrow(x)," rows   ",(vars+ncol(y))," cols  ")), " \n")
	cat("         max = ",format(paste(signif(maxx,digits=4)),justify="right", width=10), " \n")
	cat("         min = ",format(paste(signif(minx,digits=4)),justify="right", width=10), " \n")
	cat("      window = ",format(paste(window),justify="right", width=10),"\n")
	cat("       scale = ",format(paste(scale), justify="right", width=10),"\n")

	if ( nrow(unique(x)) < nrow(x)) 
		stop("\n","Warning: there are duplicate records in this dataset!","\n")
	s <- matrix(0, nrow = nrow(x)-1, ncol = ncol(x))
	if(nrow(x) != nrow(y))
		stop(paste("The lengths of 'x' and 'y' must be equal",nrow(x),length(y)))
	if(length(window) != length(scale))
		stop(paste("The lengths of 'scale' and 'window' must be equal"))
	if(ncol(ss) <= length(scale))
		stop(paste("There should be a column for your data!"))

	if (!idx) {
		s <- as.vector(rep(as.integer(1),nrow(ss)))
		ss <- cbind(ss,s)
		colnames(ss)[ncol(ss)] <- "kidx"
		ss$kidx <- 1
	}
	rm(s, x, y)

      cat(paste("Ready to go ...                 ", date()),"\n\n")

	ss <- na.omit(ss)			# initialize the data

	Wdth <- as.integer(round(window/scale))
	bias <- c(((Wdth-1)%%2)*scale*k/2, c(0,0))
	# print(bias)
	if (any(Wdth==0)) {stop("\nWindow parameter should be greater or equal to scale parameter")}
      namess <- colnames(ss)
	
	m  <- as.integer(((window/scale)-1)/2)
	if (vars>1) {
	   ss[,1:vars] <- round(sweep(sweep(ss[,1:vars],2,minx,"-"),2,scale,"/"))
	   ss[,1:vars] <- sweep(ss[,1:vars],2,k*m+1,"+")
	} else {
	   ss[,1] <- round((ss[,1]-minx)/scale)
	   ss[,1] <- ss[,1]+(k*m+1)
	}
	for (j in 1:k) {
		cat(paste("\n k = ",j," ...           ", date()),"\n")
		ss <- inks(ss, iterD=Wdth, iterS=rep(1,vars), namex=namess)
	} 
	ss[,1:vars] <- sweep(ss[,1:vars],2, k*m+1, "-")
	ss[,1:vars] <- sweep(ss[,1:vars],2, scale, "*")
	ss[,1:vars] <- sweep(ss[,1:vars],2, minx, "+")
	if (!idx) ss <- ss[,1:(vars+1)]
      if (!edges) {  
	   cf <- (m-1)*k*scale*0
	   for (v in 1:vars) { 
		ss <- ss[ss[,v]<= (maxx[v]-cf[v]) & ss[,v]>= (minx[v]+cf[v]),] 
	   } 
	}
	if (any(bias>0)) ss <- data.frame(t(apply(ss,1,"+",bias)))
	cat("\nDONE! Finshed in ")
	# cat(paste(as.numeric(gsub("Time difference of", "", round(Sys.time()-Tstart,2)))))
	cat(round(Sys.time()-Tstart,2))
	cat(' s\n\n\a\a\a')

	return(ss)
}



