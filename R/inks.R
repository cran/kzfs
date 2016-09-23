# -------------------------------------------------------------------
#  @title		
#  	   Internal Function For Spatial KZ Filter \code{kzmd} 
# 
#  @description	
#      \code{inks} is the control structure of the \code{kzmd} program.
#  It returns data frame of signal values and their grid positions.
#
#  @param	   ss	  Data frame of signal values and positions.
#  @param 	iterD	  Vector. Window size for KZ filter.
#  @param	iterS	  Vector. Scales of each dimension.
#  @param 	namex	  Names for each column of input data.
#  @export	   inks
#  @keywords   internal
# --------------------------------------------------------------------

inks <- function(ss, iterD, iterS, namex) 
{
   level <- length(iterD)
   l <- c(1:level)
   ss$kidx <- 1 
   for (i in level:1) {
	if (iterD[i]==1) { l <- c(i, l[1:(level-1)]); next }
	# cat(c("\n        Iteration level = ",namex[i]),"    ordered in (", namex[l], ")\n")
	ss <- ss[, c(namex[l], namex[(level+1):(level+2)]), drop=F]
	ss <- ss[do.call("order", ss[, namex[l], drop = FALSE]), , drop = FALSE]
	ss <- inkspark(ss, delta=iterD[i], sc=iterS[i])
	ss <- deltablur(ss, delta=iterD[i])
	ss <- ss[, c(namex), drop=F]
	l <- c(i, l[1:(level-1)])
   }
   return(data.frame(ss, row.names = NULL))
}

# -----------------------------------------------------------------
#   	Prepare data frame for calculation
#
#  \code{inkspark} prepares the expanded data frame for KZ filtering.
#  @rdname	inks
#  @param	delta   Integer. Window size of current dimension.
#  @param	sc	   	Integer. Scales of current dimension.
#  @export	 
# -----------------------------------------------------------------

inkspark <- function(ss, delta, sc) {

    len <- nrow(ss); w <- ncol(ss);   
    colnames(ss) -> name_ss
    ss <- ss[,c((w-1),w,(1:(w-2)))] 
    ss0 <- ss[0,]; ss1 <- ss[0,]
    # cat("\n inkspark ... ...  ",dim(ss),"\n")

    right <- floor((delta-1)/2)
    left <- ceiling((delta-1)/2)
    if (left>0) {
	# cat("\nleft: ")
 	for (j in 1:(left)) {
	    if (w<=4) ss[(1+j):len,w+1] <- abs(ss[(1+j):len,3]-ss[1:(len-j),3])
	    if (w>4)  ss[(1+j):len,w+1] <- rowSums(abs(ss[(1+j):len,3:(w-1)]-ss[1:(len-j),3:(w-1)])) 
	    ss[c(rep(0,j),ss[(j+1):len,w]-ss[1:(len-j),w])>(sc*j),w+1]<- 1
	    ss[1:j,w+1] <- 1 
	    tmp <- ss[ss[,w+1]>0,]
	    tmp[w] <- tmp[w] - (sc*j)
	    ss0 <- rbind(ss0,tmp[,1:w])
	    # cat(format(j,justify="right", width=3))
	}
	rm(tmp)
	ss0[1] <- NA; ss0[2] <- 0;
	# cat("  +",nrow(ss0)," ")
    }
    if (right>0) {
	# cat("\nright:")
      for (j in 1:(right)) {
	    if (w<=4) ss[1:(len-j),w+1] <- abs(ss[1:(len-j),3]-ss[(j+1):len,3])
	    if (w>4) ss[1:(len-j),w+1] <- rowSums(abs(ss[1:(len-j),3:(w-1)]-ss[(j+1):len,3:(w-1)])) 
	    ss[c(ss[(j+1):len,w]-ss[1:(len-j),w],rep(0,j))>(sc*j),w+1]<- 1
	    ss[(len-j+1):len,w+1] <- 1 
	    tmp <- ss[ss[,w+1]>0,]
	    tmp[w] <- tmp[w] + (sc*j)
	    ss1 <- rbind(ss1,tmp[,1:w])
	    # cat(format(j,justify="right", width=3))
     	}
	# if (right<left) # cat("   ") 
   	rm(tmp)
	ss1[1] <- NA; ss1[2] <- 0;
	# cat("  +",nrow(ss1)," ")
	ss0 <- rbind(ss0,ss1) 
	rm(ss1)
    }
    ss <- rbind(ss[1:w],ss0) 
    rm(ss0)
    ss <- ss[do.call("order", ss[, c(3:w,1), drop = FALSE]), , drop = FALSE]
    ss <- ss[c(c(3:w),c(1,2))]
     len <- nrow(ss)
     ss[2:len,w+1] <- rowSums(abs(ss[2:len,1:(w-2)]-ss[1:(len-1),1:(w-2)]),na.rm = F)
     ss[1,w+1] <- 1
     ss <- ss[ss[,w+1]>0 | !is.na(ss[,w-1]),1:w]
    # cat("\n\ncurrent sparks dim=",dim(ss),"       end at ",date(),"\n\n")
    return(ss)
}

# -----------------------------------------------------------------
#   	Smoothing data in the enlarged data frame
#
#  \code{deltablur} calculates the filtered value of each grid point.
#  @rdname	inks
#  @export	 
# -----------------------------------------------------------------

deltablur <- function(ss, delta) {

    len <- nrow(ss); w <- ncol(ss);   
    colnames(ss) -> name_ss
    ss <- ss[,c(w-1,w,(1:(w-2)))] 
    ss[is.na(ss[1]),2] <- 0
    # cat("\n deltablur ... ... ", dim(ss),"\n")

    ss0 <- data.frame(matrix(rep(NA,(w+2)*len),len,(w+2)))
    ss  <- data.frame(cbind(ss, ss0)) 
    rm(ss0)
    # ss[,(w+w+1)] <- NA
    ss[,(w+w+1)] <- ss[,2]
    ss[,(w+w+2)] <- ss[,1]
    ss[is.na(ss[w+w+2]),w+w+1] <- 0
 
    left  <- floor((delta-1)/2)
    right <- ceiling((delta-1)/2)
    if (left>0) {
	# cat("\nleft: ")
 	for (j in 1:left) {
	    ss[(j+1):len,(w+1):(w+w)] <- ss[1:(len-j),1:w]
	    ss[1:j,(w+1):(w+w)] <- ss[1:j,1:w]
	    ss[1:j,(w+1)] <- NA
	    ss[1:j,(w+2)] <- 0
	    ss[rowSums(abs(ss[3:(w-1)]-ss[(w+3):(w+w-1)]),na.rm = F)>0, w+2] <- 0
          ss[w+w+2] <- rowSums(cbind(ss[w+w+2]*ss[w+w+1], ss[w+1]*ss[w+2]),na.rm = T)		
          ss[w+w+1] <- ss[w+w+1] + ss[w+2]
          ss[w+w+2] <- ss[w+w+2]/ss[w+w+1]
          ss[ss[w+w+1]==0,w+w+2] <- NA
	    # cat(format(j,justify="right", width=3))
      }
    }
    if (right>0) {
	# cat("\nright:")
      for (j in 1:right) {
	    ss[1:(len-j),(w+1):(w+w)] <- ss[(j+1):len,1:w]
	    ss[(len-j+1):len,(w+1):(w+w)] <- ss[(len-j+1):len,1:w]
	    ss[(len-j+1):len,(w+1)] <- NA
	    ss[(len-j+1):len,(w+2)] <- 0
    	    ss[rowSums(abs(ss[3:(w-1)]-ss[(w+3):(w+w-1)]),na.rm = F)>0, w+2] <- 0
    	    ss[w+w+2] <- rowSums(cbind(ss[w+w+2]*ss[w+w+1], ss[w+1]*ss[w+2]),na.rm = T)		
    	    ss[w+w+1] <- ss[w+w+1] + ss[w+2]
    	    ss[w+w+2] <- ss[w+w+2]/ss[w+w+1]
     	    ss[ss[w+w+1]==0,w+w+2] <- NA
	    # cat(format(j,justify="right", width=3))
      }
    }

    # if (min(ss[w+w+1])==0) { # cat("\nX6 meets 0 value!!\n"); }
    ss <- ss[, c(c(3:w),(w+w+2),(w+w+1))]
    colnames(ss) <- c(name_ss)
    # cat("\n\ncurrent dim=",dim(ss),"              end at ",date(),"\n\n")
    return(ss)
}




