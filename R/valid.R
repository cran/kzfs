# -----------------------------------------------------------------------------------
#' @title		
#'	    Validate Estimated Wave Parameters Under Given Tolerance Setting
# 
#' @description	
#'     For a given tolerance setting, \code{kzpdr.valid} will provide cross-validation 
#' information for related results of wave parameter estimations.
# 
#' @rdname 		valid
#' @inheritParams kzpdr.spikes
#' @inheritParams kzpdr.eval
#' @param  level	\code{level} control the cross-validation process:
#'			integer number \code{k} means to run cross-validation by excluding 
#'			\code{k} pairs of directional samples each time. Default value is 1.
# 
#' @details
#'	    Due to the nosies or other reasons, there may exist fake spike signals in 
#'	the directional periodograms. Cross-validation will evaluate estimations by 
#'	excluding one or few measurements, and identify the data points that casued 
#'	inconsistent estimations. A table will be given to summarize the validation
#'	results.
#'
#'	    For a given tolerance setting, if the validation shows consistent results,
#' 	related estimation would be reliable.
#'
#' @keywords 	directional-periodogram
#' @concept 	Kolmogorov-Zurbenko periodogram
#' @concept 	directional periodogram
#' @concept 	average periodogram
#' @export
#' @seealso		\code{\link{kzpdr}}, \code{\link{kzpdr.eval}}, \code{\link{kzpdr.tol}}
#'
#' @examples
#'
#'	# load pre-saved data 
#'	data(kzpdr.demo)
#'
#'	# validation
#'	kzpdr.valid(kzpdr.demo, t.D = 2, t.F = 0.01, level = 1)
#
# -----------------------------------------------------------------------------------

kzpdr.valid <- function(rec = ls(1), t.D = 2, t.F = 0.01, level=1) {
   rec <- kzpdr.rec(rec)
   wvnmbr <- kzpdr.spikes(rec)
   fok <- kzpdr.pairs(rec)
   ag1 <- fok$ag1
   ag2 <- fok$ag2
   cmbf <- kzpdr.proj(rec)
   if (length(ag1)<4) stop("No enough data for cross-validation!\n\n")
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
   cat("  t.D =", format(t.D, width=4, justify="right"), 
	 "    t.F =", t.F, "\n")
   cat("  t.d =", format(t.d, width=4, justify="right"), 
	 "    t.f =", round(t.f,6), "\n\n")
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
	   unselected <- paste("(",ag1[notin[,i]],enc2utf8("\xB0, "),
				ag2[notin[,i]],enc2utf8("\xB0)"),sep="")
	   cat("\n", k, ".", i," Exclude", unselected, "\n\n")
	   ort <- kzpdr.eval(cmbf[cmbf$dir %in% selected,], t.D=t.D, t.F=t.F)
	   rof <- unique(ort$df[,c("grp.f","grp.d")])$grp.f
	   rod <- unique(ort$df[,c("grp.f","grp.d")])$grp.d
	   vnm <- max(length(rof), length(rod))
	   if (vnm>0) {
		sumrec0 <- data.frame(td=t.D, tf=t.F, rf=rof, rd = rod, nm = vnm, k=k, i=i)
		sumrec0$tlf <- round(ort$tolerance$freq, 6)
		sumrec0$tld <- round(ort$tolerance$direction, 6)
		sumrec0$order <- order(sumrec0$rf, sumrec0$rd)
	   } else {
		sumrec0 <- data.frame(td=t.D, tf=t.F, rf=NA, rd = NA, nm = vnm, k=k, i=i)
		sumrec0$tlf <- NA; sumrec0$tld <- NA; sumrec0$order <- 1
	   }	
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
	scnt <- aggregate(tmp0[,c("id","grd","grf")], 
			by=list(tmp0[,c("exclude")]), FUN=mypst)
	tmp0 <- merge(sref,scnt)
	scnt <- aggregate(tmp1[,c("id","grd","grf")], 
			by=list(tmp1[,c("exclude")]), FUN=mypst)
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
	   cat("\n\n  Results are inconsist on cross-validation level", lvl,"\n\n")
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
      cat("\n"); mycat(sref[,-7],SGNx$i); cat("\n");
	cat("  Note: t.f = tolerance of frequency, ")
	cat("t.d = tolerance of direction (in degree)\n")
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
		cat("  Results are consistent on cross-validation level", lvl, "\n\n")
	   } else {
		cat("    Error may exist in the excluded directions. Or it may\n")
		cat("    need to be re-grouped based on new tolerance setting.\n\n")
	   }
	}
   }}
}


# -----------------------------------------------------------------------------
#' @title
#'   Search Appropriate Tolerances Setting for Wave Parameter Estimation
#'
#' @description	
#'     \code{kzpdr.tol} will help to find the feasible tolerance settings for 
#' the wave parameter estimation. 
#'
#' @inheritParams kzpdr.spikes
#' @param	t.D	Vector for search range of direction tolerance.
#'			Default is 1:10 (in degree).
#' @param	t.F	Vector for search range of frequency tolerance. 
#'			Default value is c(0.01).
#'
#' @details
#'     Since the expected wave number is known(see \code{kzpdr.spikes}), we can  
#' search for the tolerance settings that would generate estimations with wave
#' number in this range. A table will be presented to summary feasible settings.
#'
#'     The searching process would stop when it finsihed the search range, or the 
#' increasing of the tolerance led to null result. 
#'
#' @rdname	tolerance
#' @export
#' @seealso		\code{\link{kzpdr}}, \code{\link{kzpdr.eval}}
#'  			\code{\link{kzpdr.valid}}, \code{\link{kzpdr.spikes}}
#'
#' @examples
#'
#'	# load pre-saved data 
#'	data(kzpdr.demo); 
#'
#'	# search for tolerance 
#'	kzpdr.tol(kzpdr.demo, t.D = c(1,2,3), t.F = 0.005)
#'
# -----------------------------------------------------------------------------

kzpdr.tol <- function(rec = ls(1), t.D = seq(1,10,1), t.F = 0.01) {
   rec <- kzpdr.rec(rec)
   wvnmbr <- kzpdr.spikes(rec)
   fok <- kzpdr.pairs(rec)
   ag1 <- fok$ag1
   ag2 <- fok$ag2
   cmbf <- kzpdr.proj(rec)
   if (length(ag1)<3) stop("No enough data!\n\n")
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
   cat("\nFeasible tolerance settings in searching range:\n\n")
   cat("\n"); mycat(sref,mark); cat("\n");
}



# -----------------------------------------------------------------------------
#' @title
#'  		Count Spikes For Available Directional Periodogram Records 
#'
#' @description	
#'    Function \code{kzpdr.spikes} summarizes available periodogram pattern 
#' records collected from the outputs of \code{kzpdr}, gives expected wave
#' number. This number is used by \code{kzpdr.tol} and \code{kzpdr.valid} 
#' in searching feasible tolerance setting and validation of estimated wave 
#' parameters.
#'
#' @details
#'    The expected wave number is defined as the mode of all the spike counts.
#'
#' 	If any of the sampling direction in the available directional periodogram 
#' records, say \emph{A}, happens to be orthogonal to a wave direction \emph{B}, 
#' then there will be no spike appear on related periodogram for the wave 
#' propagated in direction \emph{B}. Related spike counts will be less than the 
#' expected wave number. The absence of spike(s) in one direction can be taken
#' as the evidence for the existing of wave(s) in its orthogonal direction. 
#'
#' 	Usually, it is very rare to have a simpling direction orthogonal to a wave  
#' direction. But if we know an approximate wave direction, we can take more 
#' samplings around its orthogonal direction. Since KZ periodogram can separate 
#' wave spikes in very close frequencies, we may get more accurate estimation
#' for this wave direction with this method. It is also possible to use this way 
#' to validate estimations get by other approaches.
#'
#' @rdname spikes
#' @param  rec	Data list from the outputs of function \code{kzpdr}. It includes the
#'			data frame for the marked frequency values and corresponding directions.
#'			Defaults is searching for available records in the environment. 
#' @export
#' @seealso		\code{\link{kzpdr}}, \code{\link{kzpdr.eval}}
#'  			\code{\link{kzpdr.valid}}, \code{\link{kzpdr.tol}}
#' @examples
#'
#'	# load pre-saved data 
#'	data(kzpdr.demo); 
#'
#'	# count spikes
#'	kzpdr.spikes(kzpdr.demo)
#'
# -----------------------------------------------------------------------------------

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
   tb <- table(fc[,2])
   wvnmbr <- as.numeric(names(tb[which(tb == max(tb))]))
   mypst <- function(x){ y=paste(x[1]); 
		if (length(x)>1) { for (i in 2:length(x)) y=paste(y,x[i],sep=", ") }; y }
   if (length(sys.parents())>1) { cat("\n"); return(wvnmbr) }
   tb  <- as.numeric(names(fg))
   sgt <- tb[which(tb==wvnmbr)]-1
   sgt <- sgt[(sgt %in% tb) & !(sgt %in% wvnmbr)]
   if (length(sgt)>0) {
	angles <- as.numeric(gsub(enc2utf8("\xB0"), "", unique(unlist(fg[(sgt)])))) - 90
	angles <- ifelse(angles < -90 | angles > 180, angles%%180, angles)
	cat("\nSuggested wave direction: ", 
	     mypst(paste(angles,enc2utf8("\xB0"),sep="")), "\n")
	cat("\nSuggested wave numbers: ", mypst(unique(wvnmbr)), "\n\n")
	cat("Note: The spike counts for direction", mypst(unique(unlist(fg[sgt]))))
	cat(" is only", sgt, "\n")
	cat("      Other directions have", paste(names(fg)[[-sgt]],"spikes"))
	if (length(fg) > 2) cat(" or more")
   } else {
	cat("\nNo conclusion for wave directions.\n")
	cat("\nSuggested wave numbers: ", mypst(unique(wvnmbr)), "\n\n")
	cat("Note: Most directions have", wvnmbr)
	cat(" spike(s) on their periodograms. ")
   }
   sg2 <- tb[which(tb<wvnmbr)]
   sg2 <- sg2[(sg2 %in% tb) & !(sg2 %in% wvnmbr) & !(sg2 %in% sgt)]
   if (length(sg2)>0) { 
	angles <- as.numeric(gsub(enc2utf8("\xB0"), "", unique(unlist(fg[(sg2)])))) - 90
	angles <- ifelse(angles < -90 | angles > 180, angles%%180, angles)
	cat("\n      Please check if there are more than one waves on ")
	cat(mypst(paste(angles,enc2utf8("\xB0"),sep="")))
   }
   cat("\n\n")
}



