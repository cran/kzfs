# -----------------------------------------------------------------------------------
#' @title		
#'	Multi-Scale Motion Separation with Kolmogorov-Zurbenko Periodogram Signals
#
#' @description	
#'    Motion image identification in different types of data is very important subject  
#' in many applications. Those images may depend on time and contain different scales. 
#' The simplest example is waves in the ocean coming from two different directions. 
#' One wave can be strong long scale, and another is shorter scale wave propagating 
#' in different direction. When both are covered by strong noise, data realization 
#' could be very noisy 3D structure. Similar examples can be presented in engineering,
#' acoustics, astronomy, infection diseases developments and many other fields.
#'
#'     This package is designed for the separation of motion scales in 2D motion images 
#' on different directions. To this end, KZ periodogram is utilized to identify spatial 
#' directions and frequencies of wave signals, while KZ Fourier transform provides the
#' reconstructed signals based on identified motion parameters. 
#' 
#'     By spectral analysis of original signal in different directions, we can discover 
#' main directions in which different scale waves are propagating. Intuitively, sampling
#' along the orthogonal direction of a wave will annihilate its frequency spike on the
#' corresponding periodogram. Therefore, the presence and absence of single frequency on the 
#' periodograms of different directions can be used to identify the wave direction. 
#' This method can be enriched by finding the common projected spectral spikes that detected
#' from a series of periodograms associated with different sampling directions. Identification  
#' of wave frequencies can be done symmetrically.
#'
#'     For the task of identification, this package provides functions to check averaged 
#' periodogram for data series in a given direction or a list of directions. Averaging of 
#' these directional periodograms will help to stable the variance of spectrum. Functions are
#' also provided for automatically identifying and marking prominent spectrum spikes. The closure
#' of nearest-neighbors is used to detect the clusters formed by real waves on the frequency-
#' direction plane. The algorithm is designed to resist incorrectly identified periodogram 
#' signals caused by noises, and it gives consistent estimations when the number of sampling 
#' directions increases. The accuracy of the estimations can also be improved with the increase 
#' of the sampling number.
#'
#'     In the stage of signal reconstruction, KZ Fourier transform is a powerful tool to recover 
#' signals series. Functions are designed to reconstruct 2D spatial waves from high noise background. 
#' The reconstructed signal will be averaged along the vertical lines of its propagating
#' direction. This will significantly reduce the noises. Phase shift caused by unmatched 
#' window size and frequency is remedied with special design. This enables the KZ Fourier 
#' transform be used in cases of short data series and small window sizes.
# --------------------------------------------------------------------------------------
#' @references	
#' \itemize{
#'   \item I. G. Zurbenko, The spectral Analysis of Time Series. North-Holland, 1986.
#'   \item A. G. DiRienzo, I. G. Zurbenko, Semi-adaptive nonparametric spectral estimation, 
#'		Journal of Computational and Graphical Statistics 8(1): 41-59, 1998.
#'   \item R. Neagu, I. G. Zurbenko, Algorithm for adaptively smoothing the log-periodogram, 
#'		Journal of the Franklin Institute 340: 103-123, 2003.
#'   \item I. G. Zurbenko, M. Luo, Restoration of Time-Spatial Scales in Global Temperature 
#' 		Data, American Journal of Climate Change, 1(3): 154-163, 2012. 
#'		doi:10.4236/ajcc.2012.
#' } 
#' @seealso	 \code{\link{kzpdr}}, \code{\link{kzp2}}
#' @docType package
#' @name 	kzfs
#'
#' @importFrom   kzft       kzft     kzp 	
#' @importFrom   grDevices  dev.new  gray  rainbow
#' @importFrom   graphics   abline   plot  legend  mtext  persp  par  points  text  box  image
#' @importFrom   methods    hasArg
#' @importFrom   stats      fivenum  median  na.omit  aggregate  quantile  sd  optim
#' @importFrom   utils      combn
#' @importFrom   digest     digest
#'
NULL
## NULL

# --------------------------------------------------------------------------------------
#'       The Demo Dataset For Examples of \code{kzpdr}
#'
#' A dataset containing the spectral spike records of directional periodograms
#' output by function \code{kzpdr}. It is only used for saving running-time
#' of the examples. The following code is used to generate this dataset:
#' \itemize{
#' \item \code{kzpdr.demo <- kzpdr(a, c(0, pi/4, pi/3, -pi/3, pi/18), plot=TRUE)}
#' }
#' where \code{a} is the spatial data array for a wave field.
#'
#' @format 
#' A list for a data frame and its MD5sum value.
#' The data frame \code{rec} has 18 rows and a few variables:
#' \describe{
#'   \item{direction}{angle of sampling, in degree}
#'   \item{freq}{frequency values of spikes}
#'   \item{spg}{smoothed power periodogram values of spikes}
#'   \item{...}{...}
#' }
"kzpdr.demo"

