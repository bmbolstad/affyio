


read.celfile.header <- function(filename,info=c("basic","full"),verbose=FALSE){
  compress <- FALSE

  info <- match.arg(info)

  if (info == "basic"){
    if (verbose)
      cat("Reading", filename, "to get header information")
    headdetails <- .Call("ReadHeader", filenames[[1]], compress, PACKAGE="affyio")


    return(headdetails)
  } else {
    ### full returns greater detailed information from the header. Exact details differ depending on the file format.




  }


}
  
