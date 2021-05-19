
Read.CC.Generic <- function(filename, reduced.output=FALSE){
    
    return(.Call("Read_Generic_R_List",filename, as.integer(reduced.output), PACKAGE="affyio"))

}


Read.CYCHP <- function(filename){
    
    return(Read.CC.Generic(filename, reduced.output=TRUE))
    
}
