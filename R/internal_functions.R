
Read.CC.Generic <- function(filename, reduce.nvt=FALSE){
    
    return(.Call("Read_Generic_R_List",filename, as.integer(reduce.nvt), PACKAGE="affyio"))

}


Read.CYCHP <- function(filename){
    
    return(Read.CC.Generic(filename, reduce.nvt=TRUE))
    
}
