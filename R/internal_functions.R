
Read.CC.Generic <- function(filename){
    
    return(.Call("Read_Generic_R_List",filename,PACKAGE="affyio"))

}


Read.CYCHP <- function(filename){
    
    return(Read.CC.Generic(filename))
    
}
