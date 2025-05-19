#' Construct a rate matrix q using the output of ace 
#' This function is here because ace returns rates and an index matrix, not a matrix of rates.
#' @param aceoutput the output of ace 
#' @return A rate matrix q in the form saasi needs 
#' @export
extract_ace_q <- function(aceoutput) { 
    ind = aceoutput$index.matrix
    acerates = aceoutput$rates 
    q=ind
    rownames(q)=colnames(aceoutput$lik.anc)
    colnames(q)=rownames(q)
    for (i in 1:nrow(ind)) {
        for (j in 1:ncol(ind)) {
            q[i,j] <-  acerates[ind[i,j]]
        }
    }
    return(q)
}