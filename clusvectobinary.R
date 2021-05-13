vectobin <- function(vector){
        cols <- unique(vector)
        rows <- length(vector)
        
        mat <- matrix(data = 0, nrow = rows, ncol = length(cols))
        
        for(i in 1:rows){
                mat[i, vector[i]] <- 1
        }
        
        return(mat)
}
