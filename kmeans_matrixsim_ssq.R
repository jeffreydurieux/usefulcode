
vectobin <- function(vector){
        cols <- unique(vector)
        rows <- length(vector)
        
        mat <- matrix(data = 0, nrow = rows, ncol = length(cols))
        
        for(i in 1:rows){
                mat[i, vector[i]] <- 1
        }
        
        return(mat)
}

bin <- c(rep(1,10),rep(0,10))
P <- cbind(bin, rev(bin))

lab <- P%*%1:2

C <- matrix(data = c(-2,2,2,-2), nrow = 2)
C

X <- P %*% C
X <- X + matrix(rnorm(40, sd = 5),nrow = 20)
plot(X, col = lab)

km <- kmeans(x = X, centers = 2)


Pest <- vectobin(km$cluster)
Xhat <- Pest%*% km$centers

sum( (X-Xhat)^2 )

plot(X, col = km$cluster)
