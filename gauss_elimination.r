# Desc.: Solves a linear system with the
# Gauss elmination method.
# Group: 2
# Obs.: It assumes that the determinant
# is not zero.
# 
# I: A matrix of the coeficients(A),
# a column matrix(b),
# A column matrix with the solution
# value(X)
# 
# Ex.: n <- 3
# A <- matrix(c(1,2,3,2,1,-1,1,-1,-1), n, n)
# b <- matrix(c(3,0,-2), n, 1)
# X <- matrix(c(0,1,1), n, 1)

gauss_elimination <- function(A, b){
    n <- dim(b)[1]
    # initiates it here because of scope
    X <- matrix(rep(1,n),n,1)
    expanded <- cbind(A, b)
    
    # makes the matrix to be in 
    # echelon form(triangular form)
    for (i in 1:n-1) {
	pivot <- A[i,i]
	# Couldn't find a way to get 
	# only the multipliers we want
	mult <- A[,i]/rep(n,pivot)	
	for(j in i+1:n){
	    expanded[,j] <- expanded[,j] - mult[j]*expanded[,i]
	}
    }
    
    # The main solution
    for(i in n:1) {
	sum <- 0
	for(j in n:1){
	    sum <- sum + expanded[i,j]*X[j,1]
	}

	# takes the last column of the expanded(b)
	# subtracts de sum
	X[i] <- expanded[i,n+1] - sum
    }   

    return (X)
}
