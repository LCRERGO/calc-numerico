# Name: gauss_elimination.r
# Desc.: Solves a linear system with the
# Gauss elmination method.
# Group: 2
# Obs.: It assumes that the determinant
# is not zero.
# 
# I: A matrix of the coeficients(A),
# a column matrix(b),
# O: A column matrix with the solution
# value(X)
# 
# Ex.: n <- 3
# A <- matrix(c(1,2,3,2,1,-1,1,-1,-1), n, n)
# b <- matrix(c(3,0,-2), n, 1)
#
# X <- matrix(c(0,1,1), n, 1)

gauss_elimination <- function(A, b){
    start_time = Sys.time()
    iter = 0
    n <- dim(b)[1]
    # initiates it here because of scope
    X <- matrix(rep(0,n),n,1)
    
    # makes the matrix to be in 
    # echelon form(triangular form)
    for (k in 1:(n-1)) {
	for (i in (k+1):n){
	    iter <- iter + 1
	    mult <- A[i,k]/A[k,k]
	    A[i,k] = 0
	    for(j in (k+1):n){
		A[i,j] <- (A[i,j] - (mult*A[k,j]))
	    }
	    b[i] <- b[i] - mult*b[k]
	}
    }
    
    # The main solution
    X[n] <- b[n]/A[n,n]
    for (k in (n-1):1){
	sum <- 0
	for(j in (k+1):n){
	    sum <- sum + A[k,j]*X[j]
	    iter = iter + 1;
	}
	X[k] <- (b[k] - sum)/A[k,k]
    }

    end_time = Sys.time();
    timer = end_time - start_time;

    return (X)
}
