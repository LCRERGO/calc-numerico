
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
	    iter = iter + 1
	}
	X[k] <- (b[k] - sum)/A[k,k]
    }

    end_time = Sys.time()
    timer = end_time - start_time

    return (X)
}

# Name: lu_factoration.r
# Desc.: Solves a linear system with the
# LU factoration method.
# Group: 2
# Obs.: It assumes that the determinant
# is not zero.
# 
# I: A matrix of the coeficients(U),
# a column matrix(b),
# O: A column matrix with the solution
# value(X)
# 
# Ex.: n <- 3
# U <- matrix(c(1,2,3,2,1,-1,1,-1,-1), n, n)
# b <- matrix(c(3,0,-2), n, 1)
#
# X <- matrix(c(0,1,1), n, 1)
lu_factoration <- function(U,b){
    n <- dim(b)[1]
    L <- matrix(0,n,n)
    y <- matrix(0,n,1)
    x <- matrix(0,n,1)

    start_time <- Sys.time()
    iter <- 0

    for(i in 1:n){ # Columns
	for(j in 1:n){ # Lines

	    iter <- iter + 1
	    if(j > i ){
		m <- U[j,i]/U[i,i]
		L[j,i] <- round(m,7)
		U[j,] <- round(U[j,] - m * U[i,],7) 
	    }
	    if(j == i)
		L[i,i]=1
	}
    }

    U[j,] <- round(U[j,] - (U[j,j-1]/U[j-1,j-1]) * U[j-1],7)
    iter <- iter + 1


    y[1] <- b[1,1]/L[1,1]
    iter <- iter + 1
    for (i in 2:n){
	sum <- 0
	for (j in 1:(i-1)){
	    iter <- iter + 1
	    sum <- sum + (L[i,j] * y[j])
	}
	y[i] <- (b[i,1] - sum)/L[i,i]
    }

    x[n] <- y[n]/U[n,n]
    iter <- iter + 1
    for (k in (n-1):1){
	sum <- 0
	for (j in (k+1):n){
	    sum <- sum + U[k,j]*x[j]
	    iter <- iter + 1
	}
	x[k] <- (y[k] - sum) / U[k,k]
    }
    end_time <- Sys.time()

    timer <- end_time - start_time

    result <- c(x,timer,iter)
    return(result)
}
