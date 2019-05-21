# Name: interpolation.r
# Desc.: A set of functions to
# get the coeficients of a polynomial
# Group: 2

# It requires a library to solve
# systems of equations
source("gauss_elimination.r")

# I: A vector of the x elements(x)
# A vector of the y elements(y)
# O: A vector with the coefitients
# of the polynomial(A)
#
# Ex.:x <- c(0,1)
# y <- c(1.35,2.94)
#
# A <- c(1.35,1.59)

interpolation_vandermonde <- function(x,y){
    n <- length(x)
    # Has to be in this form for the
    # function specification
    y <- matrix(y, length(y), 1)
    # An auxiliar matrix to be passed
    # to solve function
    aux <- matrix(rep(0,n^2),n,n)

    # iterates over the x to get
    # the matrix
    for(i in 1:n){
    	for (j in 0:(n-1)) {
	    aux[i,j+1] <- x[i]^j	    
	}
    }

    A <- gauss_elimination(aux, y)

    return (A)
}

# I: A vector of the x elements(x)
# A vector of the y elements(y)
# O: A vector with the coefitients
# of the polynomial(A)
interpolation_lagrange <- function(x, y){
    n <- length(x)
    # Starts the auxiliary variables
    numerator <- rep(1, n)
    denominator <- rep(1, n)
    
    # iterates over the matrix
    # to calculate the L values
    for (j in 1:n) {
	for (i in 1:n) {
	    if (i != j){
		L[j] <- L[j] * (p - x[i])
	    }
	}
	for(i in 1:n){
	    if(i != j){
		L[j] <- L[j]/(x[j] - x[i])
	    }
	}
    }

    A <- 0
    for(j in 1:n){
	A <- A + y[j]*L[j]
    }

}
