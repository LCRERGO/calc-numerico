    
leitura_basedados <- function(arq,tag){
	dados_db <- read.csv(file=arq, header=TRUE, sep=";")
	cat("\nTudo certo com a leitura do arquivo", arq, "\n")
	return(dados_db)
}

dados <- leitura_basedados("IBGE.csv", "")

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
    cat("Gauss Elimination Method: \n")
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

    print(timer)
    cat(iter," steps to get result\n")
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
    cat("LU Factoration Method: \n")
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

    print(timer)
    cat(iter," steps to get result\n")
    return (x)
}

# Name: gauss_jacobi.r
# Desc.: Solves a linear system with the
# Gauss-Jacobi elmination method.
# Group: 2
# 
# I: A matrix of the coeficients(A),
# a column matrix(b), a column matrix(xk) and a error value.
# O: A column matrix with the solution
# value(X)
# 
# Ex.: n <- 3
# A <- matrix(c(10,1,2,2,5,3,1,1,10), n, n)
# b <- matrix(c(7,-8,6), n, 1)
# xk <- matrix(c(0.7,1.4,-1.1),n,1)
# error <- 0.0005
#
# X <- matrix(c(0.9999815,-2.000022,0.99973), n, 1)

gauss_jacobi <- function(A,b,xk,error){
	cat("Gauss-Jacobi Method: \n")
	n <- dim(b)[1]
	xk_new <- matrix(0,n,1)
	iter <- 0
	d <- error + 1
	flag <- 1

	# Test if it is possible to solve
	for (i in 1:n){
		sum <- 0
		for (j in 1:n){
			if (i!=j)
				sum <- sum + abs(A[i,j])
		}
		if (sum/abs(A[i,i])>=1)
			flag <- 0
	}


	# The main solution
	if (flag){
		start_time <- Sys.time()
		while(d>error){
			for (i in 1:n){
				sum <- 0
				for (j in 1:n){
					if (i!=j){
					sum <- sum + A[i,j] * xk[j]
					iter <- iter + 1
					}
				}
			xk_new[i] <- (1/A[i,i]) * (b[i]-sum)
			}
	
		max1 <- max(abs(xk_new - xk))
		max2 <- max(abs(xk_new))
		d <- max1/max2

		xk <- xk_new
		}
	
		end_time <- Sys.time()
		timer <- end_time - start_time
		
		print(timer)
     		cat(iter," steps to get result")
    		return (xk_new)
	}else {
	cat("This cant be solved by the method\n")
	return(NaN)
	}
}



# Name: gauss-seidel.r
# Desc.: Solves a linear system with the
# Gauss-Seidel elmination method.
# Group: 2
# 
# Imput: A matrix of the coeficients(A),
# a column matrix(b), a column matrix(xk) and a error value.
# Output: A column matrix with the solution value(X)
# 

# Ex: n <- 3
# A = matrix(c(10,1,2,2,5,3,1,1,10), n, n)
# b = matrix(c(7,-8,6), n, 1)
# xk = matrix(c(0.7,1.4,-1.1), n, 1)
# error = 0.00005
#
# X <- matrix(c(0.9999815,-2.000022,0.99973), n, 1)

gauss_seidel <- function(A,B,xk,error){
    cat("Gauss-Seidel Method: \n")
    n = dim(B)[1]    
    new_xk <- matrix(0,n,1)
    iter = 0
    k = 0
    d = error + 1
    flag = 1  

    # Test if it is possible to solve
    for (i in 1:n){
	sum = 0
	for (j in 1:n){
	    if (i != j)
		sum = sum + abs(A[i,j])
	}
	if (sum/abs(A[i,i]) >= 1)
	    flag = 0
    }

    # The main solution
    if (flag){
		start_time = Sys.time()
		while(d > error){
		    for (i in 1:n){
				sum = 0;
				for (j in 1:n){
				    if (i != j){
				    	if(j<i)
				    		sum = sum + A[i,j] * new_xk[j]
				    	else
				    		sum = sum + A[i,j] * xk[j]
				    	}
				    iter = iter + 1
					}
				new_xk[i] = (1/A[i,i]) * (B[i]-sum)
			    }

	    max1 = max(abs(new_xk-xk))
	    max2 = max(abs(new_xk))
	    d = max1/max2

		xk = new_xk
		}

		end_time = Sys.time()
		timer = end_time - start_time

		print(timer)
    		cat(iter," steps to get result\n")
    		return (new_xk)
    }else {
	cat("This cant be solved by the method\n")
	return (NaN);
	}
}

# Name: interpolation.r
# Desc.: A set of functions to
# get the coeficients of a polynomial
# Group: 2
#
# I: A vector of the x elements(x)
# A vector of the y elements(y)
# O: A vector with the coefitients
# of the polynomial(A)
#
# Ex.:x <- c(0,1)
# y <- c(1.35,2.94)
#
# A <- c(1.35,1.59)

interpolation_vandermonde <- function(x,y,option){
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
	
	if(option == 1){
      	A <- gauss_elimination(aux, y)
		return (A)
	}
	if(option == 2){
		B <- lu_factoration(aux,y)
		return (B)
	}
	if (option == 3){
		C <- gauss_jacobi(aux,y,matrix(0,n,1),0.00005)
		return (C)
	}
	if (option == 4){
		D <- gauss_seidel(aux,y,matrix(0,n,1),0.00005)
    		return (D)
	}
}

calculate <- function(x,y){
	sum <- 0
	n <- dim(y)
	for (i in 1:n)
		sum <- sum + (x^(i-1))*y[i]
	return (sum)
}