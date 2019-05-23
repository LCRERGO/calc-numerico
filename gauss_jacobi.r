##Gauss-Jacobi

gauss_jacobi <- function(A,B,xk,erro){
    new_xk <- matrix(0,3,1)

    erro <- 0.00005
    n <- dim(B)[1]
    k <- 0
    d <- erro + 1
    iter <- 0

    flag <- 1  

    for (i in 1:n){
	sum <- 0
	for (j in 1:n){
	    if (i != j)
		sum <- sum + abs(A[i,j])
	}
	if (sum/abs(A[i,i]) >= 1)
	    flag <- 0
    }


    if (flag){
	start_time <- Sys.time()
	while(d > erro){
	    k <- k + 1
	    for (i in 1:n){
		sum <- 0
		for (j in 1:n){
		    if (i != j)
			sum <- sum + A[i,j] * xk[j]
		    iter <- iter + 1
		}
		new_xk[i] <- (1/A[i,i]) * (B[i]-sum)
	    }

	    max1 <- max(abs(new_xk-xk))
	    max2 <- max(abs(new_xk))
	    d <- max1/max2

	    xk <- new_xk
	}

	end_time <- Sys.time()
	timer <- end_time - start_time

	result <- c(new_xk, timer, iter)
	return (result)
    }else 
	return (NA)
}
