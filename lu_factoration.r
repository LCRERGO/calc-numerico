##Fatoracao LU

lu_factoration <- function(U,b){
    L <- matrix(0,dim,dim);
    y <- matrix(0,dim,1);
    x <- matrix(0,dim,1);

    start_time <- Sys.time()
    iter <- 0
    n <- dim(b)[1]

    for(i in 1:dim){ ##Coluna
	for(j in 1:dim){ ##Linha

	    iter <- iter + 1
	    if(j > i ){
		m <- U[j,i]/U[i,i];
		L[j,i] <- round(m,7);
		U[j,] <- round(U[j,] - m * U[i,],7); 
	    }
	    if(j == i)
		L[i,i]=1;
	}
    }

    U[j,] <- round(U[j,] - (U[j,j-1]/U[j-1,j-1]) * U[j-1],7);
    iter <- iter + 1;


    y[1] <- b[1,1]/L[1,1];
    iter <- iter + 1;
    for (i in 2:dim){
	sum <- 0;
	for (j in 1:(i-1)){
	    iter <- iter + 1;
	    sum <- sum + (L[i,j] * y[j]);
	}
	y[i] <- (b[i,1] - sum)/L[i,i];
    }

    x[dim] <- y[dim]/U[dim,dim];
    iter <- iter + 1;
    for (k in (dim-1):1){
	sum <- 0;
	for (j in (k+1):dim){
	    sum <- sum + U[k,j]*x[j];
	    iter <- iter + 1;
	}
	x[k] <- (y[k] - sum) / U[k,k];
    }
    end_time <- Sys.time();

    timer <- end_time - start_time;

    result <- c(x,timer,iter)
    return(result)
}
