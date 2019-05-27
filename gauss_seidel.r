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
# A = matrix(c(10,1,2,2,5,3,1,1,10), n, n);
# b = matrix(c(7,-8,6), n, 1);
# xk = matrix(c(0.7,1.4,-1.1), n, 1);
# error = 0.00005;
#
# X <- matrix(c(0.9999815,-2.000022,0.99973), n, 1)

gauss_seidel <- function(A,B,xk,error){
    xk_new <- matrix(0,3,1)
    iter = 0;
    n = dim(B)[1]
    k = 0
    d = error + 1;
    flag = 1;  

    # Test if it is possible to solve
    beta <- matrix(0,3,1)

   	for (i in 2:n) {
   		beta[1] = beta[1] + ( abs(A[1,i]) / abs(A[1,1]) ); 
   	}

   	for (i in 2:n) {
   		for (j in 1:i-1) {
   			beta[i] = beta[i] + (abs(A[i,j])*beta[j]);
   		}

   		for (j in i+1:n) {
   			beta[i] = beta[i] + abs(A[i,j]);
   		}

   		beta[i] = beta[i] / abs(A[i,i]);
   	}

   	b = max(beta);

   	if(b >= 1)
   		flag = 0;

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
				    iter = iter + 1;
					}
				new_xk[i] = (1/A[i,i]) * (B[i]-sum);
			    }

	    max1 = max(abs(new_xk-xk));
	    max2 = max(abs(new_xk));
	    d = max1/max2;

		xk = new_xk
		}

		end_time = Sys.time()
		timer = end_time - start_time

		retorno = c(new_xk, timer, iter)
		return (retorno)
    }else 
	cat("This cant be solved by the method");
}

A = matrix(c(10,1,2,2,5,3,1,1,10),3,3);
B = matrix(c(7,-8,6),3,1);
xk = matrix(c(0.7,1.4,-1.1),3, 1);
error = 0.00005;

cat(gauss_seidel(A,B,xk,error));