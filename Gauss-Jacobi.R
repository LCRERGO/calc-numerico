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

jacobi <- function(A,b,xk,error){
	xk_new <- matrix(0,3,1);
	iter = 0;
	n = dim(b)[1];
	d = error + 1;
	flag = 1;  

	# Test if it is possible to solve
	for (i in 1:n){
		sum = 0;
		for (j in 1:n){
			if (i!=j)
				sum = sum + abs(A[i,j]);
		}
		if (sum/abs(A[i,i])>=1)
			flag = 0;
	}


	# The main solution
	if (flag){
		start_time = Sys.time();
		while(d>error){
			for (i in 1:n){
				sum = 0;
				for (j in 1:n){
					if (i!=j){
					sum = sum + A[i,j] * xk[j];
					iter = iter + 1;
					}
				}
			xk_new[i] = (1/A[i,i]) * (b[i]-sum);
			}
	
		max1 = max(abs(xk_new - xk));
		max2 = max(abs(xk_new));
		d = max1/max2;

		xk = xk_new;
		}
	
		end_time = Sys.time()
		timer = end_time - start_time;
		
		retorno = c(xk_new,timer,iter);
		return (retorno);
	}else 
	cat("This cant be solved by the method");
}

A = matrix(c(10,1,2,2,5,3,1,1,10),3,3);
B = matrix(c(7,-8,6),3,1);
xk = matrix(c(0.7,1.4,-1.1),3,1);
erro = 0.0005;

cat(jacobi(A,B,xk,erro));


