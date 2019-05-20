##Gauss-Jacobi

A = matrix(c(10,1,2,2,5,3,1,1,10),3,3);
B = matrix(c(7,-8,6),3,1);
xk = matrix(c(0.7,1.4,-1.1),3, 1);
xk_novo = matrix(0,3,1);
	
erro = 0.00005;
n = 3;
k = 0;
d = erro + 1;

flag = 1;  

for (i in 1:n){
	soma = 0;
	for (j in 1:n){
		if (i!=j)
			soma = soma + abs(A[i,j]);
	}
	if (soma/abs(A[i,i])>=1)
		flag = 0;
}


if (flag){
	tempoI = Sys.time();
	while(d>erro){
	k = k+1;
	for (i in 1:n){
		soma = 0;
		for (j in 1:n){
			if (i!=j)
				soma = soma + A[i,j] * xk[j];
		}
		xk_novo[i] = (1/A[i,i]) * (B[i]-soma);
	}
	
	max1 = max(abs(xk_novo-xk));
	max2 = max(abs(xk_novo));
	d = max1/max2;

	xk = xk_novo;
	}
	
	tempoF= Sys.time()
	tempo = tempoF - tempoI;
	print(xk_novo);
	cat(k," iterações para a convergência\n");
	print(tempo);
}else 
	cat("Não existe solução por Gauss-Jacobi");
