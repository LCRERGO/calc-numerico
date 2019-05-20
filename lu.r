##Fatoracao LU

iteracoes = 0;
dim = 3;
U = matrix(c(1,2,3,2,1,-1,1,-1,-1),dim,dim);
b = matrix(c(3,0,2),dim,1);
L = matrix(0,dim,dim);
y = matrix(0,dim,1);
x = matrix(0,dim,1);

tempoI = Sys.time();
for(i in 1:dim){ ##Coluna
	for(j in 1:dim){ ##Linha

		iteracoes = iteracoes + 1;
		if(j>i ){
			m = U[j,i]/U[i,i];
			L[j,i] = round(m,7);
			U[j,] = round(U[j,] - m* U[i,],7); 
		}
		if(j==i)
			L[i,i]=1;
	}
}

U[j,] = round(U[j,] - (U[j,j-1]/U[j-1,j-1]) * U[j-1],7);
iteracoes = iteracoes + 1;


y[1] = b[1,1]/L[1,1];
iteracoes = iteracoes + 1;
for (i in 2:dim){
	soma = 0;
	for (j in 1:(i-1)){
		iteracoes = iteracoes + 1;
		soma = soma + (L[i,j] * y[j]);
	}
	y[i] = (b[i,1] - soma)/L[i,i];
}

x[dim] = y[dim]/U[dim,dim];
iteracoes = iteracoes + 1;
for (k in (dim-1):1){
	soma = 0;
	for (j in (k+1):dim){
		soma = soma + U[k,j]*x[j];
		iteracoes = iteracoes + 1;
	}
	x[k] = (y[k] - soma) / U[k,k];
}
tempoF = Sys.time();

tempo = tempoF - tempoI;

print(x);
cat(iteracoes, " iteracoes\n");
print(tempo);

