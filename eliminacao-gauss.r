m=0
iteracoes = 0;
dim = 3;
a = matrix(c(1,2,3,2,1,-1,1,-1,-1),dim,dim);
b = matrix(c(3,0,-2),dim,1);
x = matrix(0,dim,1);

if(det(a)!=0){
tempoI = Sys.time();
for (k in 1:(dim-1)){
	for (i in (k+1):dim){
		m = a[i,k]/a[k,k];
		a[i,k] = 0;
		for ( j in (k+1):dim){
			iteracoes = iteracoes + 1;
			a[i,j]=(a[i,j]- (m*a[k,j]));
		}
		b[i] = b[i] - m*b[k];
	}
}

x[dim] = b[dim]/a[dim,dim];
for (k in (dim-1):1){
	soma = 0;
	for (j in (k+1):dim){
		soma = soma + a[k,j]*x[j];
		iteracoes = iteracoes + 1;
	}
	x[k] = (b[k] - soma) / a[k,k];
}
tempoF = Sys.time();
print(x);
cat(iteracoes," iteracoes\n");
tempo= tempoF - tempoI;
print(tempo);
} else{
	print("Não é possível realizar Eliminação de Gauss");
}