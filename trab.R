leitura_basedados <- function(arq,tag){
	dados_db <- read.csv(file=arq, header=TRUE, sep=tag)
	cat("\nTudo certo com a leitura do arquivo", arq, "\n")
	return(dados_db)
}

dados <- leitura_basedados("IBGE.csv", "")

print (dados)
print (dados[,1])
