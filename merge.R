library(readr)
setwd("C:/Users/Carlos M Garcia/Dropbox/wences/CODIGO/matlab_crosscov/gretton")

finais_resposta <- 
  read_csv("C:/Users/Carlos M Garcia/Dropbox/wences/CODIGO/matlab_crosscov/gretton/finais_resposta.csv")

nomesmateriais=
  read.table("C:/Users/Carlos M Garcia/Dropbox/wences/metricas/nomesmateriais.txt", quote="\"", comment.char="")


nomes=sub("\\_.*", "",nomesmateriais$V1 )
nomesdata=data.frame(Formula=nomes)
resultado=merge(x = finais_resposta, y = nomesdata,by='Formula',all=FALSE) # Equivalente

n=length(nomes)
resposta=numeric()

resposta=0

for (i in 1:n)
  resposta[i]=sample(finais_resposta$'Formation Energy (eV)'[finais_resposta$Formula==nomes[i]],1)


setdiff(finais_resposta$Formula,nomes)


write.table(resposta, file="enerxia.txt", row.names=FALSE, col.names=FALSE)


