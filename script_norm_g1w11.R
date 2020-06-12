####
if(!require(svDialogs)){install.packages("svDialogs")}
library(svDialogs)
if(!require(dplyr)){install.packages("dplyr")}
library(dplyr)
if(!require(car)){install.packages("car")}
library(car)
if(!require(tibble)){install.packages("tibble")}
library(tibble)
if(!require(ARTool)){install.packages("ARTool")}
library(ARTool)
if(!require(emmeans)){install.packages("emmeans")}
library(emmeans)





#################### 
# Funções das operações 

################### 

#funcao para operacao de divisao
div <- function(x1,x2){
  return(x1/x2)
}

#funcao para operacao de subtracao
sub <- function(x1,x2){
  return(x1-x2)
}

#funcao para operacao de radiciacao
root <- function(x1,x2){
  return(x1^(1/x2))
}

#funcao para operacao de logaritmo
logarit <- function(x1,x2){
  return(log(x1, base=x2))
}

#funcao z-score
zs <- function(x1,x2,x3){
  return((x1-x2)/x3)
}

#funcao pareto
pp <- function(x1,x2,x3){
  return((x1-x2)/sqrt(x3))
}

#funcao min-max
mm <- function(x1,x2,x3){
  return((x1-x2)/(x3-x2))
}

#funcao media-min-max
me <- function(x1,x2,x3,x4){
  return((x1-x2)/(x3-x4))
}


#função imputation
imputed.tab <- function(tabela, f.numfact) {
  tabela[tabela==0] <- 999999999999 #substitui todos os zeros por 999999999999 (um valor seguramente acima do máximo exportado pelo Progenesis )
  minimum <- min(tabela[(f.numfact+2):ncol(tabela)]) # agora que a tabela não tem mais zeros, descobre qual o mínimo
  tabela[tabela==999999999999] <- runif(sum(tabela==999999999999),min=minimum/1000, max=minimum) # substitui todos os 999999999999 por um aleatório entre o mínimo da tabela/1000 e  o mínimo da tabela.
  return(tabela)
}
  



############################
#  1A1

funcao1A1 <- function(tabela, f.numfact, f.numlevels, f.numreplic, FUN1, FUN2, FUN3, FUN4){
  primeira <- c()
  segunda <- c()
  terceira <- c()
  dados <- c()
  celulasTXT <- NULL
  for (i in 0:(f.numreplic-1)){ # primeiro loop para executar uma vez para cada replicata (cada replicata vai corresponder a um cálculo de mediana - ou do cálculo que for a vez...)
    for(k in 0:(prod(f.numlevels)-1)){ #dentro do primeiro loop faz um outro para fazer o cálculo (ex: mediana) com a "i" replicata de cada condição
      celulasTXT <- paste(celulasTXT,"tabela[", (i+(k*f.numreplic)+1), ",", (f.numfact+2),"]", sep="") # monta um texto com os endereços das células a serem usadas no cálculo
      if(k<(prod(f.numlevels)-1)){
        celulasTXT <- paste(celulasTXT,", ",sep="") # se não for a última adição ao texto, coloca vírgula no final
      }
    }
    calcobjTXT <- paste("c(",celulasTXT,")", sep="") # adiciona o "c(" para montar o comando de fazer o vetor usando as células selecionadas acima
    calcOBJ <- eval(parse(text=calcobjTXT)) # executa como comando o texto que foi gerado no loop para colocar em um vetor os valores da tabela que serão usados no cálculo
    primeira[i+1] <- FUN2(calcOBJ) # faz o cálculo e coloca o resultado num vetor para usar com o segundo cálculo depois.
    celulasTXT <- NULL # limpa a variável de montar o conjunto de células, preparando para o próximo conjunto (próxima volta no loop)
  }
  if(!is.null(FUN3)){
    for (i in 0:(f.numreplic-1)){ 
      for(k in 0:(prod(f.numlevels)-1)){ 
        celulasTXT <- paste(celulasTXT,"tabela[", (i+(k*f.numreplic)+1), ",", (f.numfact+2),"]", sep="") 
        if(k<(prod(f.numlevels)-1)){
          celulasTXT <- paste(celulasTXT,", ",sep="") 
        }
      }
      calcobjTXT <- paste("c(",celulasTXT,")", sep="") 
      calcOBJ <- eval(parse(text=calcobjTXT)) 
      segunda[i+1] <- FUN3(calcOBJ) 
      celulasTXT <- NULL 
    }
    
    if(!is.null(FUN4)){
      for (i in 0:(f.numreplic-1)){ 
        for(k in 0:(prod(f.numlevels)-1)){ 
          celulasTXT <- paste(celulasTXT,"tabela[", (i+(k*f.numreplic)+1), ",", (f.numfact+2),"]", sep="") 
          if(k<(prod(f.numlevels)-1)){
            celulasTXT <- paste(celulasTXT,", ",sep="") 
          }
        }
        calcobjTXT <- paste("c(",celulasTXT,")", sep="") 
        calcOBJ <- eval(parse(text=calcobjTXT)) 
        terceira[i+1] <- FUN4(calcOBJ) 
        celulasTXT <- NULL 
      }
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda, terceira)
    }else{
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda)
    }
  }else{
    dados <- FUN1(tabela[,(f.numfact+2)], primeira) # Aplica FUN1 utilizando a primeira coluna de proteína da tabela e o resultado da FUN2
  }
  temp <- data.frame(dados)
  nome <- names(tabela)[(f.numfact+2)]
  names(temp)[1] <- nome
  
  
  for(j in (f.numfact+3):ncol(tabela)){ #W substituído 5 por f.numfact+3 (começa a partir da segunda proteína, pois a primeira já foi feita) - loop repete para todas as outras proteínas
    for (i in 0:(f.numreplic-1)){ 
      for(k in 0:(prod(f.numlevels)-1)){ 
        celulasTXT <- paste(celulasTXT,"tabela[", (i+(k*f.numreplic)+1), ",", j,"]", sep="") 
        if(k<(prod(f.numlevels)-1)){
          celulasTXT <- paste(celulasTXT,", ",sep="") 
        }
      }
      calcobjTXT <- paste("c(",celulasTXT,")", sep="") 
      calcOBJ <- eval(parse(text=calcobjTXT)) 
      primeira[i+1] <- FUN2(calcOBJ) 
      celulasTXT <- NULL 
    }
    dados <- c()
    # caso seja necessario uma segunda estatistica
    if(!is.null(FUN3)){
      for (i in 0:(f.numreplic-1)){ 
        for(k in 0:(prod(f.numlevels)-1)){ 
          celulasTXT <- paste(celulasTXT,"tabela[", (i+(k*f.numreplic)+1), ",", j,"]", sep="") 
          if(k<(prod(f.numlevels)-1)){
            celulasTXT <- paste(celulasTXT,", ",sep="") 
          }
        }
        calcobjTXT <- paste("c(",celulasTXT,")", sep="") 
        calcOBJ <- eval(parse(text=calcobjTXT)) 
        segunda[i+1] <- FUN3(calcOBJ) 
        celulasTXT <- NULL 
      }
      
      # caso seja necessario uma terceira estatistica
      if(!is.null(FUN4)){
        for (i in 0:(f.numreplic-1)){ 
          for(k in 0:(prod(f.numlevels)-1)){ 
            celulasTXT <- paste(celulasTXT,"tabela[", (i+(k*f.numreplic)+1), ",", (f.numfact+2),"]", sep="") 
            if(k<(prod(f.numlevels)-1)){
              celulasTXT <- paste(celulasTXT,", ",sep="") 
            }
          }
          calcobjTXT <- paste("c(",celulasTXT,")", sep="") 
          calcOBJ <- eval(parse(text=calcobjTXT)) 
          terceira[i+1] <- FUN4(calcOBJ) 
          celulasTXT <- NULL 
        }
        dados <- FUN1(tabela[,j], primeira, segunda, terceira)
      }else{
        dados <- FUN1(tabela[,j], primeira, segunda)
      }
    }else{
      dados <- FUN1(tabela[,j], primeira)
    }
    temp <- data.frame(temp,dados)
    nome <- names(tabela)[j]
    names(temp)[(j-(f.numfact+1))] <- nome 
  }
  
  #W finaliza o data frame com as colunas de id e do nome de cada fator (ex: sex, activator)
  temp <- data.frame(tabela[,1:(f.numfact+1)], temp) #W substituído 3 por f.numfact+1
  return(temp)
}


############################
#   A11

funcaoA11 <- function(tabela, f.numfact, f.numlevels, f.numreplic, FUN1, FUN2, FUN3, FUN4){ #W incluídos num. fatores, níveis e replicatas. o NULL de FUN3 e FUN4 foi para a chamada da função
  # monta parte do texto para colocar no "list" da "primeira", "segunda" ou "terceira" dependendo do número de fatores
  txEvPaLIST <- character() 
  for (i in 1:f.numfact) {
    txEvPaLIST <- paste(txEvPaLIST,"tabela$",colnames(tabela)[i+1], sep="") # monta o "texto" da fórmula para aplicar a normalização, mas sem a vírgula entre os elementos
    if (i<f.numfact){
      txEvPaLIST <- paste(txEvPaLIST,", ", sep="") # acrescenta a vírgula e o espaço se não for o último elemento
    }
  }
  # primeira estatistica calculada (ex: media, mediana, soma)
  txEvPaPRI <- paste("as.vector(tapply(tabela[,", (f.numfact+2), "], list(", txEvPaLIST, "), FUN2))", sep="") #W monta o texto para "primeira" [as.vector, tapply, tabela, ver o tal do "4" (vai depender do número de fatores), txEvPaLIST montado acima, fun], coloca isso em txEvPaPRI
  primeira <- eval(parse(text = txEvPaPRI)) #W usa o texto de txEvPaPRI como um comando (como um argumento) e coloca o resultado na variável "primeira"  #ex: calcula a mediana para cada condição e coloca essas medianas na "primeira"
  #ex:                as.vector(tapply(tabela[,4], list(tabela$sex, tabela$activator), FUN2)) 
  dados <- c()
  #W constrói um vetor "conds" que contém o total de condições (produto de todos os níveis de cada fator) = no. de elementos do vetor, e número de replicatas = diferença entre os elementos do vetor 
  conds <- c()
  conds[1] <- 1
  for(i in 2:prod(f.numlevels)){
    conds[i] <- conds[i-1]+f.numreplic
  }
  # caso seja necessario uma segunda estatistica
  if(!is.null(FUN3)){
    txEvPaSEG <- paste("as.vector(tapply(tabela[,", (f.numfact+2), "], list(", txEvPaLIST, "), FUN3))", sep="") #W monta o texto para "segunda" [as.vector, tapply, tabela, ver o tal do "4" (vai depender do número de fatores), txEvPaLIST montado acima, fun], coloca isso em txEvPaSEG
    segunda <- eval(parse(text = txEvPaSEG))
    # caso seja necessario uma terceira estatistica
    if(!is.null(FUN4)){
      txEvPaTER <- paste("as.vector(tapply(tabela[,", (f.numfact+2), "], list(", txEvPaLIST, "), FUN4))", sep="") #W monta o texto para "terceira" [as.vector, tapply, tabela, ver o tal do "4" (vai depender do número de fatores), txEvPaLIST montado acima, fun], coloca isso em txEvPaTER
      terceira <- eval(parse(text = txEvPaTER))
      for(i in conds){
        dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), f.numfact+2], primeira[ceiling(i/f.numreplic)], segunda[ceiling(i/f.numreplic)], terceira[ceiling(i/f.numreplic)]) #W substituídos: 7 por f.numreplic-1; 8 por f.numreplic; 4 por f.numfact+2
      }
    }else{
      for(i in conds){
        dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), f.numfact+2], primeira[ceiling(i/f.numreplic)], segunda[ceiling(i/f.numreplic)])
      }
    }
  }else{
    for(i in conds){
      dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), f.numfact+2], primeira[ceiling(i/f.numreplic)]) #ex: divide o valor original (tabela) pela mediana(primeira), repete isso para as 8 replicatas ([i:(i+7), 4]) usando a mesma mediana ([ceiling(i/8)]), depois repete esse procedimento para as 4 condições (for(i in c(1,9,17,25)))
    }
  }
  # iniciando um data frame
  temp <- data.frame(dados)
  # passando o nome da proteina para o novo data frame
  nome <- names(tabela)[f.numfact+2]
  names(temp)[1] <- nome
  
  # agora que ja tem o data frame criado, loop para fazer para todas as proteinas
  for(j in (f.numfact+3):ncol(tabela)){ #W substituído 5 por f.numfact+3 (começa a partir da segunda proteína, pois a primeira já foi feita)
    #   if(j>ncol(tabela)){break} #W ---> VERIFICAR. O loop "for" está passando do limite estabelecido. Por exemplo, se o número de colunas é 7, ele executa mais uma vez com j=8. Com eta linha resolveu, mas é preciso ver por que. Novo teste, o problema não apareceu. Colocado parêntesis em volta de f.numfact+3.
    txEvPaPRI <- paste("as.vector(tapply(tabela[,", j, "], list(", txEvPaLIST, "), FUN2))", sep="") #W monta o texto para "primeira" [as.vector, tapply, tabela, ver o tal do "4" (vai depender do número de fatores), txEvPaLIST montado acima, fun], coloca isso em txEvPaPRI
    primeira <- eval(parse(text = txEvPaPRI)) #W usa o texto de txEvPaPRI como um comando (como um argumento) e coloca o resultado na variável "primeira"  #ex: calcula a mediana para cada condição e coloca essas medianas na "primeira"
    dados <- c()
    #W constrói um vetor "conds" que contém o total de condições (produto de todos os níveis de cada fator) = no. de elementos do vetor, e número de replicatas = diferença entre os elementos do vetor 
    conds <- c()
    conds[1] <- 1
    for(i in 2:prod(f.numlevels)){
      conds[i] <- conds[i-1]+f.numreplic
    }
    if(!is.null(FUN3)){
      txEvPaSEG <- paste("as.vector(tapply(tabela[,", j, "], list(", txEvPaLIST, "), FUN3))", sep="") #W monta o texto para "segunda" [as.vector, tapply, tabela, ver o tal do "4" (vai depender do número de fatores), txEvPaLIST montado acima, fun], coloca isso em txEvPaSEG
      segunda <- eval(parse(text = txEvPaSEG))
      if(!is.null(FUN4)){
        txEvPaTER <- paste("as.vector(tapply(tabela[,", j, "], list(", txEvPaLIST, "), FUN4))", sep="") #W monta o texto para "terceira" [as.vector, tapply, tabela, ver o tal do "4" (vai depender do número de fatores), txEvPaLIST montado acima, fun], coloca isso em txEvPaTER
        terceira <- eval(parse(text = txEvPaTER))
        for(i in conds){
          dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), j], primeira[ceiling(i/f.numreplic)], segunda[ceiling(i/f.numreplic)], terceira[ceiling(i/f.numreplic)]) #W substituídos: 7 por f.numreplic-1; 8 por f.numreplic
        }
      }else{
        for(i in conds){
          dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), j], primeira[ceiling(i/f.numreplic)], segunda[ceiling(i/f.numreplic)])
        }
      }
    }else{
      for(i in conds){
        dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), j], primeira[ceiling(i/f.numreplic)])
      }
    }
    temp <- data.frame(temp,dados)
    nome <- names(tabela)[j]
    names(temp)[(j-(f.numfact+1))] <- nome #W substituído 3 por f.numfact+1
  }
  
  # finaliza o data frame com as colunas de id, sex e activator
  temp <- data.frame(tabela[,1:(f.numfact+1)], temp)  #W substituído 3 por f.numfact+1
  return(temp)
}


############################
#  11A

funcao11A <- function(tabela, f.numfact, f.numlevels, f.numreplic, FUN1, FUN2, FUN3, FUN4){
  primeira <- c()
  segunda <- c()
  terceira <- c()
  dados <- c()
  for (i in 1:nrow(tabela)){
    primeira[i] <- FUN2(unlist(tabela[i,(f.numfact+2):ncol(tabela)])) # loop para fazer a mediana de cada linha e colocar num verot (primeira) # linhas contêm replicatas e condições; colunas contêm proteínas. # agrupamentos que terminam com A (11A, 1AA, A1A, AAA) devem usar (f.numfact+2):ncol(tabela)
  }
  if(!is.null(FUN3)){
    for (i in 1:nrow(tabela)){
      segunda[i] <- FUN3(unlist(tabela[i,(f.numfact+2):ncol(tabela)]))
    }
    if(!is.null(FUN4)){
      for (i in 1:nrow(tabela)){
        terceira[i] <- FUN4(unlist(tabela[i,(f.numfact+2):ncol(tabela)]))
      }
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda, terceira)
    }else{
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda)
    }
  }else{
    dados <- FUN1(tabela[,(f.numfact+2)], primeira) # Aplica FUN1 utilizando a primeira coluna de proteína da tabela e o resultado da FUN2
  }
  temp <- data.frame(dados)
  nome <- names(tabela)[(f.numfact+2)]
  names(temp)[1] <- nome
  
  
  for(j in (f.numfact+3):ncol(tabela)){ #W substituído 5 por f.numfact+3 (começa a partir da segunda proteína, pois a primeira já foi feita) - loop repete para todas as outras proteínas
    for (i in 1:nrow(tabela)){
      primeira[i] <- FUN2(unlist(tabela[i,(f.numfact+2):ncol(tabela)]))
    }
    dados <- c()
    # caso seja necessario uma segunda estatistica
    if(!is.null(FUN3)){
      for (i in 1:nrow(tabela)){
        segunda[i] <- FUN3(unlist(tabela[i,(f.numfact+2):ncol(tabela)]))
      }
      # caso seja necessario uma terceira estatistica
      if(!is.null(FUN4)){
        for (i in 1:nrow(tabela)){
          terceira[i] <- FUN4(unlist(tabela[i,(f.numfact+2):ncol(tabela)]))
        }
        dados <- FUN1(tabela[,j], primeira, segunda, terceira)
      }else{
        dados <- FUN1(tabela[,j], primeira, segunda)
      }
    }else{
      dados <- FUN1(tabela[,j], primeira)
    }
    temp <- data.frame(temp,dados)
    nome <- names(tabela)[j]
    names(temp)[(j-(f.numfact+1))] <- nome #W substituído 3 por f.numfact+1
  }
  
  #W finaliza o data frame com as colunas de id e do nome de cada fator (ex: sex, activator)
  temp <- data.frame(tabela[,1:(f.numfact+1)], temp) #W substituído 3 por f.numfact+1
  return(temp)
}


############################
#  1AA

funcao1AA <- function(tabela, f.numfact, f.numlevels, f.numreplic, FUN1, FUN2, FUN3, FUN4){
  conds <- c()
  primeira <- c()
  segunda <- c()
  terceira <- c()
  dados <- c()
  calcOBJ <- c()
  conds[1] <- 1
  for(i in 2:prod(f.numlevels)){
    conds[i] <- conds[i-1]+f.numreplic # define a linha de início de cada condição (produto dos números de níveis de todos os fatores)
  }
  
  for (k in 0:(f.numreplic-1)) {  # loop para repetir p/ cada replicata
    for (i in conds) { # loop para usar todas as condições de cada replicata
      calcOBJ <- c(calcOBJ, unlist(tabela[(i+k),(f.numfact+2):ncol(tabela)])) # coloca todos os valores que serão usdos para calcular a função (ex:mediana) no vetor calcOBJ
    }
    primeira[k+1] <- FUN2(calcOBJ) # calcula a FUN2 (ex:mediana) e coloca no vetor primeira. Cada mediana colocada aqui representa os valores de todas as proteínas, todas as condiçòes e 1 replicata. 
    calcOBJ <- c()
  }
  
  if(!is.null(FUN3)){
    for (k in 0:(f.numreplic-1)) {  
      for (i in conds) {
        calcOBJ <- c(calcOBJ, unlist(tabela[(i+k),(f.numfact+2):ncol(tabela)]))
      }
      segunda[k+1] <- FUN3(calcOBJ) 
      calcOBJ <- c()
    }
    
    if(!is.null(FUN4)){
      for (k in 0:(f.numreplic-1)) {  
        for (i in conds) {
          calcOBJ <- c(calcOBJ, unlist(tabela[(i+k),(f.numfact+2):ncol(tabela)]))
        }
        terceira[k+1] <- FUN4(calcOBJ) 
        calcOBJ <- c()
      }
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda, terceira)
    }else{
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda)
    }
  }else{
    dados <- FUN1(tabela[,(f.numfact+2)], primeira) # Aplica FUN1 utilizando a primeira coluna de proteína da tabela e o resultado da FUN2
  }
  temp <- data.frame(dados)
  nome <- names(tabela)[(f.numfact+2)]
  names(temp)[1] <- nome
  
  
  for(j in (f.numfact+3):ncol(tabela)){ #W substituído 5 por f.numfact+3 (começa a partir da segunda proteína, pois a primeira já foi feita) - loop repete para todas as outras proteínas
    dados <- c()
    # caso seja necessario uma segunda estatistica
    if(!is.null(FUN3)){ # não precisa calcular novamente, pois primeira, segunda e terceira são as mesmas para qualquer coluna, portanto não dependem de j
      # caso seja necessario uma terceira estatistica
      if(!is.null(FUN4)){
        dados <- FUN1(tabela[,j], primeira, segunda, terceira)
      }else{
        dados <- FUN1(tabela[,j], primeira, segunda)
      }
    }else{
      dados <- FUN1(tabela[,j], primeira)
    }
    temp <- data.frame(temp,dados)
    nome <- names(tabela)[j]
    names(temp)[(j-(f.numfact+1))] <- nome 
  }
  
  #W finaliza o data frame com as colunas de id e do nome de cada fator (ex: sex, activator)
  temp <- data.frame(tabela[,1:(f.numfact+1)], temp) #W substituído 3 por f.numfact+1
  return(temp)
}


############################
#   A1A

funcaoA1A <- function(tabela, f.numfact, f.numlevels, f.numreplic, FUN1, FUN2, FUN3, FUN4){
  primeira <- c()
  segunda <- c()
  terceira <- c()
  dados <- c()
  conds <- c()
  for (i in 1:prod(f.numlevels)){
    primeira[i] <- FUN2(unlist(tabela[(f.numreplic*(i-1)+1):(f.numreplic*i),(f.numfact+2):ncol(tabela)]))
  }
  conds[1] <- 1
  for(i in 2:prod(f.numlevels)){
    conds[i] <- conds[i-1]+f.numreplic
  }
  # caso seja necessario uma segunda estatistica
  if(!is.null(FUN3)){
    for (i in 1:prod(f.numlevels)){
      segunda[i] <- FUN3(unlist(tabela[(f.numreplic*(i-1)+1):(f.numreplic*i),(f.numfact+2):ncol(tabela)]))
    }
    # caso seja necessario uma terceira estatistica
    if(!is.null(FUN4)){
      for (i in 1:prod(f.numlevels)){
        terceira[i] <- FUN4(unlist(tabela[(f.numreplic*(i-1)+1):(f.numreplic*i),(f.numfact+2):ncol(tabela)]))
      }
      for(i in conds){
        dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), f.numfact+2], primeira[ceiling(i/f.numreplic)], segunda[ceiling(i/f.numreplic)], terceira[ceiling(i/f.numreplic)]) #W substituídos: 7 por f.numreplic-1; 8 por f.numreplic; 4 por f.numfact+2
      }
    }else{
      for(i in conds){
        dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), f.numfact+2], primeira[ceiling(i/f.numreplic)], segunda[ceiling(i/f.numreplic)])
      }
    }
  }else{
    for(i in conds){
      dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), f.numfact+2], primeira[ceiling(i/f.numreplic)]) #ex: divide o valor original (tabela) pela mediana(primeira), repete isso para as 8 replicatas ([i:(i+7), 4]) usando a mesma mediana ([ceiling(i/8)]), depois repete esse procedimento para as 4 condições (for(i in c(1,9,17,25)))
    }
  }
  # iniciando um data frame
  temp <- data.frame(dados)
  # passando o nome da proteina para o novo data frame
  nome <- names(tabela)[f.numfact+2]
  names(temp)[1] <- nome
  
  # agora que ja tem o data frame criado, loop para fazer para todas as proteinas
  for(j in (f.numfact+3):ncol(tabela)){ #W substituído 5 por f.numfact+3 (começa a partir da segunda proteína, pois a primeira já foi feita)
    #   if(j>ncol(tabela)){break} #W ---> VERIFICAR. O loop "for" está passando do limite estabelecido. Por exemplo, se o número de colunas é 7, ele executa mais uma vez com j=8. Com eta linha resolveu, mas é preciso ver por que. Novo teste, o problema não apareceu. Colocado parêntesis em volta de f.numfact+3.
    for (i in 1:prod(f.numlevels)){
      primeira[i] <- FUN2(unlist(tabela[(f.numreplic*(i-1)+1):(f.numreplic*i),(f.numfact+2):ncol(tabela)]))
    }
    dados <- c()
    #W constrói um vetor "conds" que contém o total de condições (produto de todos os níveis de cada fator) = no. de elementos do vetor, e número de replicatas = diferença entre os elementos do vetor 
    conds <- c()
    conds[1] <- 1
    for(i in 2:prod(f.numlevels)){
      conds[i] <- conds[i-1]+f.numreplic
    }
    if(!is.null(FUN3)){
      for (i in 1:prod(f.numlevels)){
        segunda[i] <- FUN3(unlist(tabela[(f.numreplic*(i-1)+1):(f.numreplic*i),(f.numfact+2):ncol(tabela)]))
      }
      if(!is.null(FUN4)){
        for (i in 1:prod(f.numlevels)){
          terceira[i] <- FUN4(unlist(tabela[(f.numreplic*(i-1)+1):(f.numreplic*i),(f.numfact+2):ncol(tabela)]))
        }
        for(i in conds){
          dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), j], primeira[ceiling(i/f.numreplic)], segunda[ceiling(i/f.numreplic)], terceira[ceiling(i/f.numreplic)]) #W substituídos: 7 por f.numreplic-1; 8 por f.numreplic
        }
      }else{
        for(i in conds){
          dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), j], primeira[ceiling(i/f.numreplic)], segunda[ceiling(i/f.numreplic)])
        }
      }
    }else{
      for(i in conds){
        dados[i:(i+f.numreplic-1)] <- FUN1(tabela[i:(i+f.numreplic-1), j], primeira[ceiling(i/f.numreplic)])
      }
    }
    temp <- data.frame(temp,dados)
    nome <- names(tabela)[j]
    names(temp)[(j-(f.numfact+1))] <- nome #W substituído 3 por f.numfact+1
  }
  
  # finaliza o data frame com as colunas de id, sex e activator
  temp <- data.frame(tabela[,1:(f.numfact+1)], temp)  #W substituído 3 por f.numfact+1
  return(temp)
}


#######################################
#   AA1

funcaoAA1 <- function(tabela, f.numfact, f.numlevels, f.numreplic, FUN1, FUN2, FUN3, FUN4){ #W incluídos num. fatores, níveis e replicatas. o NULL de FUN3 e FUN4 foi para a chamada da função
  # primeira estatistica calculada (ex: media, mediana, soma)
  primeira <- FUN2(tabela[,(f.numfact+2)]) #W o "4" (primeira coluna contendo dados de ptn) foi substituído pelo número de fatores+2 (ex: 1a col = rep id, 2a col = 1o fator, 3a col = 2o fator, 4a col = 1a ptn). Aplica FUN2 à primeira coluna que tenha dados de proteínas.
  dados <- c()
  # caso seja necessario uma segunda estatistica
  if(!is.null(FUN3)){
    segunda <- FUN3(tabela[,(f.numfact+2)])
    # caso seja necessario uma terceira estatistica
    if(!is.null(FUN4)){
      terceira <- FUN4(tabela[,(f.numfact+2)])
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda, terceira)
    }else{
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda)
    }
  }else{
    dados <- FUN1(tabela[,(f.numfact+2)], primeira) # Aplica FUN1 utilizando a primeira coluna de proteína da tabela e o resultado da FUN2
  }
  # iniciando um data frame
  temp <- data.frame(dados)
  # passando o nome da proteina para o novo data frame
  nome <- names(tabela)[(f.numfact+2)]
  names(temp)[1] <- nome
  
  for(j in (f.numfact+3):ncol(tabela)){ #W substituído 5 por f.numfact+3 (começa a partir da segunda proteína, pois a primeira já foi feita) - loop repete para todas as outras proteínas
    primeira <- FUN2(tabela[,j])
    dados <- c()
    # caso seja necessario uma segunda estatistica
    if(!is.null(FUN3)){
      segunda <- FUN3(tabela[,j])
      # caso seja necessario uma terceira estatistica
      if(!is.null(FUN4)){
        terceira <- FUN4(tabela[,j])
        dados <- FUN1(tabela[,j], primeira, segunda, terceira)
      }else{
        dados <- FUN1(tabela[,j], primeira, segunda)
      }
    }else{
      dados <- FUN1(tabela[,j], primeira)
    }
    temp <- data.frame(temp,dados)
    nome <- names(tabela)[j]
    names(temp)[(j-(f.numfact+1))] <- nome #W substituído 3 por f.numfact+1
  }
  
  #W finaliza o data frame com as colunas de id e do nome de cada fator (ex: sex, activator)
  temp <- data.frame(tabela[,1:(f.numfact+1)], temp) #W substituído 3 por f.numfact+1
  return(temp)
}


#######################################
#   AAA

funcaoAAA <- function(tabela, f.numfact, FUN1, FUN2, FUN3, FUN4){ 
  
  # primeira estatistica calculada (ex: media, mediana, soma)
  primeira <- FUN2(unlist(tabela[(f.numfact+2):ncol(tabela)])) #W  Aplica FUN2 a todas as colunas que tenha dados de proteínas  (e todas as respectivas linhas). Ex: faz a mediana de todos os valores e coloca em "primeira"
  dados <- c()
  # caso seja necessario uma segunda estatistica
  if(!is.null(FUN3)){
    segunda <- FUN3(unlist(tabela[(f.numfact+2):ncol(tabela)]))
    # caso seja necessario uma terceira estatistica
    if(!is.null(FUN4)){
      terceira <- FUN4(unlist(tabela[(f.numfact+2):ncol(tabela)]))
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda, terceira)
    }else{
      dados <- FUN1(tabela[,(f.numfact+2)], primeira, segunda)
    }
  }else{
    dados <- FUN1(tabela[,(f.numfact+2)], primeira) # Aplica FUN1 utilizando a primeira coluna de proteína da tabela e o resultado da FUN2
  }
  # iniciando um data frame
  temp <- data.frame(dados)
  # passando o nome da proteina para o novo data frame
  nome <- names(tabela)[(f.numfact+2)]
  names(temp)[1] <- nome
  
  for(j in (f.numfact+3):ncol(tabela)){ #W substituído 5 por f.numfact+3 (começa a partir da segunda proteína, pois a primeira já foi feita) - loop repete para todas as outras proteínas
    primeira <- FUN2(unlist(tabela[(f.numfact+2):ncol(tabela)]))
    dados <- c()
    # caso seja necessario uma segunda estatistica
    if(!is.null(FUN3)){
      segunda <- FUN3(unlist(tabela[(f.numfact+2):ncol(tabela)]))
      # caso seja necessario uma terceira estatistica
      if(!is.null(FUN4)){
        terceira <- FUN4(unlist(tabela[(f.numfact+2):ncol(tabela)]))
        dados <- FUN1(tabela[,j], primeira, segunda, terceira)
      }else{
        dados <- FUN1(tabela[,j], primeira, segunda)
      }
    }else{
      dados <- FUN1(tabela[,j], primeira)
    }
    temp <- data.frame(temp,dados)
    nome <- names(tabela)[j]
    names(temp)[(j-(f.numfact+1))] <- nome #W substituído 3 por f.numfact+1
  }
  
  #W finaliza o data frame com as colunas de id e do nome de cada fator (ex: sex, activator)
  temp <- data.frame(tabela[,1:(f.numfact+1)], temp) #W substituído 3 por f.numfact+1
  return(temp)
}




############################
#W funcao para alimentar o dataframe de estatística comparativa entre as normalizações
#W  o dataframe inicial (só as colunas vazias) tem que ser criado no código principal, fora da função.
#W a cada método de normalização calculado, esta função é chamada novamente e acrescenta linhas ao dataframe
funcaoCompar <- function(f.data, f.numfactors, f.numlevels, f.numreplic, f.tst_thrsh, f.NormStat, f.empty.NormStat, f.bal.design, f.rep.meas, f.norm.meth) {
  
  #W A cada teste calculado, as respectivas linhas contendo os resultados são acrescentadas ao dataframe f.NormStat
  #W Calcula teste t, 1-way ou 2-way ANOVA p/ cada ptn (depende do desenho experimental, se tiver 2 fatores, 2-way, se tiver 1 fator c/ mais que 2 níveis, 1-way, se tiver 1 fator e 2 níveis, teste t)
  #W Calcula também testes de Whapiro-Wilk (normalidade) e de Levene (homogeneidade de variâncias) p/ cada ptn para se avaliar se ANOVA ou t são adequados ou se convém usar outros testes
  #W Calcula 2-way ANOVA não paramétrico (adequado para amostras sem distribuição normal), para desenhos com medidas isoladas ou com medidas repetidas
  
  
  #W dados para teste - DELETAR!!
  #W
  # f.data <- raw_data
  # f.numfactors <- 1
  # f.numlevels <- 4
  # f.numreplic <- 8
  # f.tst_thrsh <- 0.05
  # f.empty.NormStat <- empty.NormStat
  # f.bal.design <- bal.design
  # f.norm.meth <- "TESTE"
  #W
  #W
  
  
  
  if (f.numfactors == 1) {
    prg <- f.data #recebe os dados normalizados por alguma das funções de normalização (ou um dos dados direto do Progenesis)
    factor1 <- colnames(prg)[2] #recupera o nome do fator 1 
    nrep <- nrow(prg) # considera como número de repetições o número de linhas da análise atual
    nreptest <- (f.numlevels[[1]])*f.numreplic # testa se esse número de linhas é coerente com replictas * condições 
    if (nrep != f.numreplic*f.numlevels[[1]]){
      print("Check the experimental design. The number of rows in the file does not match the number of replicates * number of levels in each factor")
    }
    f.paired.t <- FALSE
    if (f.rep.meas == "Y" | f.rep.meas == "y") {
      f.paired.t <- TRUE
    }
    accs <- colnames(prg)[-(1:2)] # recupera o código de acesso das proteínas (nomes das colunas, exceto as 2 primeiras, que são ID da amostra e condição)
    if (f.numlevels[1] > 2) {
      #
      #
      ############# ONE WAY ANOVA   +Tukey HSD   +Shapiro-Wilk  +Levene       (http://www.sthda.com/english/wiki/wiki.php?title=one-way-anova-test-in-r)
      #
      #
      aovs <- lapply(accs, function(prot_acc) {  #cada termo de accs vai ser tratado como prot_acc neste grupo lapply
        formANOVA1 <- formula(paste(prot_acc, " ~ ", factor1)) # monta o "texto" da fórmula do 1-way ANOVA
        res.aov1 <- aov(formANOVA1, data = prg) # aplica o teste 1-way ANOVA ao dataframe prg
        pvalues.res.aov1 <- summary(res.aov1)[[1]][5] #extrai a parte do summary que tem o p-valor (a coluna 5 do primeiro elemento)
        pvalues.res.aov1 <- pvalues.res.aov1[-c(2), , drop = FALSE] #tira a linha excedente (segunda linha), para ficar só com o p-valor
        names(pvalues.res.aov1)[1] <- prot_acc # coloca o respectivo nome da proteína no nome da coluna referente ao p-valor
        pvalues.res.aov1 <- rownames_to_column(pvalues.res.aov1)
        res.Tukey <- TukeyHSD(res.aov1)
        pval.res.Tukey <- data.frame(res.Tukey [[1]][,4])
        names(pval.res.Tukey)[1] <- prot_acc # coloca o respectivo nome da proteína no nome da coluna referente ao p-valor
        pval.res.Tukey <- rownames_to_column(pval.res.Tukey)
        if(pvalues.res.aov1<f.tst_thrsh){  # Os p-valores do teste Tukey HSD só são incorporados à tabela de resultados se o ANOVA fro significativo, já que o Tukey é um pós-teste. Caso contrário, é colocado o valor 999
          #just continue
        }else{
          pval.res.Tukey [,2] <- 999
        }
        pvalues.res.aov1 <- bind_rows(pvalues.res.aov1, pval.res.Tukey)
        pvalues.res.aov1 <- column_to_rownames(pvalues.res.aov1)
        pvalues.res.aov1 # coloca isso no aovs
      })
      pvals <- do.call(cbind, aovs) %>% add_rownames("subdiv_test") #ajusta a estrutura do aovs para um dataframe
      part.NormStat <- bind_rows(f.empty.NormStat,pvals) #coloca  a linha dos p-valores na estrutura do f.NormStat - o dataframe geral das comparações deste script, mas ainda um arquivo parcial
      part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]<f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor menor que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
      part.NormStat$normaliz_meth <- f.norm.meth #W Achar um jeito de mandar o nome do método automaticamente para a função #coloca o nome da estratégia de normalização em todas as linhas desta parcial
      part.NormStat$stat_test <- as.character(part.NormStat$stat_test)
      part.NormStat$stat_test [1] <- "1w_AOV" # coloca o nome do teste na primeira linha desta parcial (que tem o p-valor do ANOVA)
      part.NormStat$stat_test [2:nrow(part.NormStat)] <- "TukeyHSD" # coloca o nome do teste nas outras linhas desta parcial (que tem os p-valores do Tukey HSD)
      part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
      score.sum.par <- sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes para o método de normalização que estiver sendo avaliado.
      f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
      part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
      #
      #W teste de Shapiro Wilk (normalidade)       p/ cada ptn --- If p-value > signif threshold, then normality is not likely to be violated.
      #
      Shapiro <- lapply(accs, function(prot_acc) {
        formANOVA1 <- formula(paste(prot_acc, " ~ ", factor1))
        res.aov1 <- aov(formANOVA1, data = prg)
        aov_residuals <- residuals(object = res.aov1)
        res.shap <- shapiro.test(x=aov_residuals)
        pvalues.res.shap <- data.frame(res.shap[[2]][1])
        names(pvalues.res.shap)[1] <- prot_acc
        pvalues.res.shap
      })
      Spvals <- do.call(cbind, Shapiro) %>% add_rownames("stat_test")
      Spvals[1,1] <-  "Shapiro-Wilk test" 
      part.NormStat <- bind_rows(f.empty.NormStat,Spvals) #coloca  as linhas dos p-valores na estrutura do f.NormStat - o dataframe geral das comparações deste script, mas ainda um arquivo parcial
      part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]>f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor maior que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
      part.NormStat$normaliz_meth <- f.norm.meth #W Achar um jeito de mandar o nome do método automaticamente para a função #coloca o nome da estratégia de normalização em todas as linhas desta parcial
      part.NormStat$stat_test  <- "Shap-Wilk" # coloca o nome do teste em todas as linhas desta parcial
      part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
      score.sum.par <- score.sum.par + sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes para o método de normalização que estiver sendo avaliado.
      f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
      part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
      #
      #W  teste de Levene (homogeneidade de variâncias)              p/ cada ptn  --- If p-value > signif threshold, then variances are likely to be the same across the conditions.
      # This test asumes normality, so the normality (Shapiro-Wilk) test should be performed before this test.
      #
      Levs <- lapply(accs, function(prot_acc) {
        formLevene <- formula(paste(prot_acc, " ~ ", factor1))
        res.Levene <- leveneTest(formLevene, data = prg)
        pvalues.res.Levene <- data.frame(res.Levene[[3]][1])
        names(pvalues.res.Levene)[1] <- prot_acc
        pvalues.res.Levene
      })
      Lpvals <- do.call(cbind, Levs) %>% add_rownames("stat_test")
      Lpvals[1,1] <-  "Levene test" 
      part.NormStat <- bind_rows(f.empty.NormStat,Lpvals) #coloca  as linhas dos p-valores na estrutura do f.NormStat - o dataframe geral das comparações deste script, mas ainda um arquivo parcial
      part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]>f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor maior que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
      part.NormStat$normaliz_meth <- f.norm.meth #W Achar um jeito de mandar o nome do método automaticamente para a função #coloca o nome da estratégia de normalização em todas as linhas desta parcial
      part.NormStat$stat_test  <- "Levene" # coloca o nome do teste em todas as linhas desta parcial
      part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
      score.sum.par <- score.sum.par + sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes para o método de normalização que estiver sendo avaliado.
      part.NormStat$score_sum_par[1] <- score.sum.par
      f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
      part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
      #W
      # Welch test (one way) não assume homogeneidade de variâncias  (para ser usado se não passou no teste de homogeneidade de variâncias, Levene)
      # Friedman test (one way repeated measures) se uma das condições do ANOVA não foi atendida (ausência de outliers significantes, normalidade, esfericidade (homogeneidade de variâncias). )
      # Unfortunately, there are no non-parametric alternatives to the two-way and the three-way repeated measures ANOVA. Thus, in the situation where the assumptions are not met, you could consider running the two-way/three-way repeated measures ANOVA on the transformed and non-transformed data to see if there are any meaningful differences. https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
      # Kruskal-Wallis rank sum test (não paramétrico) não assume normalidade (para ser usado se não passou no teste de normalidade, Shapiro-Wilk)
      # Desenho experimental desbalanceado (diferentes números de replicatas em cada condição) p/ 1-way ANOVA não é problemático, desde que não afete a homogeneidade de variâncias. Como isso está sendo testado, não é necessário inserir outro teste. (https://www.theanalysisfactor.com/when-unequal-sample-sizes-are-and-are-not-a-problem-in-anova/)
    } else {
      #
      #
      ############# T-TEST   
      #
      #
      #W se tem só 1 fator e 2 níveis, então usa o teste t
      Ttst <- lapply(accs, function(prot_acc) {
        formTtest <- formula(paste(prot_acc, " ~ ", factor1))
        res.Ttest <- t.test(formTtest, data = prg, var.equal = TRUE, paired = f.paired.t)  #W --- AJUSTAR --- mudar o "var.equal" em função do resultado do teste de variância. Se for FALSE, o teste feito será o de Welch. A resolver: Faz teste t para todas, depois Welch para todas ou deixa o script decidir e faz t para variância homogênea e Welch para variância não homogênea? 
        pvalues.res.Ttest <- data.frame(res.Ttest$p.value)
        names(pvalues.res.Ttest)[1] <- prot_acc
        pvalues.res.Ttest
      })
      Tpvals <- do.call(cbind, Ttst) %>% add_rownames("stat_test")
      Tpvals[1,1] <-  "t test" 
      part.NormStat <- bind_rows(f.empty.NormStat,Tpvals) #coloca  as linhas dos p-valores na estrutura do f.NormStat - o dataframe geral das comparações deste script, mas ainda um arquivo parcial
      part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]<f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor maior que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
      part.NormStat$normaliz_meth <- f.norm.meth #W Achar um jeito de mandar o nome do método automaticamente para a função #coloca o nome da estratégia de normalização em todas as linhas desta parcial
      part.NormStat$stat_test  <- "t-test" # coloca o nome do teste em todas as linhas desta parcial
      part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
      score.sum.par <- sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes para o método de normalização que estiver sendo avaliado.
      f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
      part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
      #
      #W  teste de Shapiro Wilk (normalidade)                  p/ cada ptn --- If p-value > signif threshold, then normality is not likely to be violated.
      #
      Shapiro <- lapply(accs, function(prot_acc) {
        tx.shap <- paste("shapiro.test(prg$",prot_acc,")")
        res.shap <- eval(parse(text = tx.shap))
        pvalues.res.shap <- data.frame(data.frame(res.shap$p.value))
        names(pvalues.res.shap)[1] <- prot_acc
        pvalues.res.shap
      })
      Spvals <- do.call(cbind, Shapiro) %>% add_rownames("stat_test")
      Spvals[1,1] <-  "Shapiro-Wilk test" 
      part.NormStat <- bind_rows(f.empty.NormStat,Spvals) #coloca  as linhas dos p-valores na estrutura do f.NormStat - o dataframe geral das comparações deste script, mas ainda um arquivo parcial
      part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]>f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor maior que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
      part.NormStat$normaliz_meth <- f.norm.meth #W Achar um jeito de mandar o nome do método automaticamente para a função #coloca o nome da estratégia de normalização em todas as linhas desta parcial
      part.NormStat$stat_test  <- "Shap-Wilk" # coloca o nome do teste em todas as linhas desta parcial
      part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
      score.sum.par <- score.sum.par + sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes para o método de normalização que estiver sendo avaliado.
      f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
      part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
      #
      #W Calcula teste de Levene (homogeneidade de variâncias)              p/ cada ptn  --- If p-value > signif threshold, then variances are likely to be the same across the conditions.
      # This test asumes normality, so the normality (Shapiro-Wilk) test should be performed before this test.
      #
      Levs <- lapply(accs, function(prot_acc) {
        formLevene <- formula(paste(prot_acc, " ~ ", factor1))
        res.Levene <- leveneTest(formLevene, data = prg)
        pvalues.res.Levene <- data.frame(res.Levene[[3]][1])
        names(pvalues.res.Levene)[1] <- prot_acc
        pvalues.res.Levene
      })
      Lpvals <- do.call(cbind, Levs) %>% add_rownames("stat_test")
      Lpvals[1,1] <-  "Levene test" 
      part.NormStat <- bind_rows(f.empty.NormStat,Lpvals) #coloca  as linhas dos p-valores na estrutura do f.NormStat - o dataframe geral das comparações deste script, mas ainda um arquivo parcial
      part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,6:ncol(part.NormStat)]>f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor maior que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
      part.NormStat$normaliz_meth <- f.norm.meth #W Achar um jeito de mandar o nome do método automaticamente para a função #coloca o nome da estratégia de normalização em todas as linhas desta parcial
      part.NormStat$stat_test  <- "Levene" # coloca o nome do teste em todas as linhas desta parcial
      part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
      score.sum.par <- score.sum.par + sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes para o método de normalização que estiver sendo avaliado.
      part.NormStat$score_sum_par[1] <- score.sum.par
      f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
      part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
    }
  } else {
    if(f.numfactors == 2) {
      #
      #
      ############# TWO-WAY ANOVA   +Tukey  +Shapiro-Wilk  +Levene       (http://www.sthda.com/english/wiki/two-way-anova-test-in-r#compute-two-way-anova-test)
      #
      #
      if(f.bal.design == "n"|f.bal.design ==  "N") {
        print("Type-III sum of squares Not yet implemented.") #W <<###########################################
        #W
        #W  INSERIR MAIS UM teste: se o desenho experimental não for balanceado (mesmo número de análises em cada grupo), fazer o Type-III sums of squares (http://www.sthda.com/english/wiki/two-way-anova-test-in-r#compute-two-way-anova-test)
        #W
      } else {
        prg <- f.data #recebe os dados normalizados por alguma das funções de normalização (ou um dos dados direto do Progenesis)
        factor1 <- colnames(prg)[2] #recupera o nome do fator 1 e do fator 2
        factor2 <- colnames(prg)[3]
        nrep <- nrow(prg) # considera como número de repetições o número de linhas da análise atual
        nreptest <- (f.numlevels[[1]])*(f.numlevels[[2]])*f.numreplic # testa se esse número de linhas é coerente com replictas * condições de cada fator
        if (nrep != f.numreplic*f.numlevels[[1]]*f.numlevels[[2]]){
          print("Check the experimental design. The number of rows in the file does not match the number of replicates * number of levels in each factor")
        }
        accs <- colnames(prg)[-(1:3)] # recupera o código de acesso das proteínas (nomes das colunas, exceto as 3 primeiras)
        aovs <- lapply(accs, function(prot_acc) {  #cada termo de accs vai ser tratado como prot_acc neste grupo lapply
          formANOVA2 <- formula(paste(prot_acc, " ~ ", factor1," * ", factor2)) # monta o "texto" da fórmula do 2-way ANOVA
          res.aov2 <- aov(formANOVA2, data = prg) # aplica o teste 2-way ANOVA ao dataframe prg
          pvalues.res.aov2 <- summary(res.aov2)[[1]][5] #extrai a parte do summary que tem os p-valores
          pvalues.res.aov2 <- pvalues.res.aov2[-c(4), , drop = FALSE] #tira a coluna excedente, para ficar só com os p-valores
          names(pvalues.res.aov2)[1] <- prot_acc # coloca o respectivo nome da proteína no nome da coluna referente ao p-valor
          if(f.numlevels[[1]]>2){  # calcula Tukey HSD para o fator1
            pvalues.res.aov2 <- rownames_to_column(pvalues.res.aov2)
            res.Tukey <- TukeyHSD(res.aov2,which = factor1)
            pval.res.Tukey <- data.frame(res.Tukey [[1]][,4])
            names(pval.res.Tukey)[1] <- prot_acc # coloca o respectivo nome da proteína no nome da coluna referente ao p-valor
            pval.res.Tukey <- rownames_to_column(pval.res.Tukey)
            if(pvalues.res.aov2[2,1]<f.tst_thrsh){  # Os p-valores do teste Tukey HSD só são incorporados à tabela de resultados se o ANOVA fro significativo, já que o Tukey é um pós-teste. Caso contrário, é colocado o valor 999
              #just continue
            }else{
              pval.res.Tukey [,2] <- 999
            }
            pvalues.res.aov2 <- bind_rows(pvalues.res.aov2, pval.res.Tukey)
            pvalues.res.aov2 <- column_to_rownames(pvalues.res.aov2)
          }
          if(f.numlevels[[2]]>2){  # calcula Tukey HSD para o fator2
            pvalues.res.aov2 <- rownames_to_column(pvalues.res.aov2)
            res.Tukey <- TukeyHSD(res.aov2,which = factor2)
            pval.res.Tukey <- data.frame(res.Tukey [[1]][,4])
            names(pval.res.Tukey)[1] <- prot_acc # coloca o respectivo nome da proteína no nome da coluna referente ao p-valor
            pval.res.Tukey <- rownames_to_column(pval.res.Tukey)
            if(pvalues.res.aov2[2,2]<f.tst_thrsh){  # Os p-valores do teste Tukey HSD só são incorporados à tabela de resultados se o ANOVA fro significativo, já que o Tukey é um pós-teste. Caso contrário, é colocado o valor 999
              #just continue
            }else{
              pval.res.Tukey [,2] <- 999
            }
            pvalues.res.aov2 <- bind_rows(pvalues.res.aov2, pval.res.Tukey)
            pvalues.res.aov2 <- column_to_rownames(pvalues.res.aov2)
          }
          pvalues.res.aov2 # coloca isso no aovs
        })
        pvals <- do.call(cbind, aovs) %>% add_rownames("subdiv_test") #ajusta a estrutura do aovs para um dataframe
        part.NormStat <- bind_rows(f.empty.NormStat,pvals) #coloca  as linhas dos p-valores na estrutura do f.NormStat - o dataframe geral das comparações deste script, mas ainda um arquivo parcial
        part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]<f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor menor que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
        part.NormStat$normaliz_meth <- f.norm.meth #coloca o nome da estratégia de normalização em todas as linhas desta parcial
        #part.NormStat$stat_test  <- "2w_AOV" # coloca o nome do teste em todas as linhas desta parcial
        part.NormStat$stat_test <- as.character(part.NormStat$stat_test)
        part.NormStat$stat_test [1:3] <- "2w_AOV" # coloca o nome do teste nas 3 primeiras linha desta parcial (que tem o p-valor do ANOVA de cada fator e da interação)
        if(f.numlevels[[1]]>2 | f.numlevels[[2]]>2){  # se tiver dados de Tukey, inclui na part.NormStat
          part.NormStat$stat_test [4:nrow(part.NormStat)] <- "TukeyHSD" # coloca o nome do teste nas outras linhas desta parcial (que tem os p-valores do Tukey HSD)
        }
        part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
        score.sum.par <- sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes para o método de normalização que estiver sendo avaliado.
        f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
        part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
        #W para 2 fatores
        #W Calcula teste de Shapiro Wilk (normalidade) p/ cada ptn --- If p-value > signif threshold, then normality is not likely to be violated.
        Shapiro <- lapply(accs, function(prot_acc) {
          formANOVA2 <- formula(paste(prot_acc, " ~ ", factor1," * ", factor2))
          res.aov2 <- aov(formANOVA2, data = prg)
          aov_residuals <- residuals(object = res.aov2)
          res.shap <- shapiro.test(x=aov_residuals)
          pvalues.res.shap <- data.frame(res.shap[[2]][1])
          names(pvalues.res.shap)[1] <- prot_acc
          pvalues.res.shap
        })
        Spvals <- do.call(cbind, Shapiro) %>% add_rownames("stat_test")
        Spvals[1,1] <-  "Shapiro-Wilk test" 
        part.NormStat <- bind_rows(f.empty.NormStat,Spvals) #coloca  as linhas dos p-valores na estrutura do f.NormStat - o dataframe geral das comparações deste script, mas ainda um arquivo parcial
        part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]>f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor maior que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
        part.NormStat$normaliz_meth <- f.norm.meth #W Achar um jeito de mandar o nome do método automaticamente para a função #coloca o nome da estratégia de normalização em todas as linhas desta parcial
        part.NormStat$stat_test  <- "Shap-Wilk" # coloca o nome do teste em todas as linhas desta parcial
        part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
        score.sum.par <- score.sum.par + sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes para o método de normalização que estiver sendo avaliado.
        f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
        part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
        #W para 2 fatores
        #W Calcula teste de Levene (homogeneidade de variâncias) p/ cada ptn  --- If p-value > signif threshold, then variances are likely to be the same across the conditions.
        # This test asumes normality, so the normality (Shapiro-Wilk) test should be performed before this test.
        Levs <- lapply(accs, function(prot_acc) {
          formLevene <- formula(paste(prot_acc, " ~ ", factor1," * ", factor2))
          res.Levene <- leveneTest(formLevene, data = prg)
          pvalues.res.Levene <- data.frame(res.Levene[[3]][1])
          names(pvalues.res.Levene)[1] <- prot_acc
          pvalues.res.Levene
        })
        Lpvals <- do.call(cbind, Levs) %>% add_rownames("stat_test")
        Lpvals[1,1] <-  "Levene test" 
        part.NormStat <- bind_rows(f.empty.NormStat,Lpvals) #coloca  as linhas dos p-valores na estrutura do f.NormStat - o dataframe geral das comparações deste script, mas ainda um arquivo parcial
        part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]>f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor maior que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
        part.NormStat$normaliz_meth <- f.norm.meth #W Achar um jeito de mandar o nome do método automaticamente para a função #coloca o nome da estratégia de normalização em todas as linhas desta parcial
        part.NormStat$stat_test  <- "Levene" # coloca o nome do teste em todas as linhas desta parcial
        part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
        score.sum.par <- score.sum.par + sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes para o método de normalização que estiver sendo avaliado.
        part.NormStat$score_sum_par[1] <- score.sum.par
        f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
        part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
        
        #W
        #W ART Align-and-rank data for a nonparametric ANOVA (equivalente ao não paramétrico Kruskal Wallis) para ser usado se a distribuição (dos "residuals") não for normal e tiver 2 fatores (http://depts.washington.edu/acelab/proj/art/index.html). Ver também: https://stats.stackexchange.com/questions/12151/is-there-an-equivalent-to-kruskal-wallis-one-way-test-for-a-two-way-model (outras possibilidades)
        #W      Implementadas duas versões: medidas isoladas e medidas repetidas 
        #W      Implementado post-hoc emmeans
        #W

        prgART <- f.data #recebe os dados normalizados por alguma das funções de normalização (ou um dos dados direto do Progenesis)
        prgART[1] <-factor(prgART[1]) # transforma rep_id em fator - necessário para análises de medidas repetidas
        prgART[2] <-factor(prgART[2]) # transforma factor 1 e factor 2 em fator - necessário para criar o modelo ART
        prgART[3] <-factor(prgART[3])
        factor1 <- colnames(prgART)[2] #recupera o nome do fator 1 e do fator 2
        factor2 <- colnames(prgART)[3]
        accs <- colnames(prgART)[-(1:3)] # recupera o código de acesso das proteínas (nomes das colunas, exceto as 3 primeiras)
        
        
        if (f.rep.meas == "Y" | f.rep.meas == "y") {
          
          # teste para dados não paramétricos (sem distribuição normal) ***E com medidas repetidas*** (várias análises do mesmo indivíduo) (é obrigatório ter rep_id repetidos)  
          ARTs.rep.meas <- lapply(accs, function(prot_acc) {  #cada termo de accs vai ser tratado como prot_acc neste grupo lapply
            formART <- formula(paste(prot_acc, " ~ ", factor1," * ", factor2, "+ (1|",IDs,")")) # monta o "texto" da fórmula do modelo (m) para o ART
            m.ART <- art(formART, data = prgART) # cria um modelo para ser usado no ANOVA
            anm <- anova(m.ART) # aplica o teste 2-way ANOVA ao modelo
            rownames(anm) <- c() # apaga os nomes de linhas (eram só números mas impediam a operação)
            anm <- data.frame(anm[-c(2:6)]) # mantém só as colunas que descrevem as condições e de p-valores, apaga as outras
            names(anm)[1] <- "rowname" #a coluna que descreve as condições é preparada para ser colocada como nomes das linahs
            anm <- column_to_rownames(anm) # a descrição das condições passa para nome das linhas
            names(anm)[1] <- prot_acc # o nome da proteína é usado como nome da coluna de p-valores
            if(f.numlevels[[1]]>2){  # calcula emmeans post-hoc para o fator1 se tiver mais que 2 níveis
              anm <- rownames_to_column(anm) # retira novamente o nome das linhas e guarda numa coluna
              emm.TXT <- paste("emmeans(artlm(m.ART, \"",factor1,"\"), pairwise ~ ",factor1,")", sep="") #prepara o comando do emmeans em forma de texto
              res.emm <- eval(parse(text=emm.TXT)) # converte o texto em comando e executa
              res.emm <- data.frame(res.emm$contrasts) # cria um dataframe com o subset que contém p-valores e descrição
              res.emm <- res.emm[-c(2:5)] # deixa só as colunas de p-valores e descrição das condições
              names(res.emm) <- c("rowname", prot_acc) # coloca o respectivo nome da proteína no nome da coluna referente ao p-valor
              if(anm[2,1]<f.tst_thrsh){  # Os p-valores do teste emmeans só são incorporados à tabela de resultados se o ANOVA for significativo, já que o emmeans é um pós-teste. Caso contrário, é colocado o valor 999
                #just continue
              }else{
                res.emm [,2] <- 999
              }
              anm <- bind_rows(anm, res.emm) # acrescenta as linhas de resultado do emmeans ao resultado do ART anova
              anm <- column_to_rownames(anm) # devolve o nome das linhas
            }
            if(f.numlevels[[2]]>2){  # calcula emmeans post-hoc para o fator2, se tiver mais que dois níveis (repete o que fez para o fator 1)
              anm <- rownames_to_column(anm)
              emm.TXT <- paste("emmeans(artlm(m.ART, \"",factor2,"\"), pairwise ~ ",factor2,")", sep="")
              res.emm <- eval(parse(text=emm.TXT))
              res.emm <- data.frame(res.emm$contrasts) 
              res.emm <- res.emm[-c(2:5)] 
              names(res.emm) <- c("rowname", prot_acc) 
              if(anm[2,2]<f.tst_thrsh){  
                #just continue
              }else{
                res.emm [,2] <- 999
              }
              anm <- bind_rows(anm, res.emm)
              anm <- column_to_rownames(anm)
            }
            anm # transfere o valor de anm para ARTS.indep e finaliza o loop do lapply
          })
          pvals <- do.call(cbind, ARTs.indep) %>% add_rownames("subdiv_test") #ajusta a estrutura do ARTS.indep para um dataframe
          part.NormStat <- bind_rows(f.empty.NormStat,pvals) #coloca  as linhas dos p-valores num arquivo parcial que tem a mesma estrutura do f.NormStat - o dataframe geral das comparações deste script. 
          part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]<f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor menor que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
          part.NormStat$normaliz_meth <- f.norm.meth #coloca o nome da estratégia de normalização em todas as linhas desta parcial
          part.NormStat$stat_test <- as.character(part.NormStat$stat_test)
          part.NormStat$stat_test [1:3] <- "npar_2w_AOV" # coloca o nome do teste nas 3 primeiras linhas desta parcial (que tem o p-valor do ANOVA de cada fator e da interação)
          if(f.numlevels[[1]]>2 | f.numlevels[[2]]>2){  # se tiver dados de post-hoc (emmeans), inclui na part.NormStat
            part.NormStat$stat_test [4:nrow(part.NormStat)] <- "emmeans" # coloca o nome do teste nas outras linhas desta parcial (que tem os p-valores do emmeans)
          }
          part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
          score.sum.npar <- sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes não paramétricos para o método de normalização que estiver sendo avaliado.
          f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
          part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
          
        }else{
          
          # teste para dados não paramétricos (sem distribuição normal) ***E sem medidas repetidas*** (cada análise é independente) (não pode ter rep_id repetidos)  
          ARTs.indep <- lapply(accs, function(prot_acc) {  #cada termo de accs vai ser tratado como prot_acc neste grupo lapply
            formART <- formula(paste(prot_acc, " ~ ", factor1," * ", factor2)) # monta o "texto" da fórmula do modelo (m) para o ART
            m.ART <- art(formART, data = prgART) # cria um modelo para ser usado no ANOVA
            anm <- anova(m.ART) # aplica o teste 2-way ANOVA ao modelo
            rownames(anm) <- c() # apaga os nomes de linhas (eram só números mas impediam a operação)
            anm <- data.frame(anm[-c(2:6)]) # mantém só as colunas que descrevem as condições e de p-valores, apaga as outras
            names(anm)[1] <- "rowname" #a coluna que descreve as condições é preparada para ser colocada como nomes das linahs
            anm <- column_to_rownames(anm) # a descrição das condições passa para nome das linhas
            names(anm)[1] <- prot_acc # o nome da proteína é usado como nome da coluna de p-valores
            if(f.numlevels[[1]]>2){  # calcula emmeans post-hoc para o fator1 se tiver mais que 2 níveis
              anm <- rownames_to_column(anm) # retira novamente o nome das linhas e guarda numa coluna
              emm.TXT <- paste("emmeans(artlm(m.ART, \"",factor1,"\"), pairwise ~ ",factor1,")", sep="") #prepara o comando do emmeans em forma de texto
              res.emm <- eval(parse(text=emm.TXT)) # converte o texto em comando e executa
              res.emm <- data.frame(res.emm$contrasts) # cria um dataframe com o subset que contém p-valores e descrição
              res.emm <- res.emm[-c(2:5)] # deixa só as colunas de p-valores e descrição das condições
              names(res.emm) <- c("rowname", prot_acc) # coloca o respectivo nome da proteína no nome da coluna referente ao p-valor
              if(anm[2,1]<f.tst_thrsh){  # Os p-valores do teste emmeans só são incorporados à tabela de resultados se o ANOVA for significativo, já que o emmeans é um pós-teste. Caso contrário, é colocado o valor 999
                #just continue
              }else{
                res.emm [,2] <- 999
              }
              anm <- bind_rows(anm, res.emm) # acrescenta as linhas de resultado do emmeans ao resultado do ART anova
              anm <- column_to_rownames(anm) # devolve o nome das linhas
            }
            if(f.numlevels[[2]]>2){  # calcula emmeans post-hoc para o fator2, se tiver mais que dois níveis (repete o que fez para o fator 1)
              anm <- rownames_to_column(anm)
              emm.TXT <- paste("emmeans(artlm(m.ART, \"",factor2,"\"), pairwise ~ ",factor2,")", sep="")
              res.emm <- eval(parse(text=emm.TXT))
              res.emm <- data.frame(res.emm$contrasts) 
              res.emm <- res.emm[-c(2:5)] 
              names(res.emm) <- c("rowname", prot_acc) 
              if(anm[2,2]<f.tst_thrsh){  
                #just continue
              }else{
                res.emm [,2] <- 999
              }
              anm <- bind_rows(anm, res.emm)
              anm <- column_to_rownames(anm)
            }
            anm # transfere o valor de anm para ARTS.indep e finaliza o loop do lapply
          })
          pvals <- do.call(cbind, ARTs.indep) %>% add_rownames("subdiv_test") #ajusta a estrutura do ARTS.indep para um dataframe
          part.NormStat <- bind_rows(f.empty.NormStat,pvals) #coloca  as linhas dos p-valores num arquivo parcial que tem a mesma estrutura do f.NormStat - o dataframe geral das comparações deste script. 
          part.NormStat$number_of_ptns_pass <- rowSums(part.NormStat[,7:ncol(part.NormStat)]<f.tst_thrsh,na.rm = TRUE) # conta o número de linhas com p-valor menor que o limite e coloca na respectiva coluna (conta só nas colunas de proteínas, a partir da 6a coluna)
          part.NormStat$normaliz_meth <- f.norm.meth #coloca o nome da estratégia de normalização em todas as linhas desta parcial
          part.NormStat$stat_test <- as.character(part.NormStat$stat_test)
          part.NormStat$stat_test [1:3] <- "npar_2w_AOV_rep" # coloca o nome do teste nas 3 primeiras linhas desta parcial (que tem o p-valor do ANOVA de cada fator e da interação)
          if(f.numlevels[[1]]>2 | f.numlevels[[2]]>2){  # se tiver dados de post-hoc (emmeans), inclui na part.NormStat
            part.NormStat$stat_test [4:nrow(part.NormStat)] <- "emmeans" # coloca o nome do teste nas outras linhas desta parcial (que tem os p-valores do emmeans)
          }
          part.NormStat$test_threshold <- f.tst_thrsh # coloca o valor de limite do teste em todas as linhas desta parcial
          score.sum.npar <- sum(part.NormStat$number_of_ptns_pass) # cria um score com a soma de todos os testes não paramétricos para o método de normalização que estiver sendo avaliado.
          f.NormStat <- bind_rows(f.NormStat,part.NormStat) # acrescenta as linhas desta parcial ao dataframe global
          part.NormStat <- part.NormStat[-c(1:nrow(part.NormStat)),] # deleta todas as linhas do dataframe parcial, restaura a condição original dele para poder receber dados de outros testes
          
        }
        
      }
    } else {
      #factorial anova (separar three-way?)
      print("Not yet implemented.")
    }
  }
  
  
  return(f.NormStat)
}
# Final da função de estatística comparativa 
#--------------------------------------------------------------------------------------------------------------        






############################
# carregando dados e fazendo log


# carregando data_raw
raw_data <- read.csv(dlg_open(default = getwd(), title='Choose the Progenesis RAW data', multiple = FALSE, filters = dlg_filters["All", ])$res, as.is=1, dec=".")
#raw_data <- read.csv(file.choose(), as.is=1, dec=".")

# carregando data_prog
Progenesis_data <- read.csv(dlg_open(default = getwd(), title='Choose the Progenesis NORMALIZED data', multiple = FALSE, filters = dlg_filters["All", ])$res, as.is=1, dec=".")
#Progenesis_data <- read.csv(file.choose(), as.is=1)




#W
#W Parâmetros de definição desenho experimental (fatores, níveis, replicatas, balanceamento, limite de significância)
#W
num.factors <- as.numeric(dlgInput("Enter number of factors", default=2, Sys.info()["user"])$res)
num.levels <- c()
for (i in 1:num.factors) {
  num.levels[i] <- as.numeric(dlgInput(paste("Enter number of levels for factor ",i,"(",colnames(raw_data)[i+1],"): "), default=2, Sys.info()["user"])$res)
}
#if (num.factors == 1 && num.levels == 2) {  #W isso provavelmente precisará ser mudado para incluir a opção de ANOVA de medidas repetidas (repeated measures ANOVA) p/ 1-way, 2-way e 3-way (https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/)
#  paired <- dlgInput("Are the samples paired? (e.g. two conditions observed in the same individual) (y/n)", default="y", Sys.info()["user"])$res
#}else{
#  paired <- NULL
#}  
rep.meas <- dlgInput("Is the experimental design based on repeated measures? (y/n)", default="n", Sys.info()["user"])$res
bal.design <- dlgInput("Is this experimental design balanced? (same number of replicates for all conditions) (y/n)", default="y", Sys.info()["user"])$res
num.replic <- as.numeric(dlgInput("Enter number of replicates", default=8, Sys.info()["user"])$res)
tst_thrsh <- as.numeric(dlgInput("Enter the threshold for statistical significance:", default=0.05, Sys.info()["user"])$res)
imputation <- dlgInput("Perform imputation? (y/n)", default="y", Sys.info()["user"])$res


if(imputation=="y"|imputation=="Y") {
  raw_data <- imputed.tab(raw_data, num.factors)
  Progenesis_data <- imputed.tab(Progenesis_data, num.factors)
}


# fazendo dados raw com log base 2
dados_log2 <- log(raw_data[,(num.factors+2):ncol(raw_data)], 2)
raw_data_log <- cbind(raw_data[,1:3], dados_log2)

# fazendo dados progenesis com log base 2
dados_log <- log(Progenesis_data[,(num.factors+2):ncol(Progenesis_data)], 2)
Progenesis_data_log <- cbind(Progenesis_data[,1:3], dados_log)





#
#W inicializa o dataframe de estatística comparativa entre as normalizações (para ser usado na função Compar)
#
#
empty.NormStat <- data.frame(normaliz_meth=character(1),stat_test=character(1),subdiv_test=character(1),test_threshold=numeric(1),number_of_ptns_pass=numeric(1),score_sum_par=numeric(1))
empty.NormStat <- data.frame(empty.NormStat, raw_data[1,4:ncol(raw_data)]) [-1,]
#W --- ATENÇÃO --- ISTO ZERA O DATAFRAME NormStat, portanto só deve ser usado uma vez para criar. NÃO USAR ENTRE UMA NORMALIZAÇÃO E OUTRA, POIS VAI DEIXAR SÓ OS DADOS DA ÚLTIMA
NormStat <- empty.NormStat



#W
#W Nomenclatura dos métodos de normalização
#W
#W 1- Raw ou Progenesis (R ou P) - método de saída do Progenesis, Raw = sem normalização, Progenesis = normalizado pelo Progenesis
#W 2- Log ou Não-log    (L ou N) - transformação log feita aqui
#W 3- Divisão, Subtração, Radiciação, Log, Z-score, Pareto, Min-max, Média-min-max  (D, S, R, L, ZS, PP, MM ou ME) - cálculo feito entre o valor original e o valor do grupo 
#W 4- Mediana, Média, Soma (M, A ou S) - cálculo feito entre os valores do grupo - não se aplica a Z-score, Pareto, Min-max, Média-min-max, pois já incluem o cálculo dentro do grupo. Por isso  Z-score, Pareto, Min-max, Média-min-max já têm 2 caracteres na nomenclatura, correnpondente ao 3o e 4o caractere.
#W 5- All ou 1 (A ou 1) - Agrupamento de replicatas - O grupo inclui todas as replicatas ou só a replicata do valor original
#W 6- All ou 1 (A ou 1) - Agrupamento de condições - O grupo inclui todas as condições ou só a condição do valor original
#W 7- All ou 1 (A ou 1) - Agrupamento de proteínas - O grupo inclui todas as proteínas ou só a proteína do valor original
#W
#W EX: RNDMA11
#W     |||||||
#W     ||||||-o grupo inclui SÓ a proteína do valor original
#W     |||||-o grupo inclui SÓ a condição do valor original
#W     ||||-o grupo inclui TODAS as replicatas
#W     |||-será calculada a MEDIANA do grupo
#W     ||-cada valor original será substituido por esse valor original DIVIDIDO pela mediana do grupo
#W     |-não é usada a transformação LOG dos valores originais
#W     -não é usada a normalização inicial do PROGENESIS
#W     
#W     
#W     Os dois primeiros caracteres se referem a transformações de todos os valores antes da normalização neste script
#W     Os cinco caracteres seguintes descrevem o método usado neste script
#W     Os métodos que têm "-" nesses cinco caracteres não passam por normalização neste script. Recebem, no máximo, a transformação log.
#W     


#######################################
# no norm (controls)
#######################################

#W raw no log - no norm (RN.....)
#W Não usa função de normalização , pois este é o controle1, raw do Progenesis, sem normalizar e sem log. Faz só a estatiística para comparar
RN..... <- raw_data
NormStat <- funcaoCompar(raw_data, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RN.....") # ver o que mais precisa quando for implementar a função.

#W raw log - no norm (RL.....)
#W Não usa função de normalização , pois este é o controle2, raw do Progenesis, sem normalizar, mas com log. Faz só a estatiística para comparar
RL..... <- raw_data_log
NormStat <- funcaoCompar(raw_data_log, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RL.....") # ver o que mais precisa quando for implementar a função.

#W Progenesis no log - no norm (PN.....)
#W Não usa função de normalização , pois este é o controle3, normalizado pelo Progenesis, sem log. Faz só a estatiística para comparar
PN..... <- Progenesis_data
NormStat <- funcaoCompar(Progenesis_data, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PN.....") # ver o que mais precisa quando for implementar a função.

#W Progenesis log - no norm (PL.....)
#W Não usa função de normalização , pois este é o controle4, normalizado pelo Progenesis e com log. Faz só a estatiística para comparar
PL..... <- Progenesis_data_log
NormStat <- funcaoCompar(Progenesis_data_log, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PL.....") # ver o que mais precisa quando for implementar a função.





#######################################
# 1A1
#######################################

# raw no log
RNDM1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, div, median, NULL, NULL) 
NormStat <- funcaoCompar(RNDM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDM1A1") 
RNDA1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RNDA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDA1A1")
RNDS1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RNDS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDS1A1")
RNSM1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RNSM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSM1A1")
RNSA1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RNSA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSA1A1")
RNSS1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RNSS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSS1A1")
RNRM1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RNRM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRM1A1")
RNRA1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RNRA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRA1A1")
RNRS1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RNRS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRS1A1")
RNLM1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RNLM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLM1A1")
RNLA1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RNLA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLA1A1")
RNLS1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RNLS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLS1A1")
RNZS1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RNZS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNZS1A1")
RNPP1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RNPP1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNPP1A1")
RNMM1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RNMM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMM1A1")
RNME1A1 <- funcao1A1(raw_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RNME1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNME1A1")

# raw log
RLDM1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(RLDM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDM1A1")
RLDA1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RLDA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDA1A1")
RLDS1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RLDS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDS1A1")
RLSM1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RLSM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSM1A1")
RLSA1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RLSA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSA1A1")
RLSS1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RLSS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSS1A1")
RLRM1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RLRM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRM1A1")
RLRA1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RLRA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRA1A1")
RLRS1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RLRS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRS1A1")
RLLM1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RLLM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLM1A1")
RLLA1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RLLA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLA1A1")
RLLS1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RLLS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLS1A1")
RLZS1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RLZS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLZS1A1")
RLPP1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RLPP1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLPP1A1")
RLMM1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RLMM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMM1A1")
RLME1A1 <- funcao1A1(raw_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RLME1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLME1A1")

# prog no log
PNDM1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PNDM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDM1A1")
PNDA1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PNDA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDA1A1")
PNDS1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PNDS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDS1A1")
PNSM1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PNSM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSM1A1")
PNSA1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PNSA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSA1A1")
PNSS1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PNSS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSS1A1")
PNRM1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PNRM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRM1A1")
PNRA1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PNRA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRA1A1")
PNRS1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PNRS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRS1A1")
PNLM1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PNLM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLM1A1")
PNLA1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PNLA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLA1A1")
PNLS1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PNLS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLS1A1")
PNZS1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PNZS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNZS1A1")
PNPP1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PNPP1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNPP1A1")
PNMM1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PNMM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMM1A1")
PNME1A1 <- funcao1A1(Progenesis_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PNME1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNME1A1")

# prog log
PLDM1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PLDM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDM1A1")
PLDA1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PLDA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDA1A1")
PLDS1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PLDS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDS1A1")
PLSM1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PLSM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSM1A1")
PLSA1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PLSA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSA1A1")
PLSS1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PLSS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSS1A1")
PLRM1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PLRM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRM1A1")
PLRA1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PLRA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRA1A1")
PLRS1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PLRS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRS1A1")
PLLM1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PLLM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLM1A1")
PLLA1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PLLA1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLA1A1")
PLLS1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PLLS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLS1A1")
PLZS1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PLZS1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLZS1A1")
PLPP1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PLPP1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLPP1A1")
PLMM1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PLMM1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMM1A1")
PLME1A1 <- funcao1A1(Progenesis_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PLME1A1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLME1A1")


#alguns tem warnings pois tem algumas proteinas que todas as observacoes sao 0
#isso faz com que algumas contas sejam divisoes por 0, que causa resultados NaN (not a number)
#ou log(0) que da -Inf





#######################################
# A11
#######################################

# raw no log
RNDMA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, div, median, NULL, NULL) #W #################### W -> CORRIGIDA para funcionar com 1, 2, 3 ou mais fatores e trazendo o NULL para cá, em vez de deixar na função
NormStat <- funcaoCompar(RNDMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDMA11") # chama a função para incluir este método de normalização na tabela de estatística comparativa entre os métodos.
RNDAA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RNDAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDAA11")
RNDSA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RNDSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDSA11")
RNSMA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RNSMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSMA11")
RNSAA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RNSAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSAA11")
RNSSA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RNSSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSSA11")
RNRMA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RNRMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRMA11")
RNRAA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RNRAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRAA11")
RNRSA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RNRSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRSA11")
RNLMA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RNLMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLMA11")
RNLAA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RNLAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLAA11")
RNLSA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RNLSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLSA11")
RNZSA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RNZSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNZSA11")
RNPPA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RNPPA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNPPA11")
RNMMA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RNMMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMMA11")
RNMEA11 <- funcaoA11(raw_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RNMEA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMEA11")

# raw log
RLDMA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(RLDMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDMA11")
RLDAA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RLDAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDAA11")
RLDSA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RLDSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDSA11")
RLSMA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RLSMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSMA11")
RLSAA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RLSAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSAA11")
RLSSA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RLSSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSSA11")
RLRMA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RLRMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRMA11")
RLRAA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RLRAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRAA11")
RLRSA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RLRSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRSA11")
RLLMA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RLLMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLMA11")
RLLAA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RLLAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLAA11")
RLLSA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RLLSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLSA11")
RLZSA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RLZSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLZSA11")
RLPPA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RLPPA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLPPA11")
RLMMA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RLMMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMMA11")
RLMEA11 <- funcaoA11(raw_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RLMEA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMEA11")

# prog no log
PNDMA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PNDMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDMA11")
PNDAA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PNDAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDAA11")
PNDSA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PNDSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDSA11")
PNSMA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PNSMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSMA11")
PNSAA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PNSAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSAA11")
PNSSA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PNSSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSSA11")
PNRMA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PNRMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRMA11")
PNRAA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PNRAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRAA11")
PNRSA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PNRSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRSA11")
PNLMA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PNLMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLMA11")
PNLAA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PNLAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLAA11")
PNLSA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PNLSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLSA11")
PNZSA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PNZSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNZSA11")
PNPPA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PNPPA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNPPA11")
PNMMA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PNMMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMMA11")
PNMEA11 <- funcaoA11(Progenesis_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PNMEA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMEA11")

# prog log
PLDMA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PLDMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDMA11")
PLDAA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PLDAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDAA11")
PLDSA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PLDSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDSA11")
PLSMA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PLSMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSMA11")
PLSAA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PLSAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSAA11")
PLSSA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PLSSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSSA11")
PLRMA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PLRMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRMA11")
PLRAA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PLRAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRAA11")
PLRSA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PLRSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRSA11")
PLLMA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PLLMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLMA11")
PLLAA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PLLAA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLAA11")
PLLSA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PLLSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLSA11")
PLZSA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PLZSA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLZSA11")
PLPPA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PLPPA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLPPA11")
PLMMA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PLMMA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMMA11")
PLMEA11 <- funcaoA11(Progenesis_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PLMEA11, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMEA11")


#alguns tem warnings pois tem algumas proteinas que todas as observacoes sao 0
#isso faz com que algumas contas sejam divisoes por 0, que causa resultados NaN (not a number)
#ou log(0) que da -Inf






#######################################
# 11A
#######################################

# raw no log
RNDM11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, div, median, NULL, NULL) 
NormStat <- funcaoCompar(RNDM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDM11A") 
RNDA11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RNDA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDA11A")
RNDS11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RNDS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDS11A")
RNSM11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RNSM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSM11A")
RNSA11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RNSA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSA11A")
RNSS11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RNSS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSS11A")
RNRM11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RNRM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRM11A")
RNRA11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RNRA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRA11A")
RNRS11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RNRS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRS11A")
RNLM11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RNLM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLM11A")
RNLA11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RNLA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLA11A")
RNLS11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RNLS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLS11A")
RNZS11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RNZS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNZS11A")
RNPP11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RNPP11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNPP11A")
RNMM11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RNMM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMM11A")
RNME11A <- funcao11A(raw_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RNME11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNME11A")

# raw log
RLDM11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(RLDM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDM11A")
RLDA11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RLDA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDA11A")
RLDS11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RLDS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDS11A")
RLSM11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RLSM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSM11A")
RLSA11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RLSA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSA11A")
RLSS11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RLSS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSS11A")
RLRM11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RLRM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRM11A")
RLRA11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RLRA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRA11A")
RLRS11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RLRS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRS11A")
RLLM11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RLLM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLM11A")
RLLA11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RLLA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLA11A")
RLLS11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RLLS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLS11A")
RLZS11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RLZS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLZS11A")
RLPP11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RLPP11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLPP11A")
RLMM11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RLMM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMM11A")
RLME11A <- funcao11A(raw_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RLME11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLME11A")

# prog no log
PNDM11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PNDM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDM11A")
PNDA11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PNDA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDA11A")
PNDS11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PNDS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDS11A")
PNSM11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PNSM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSM11A")
PNSA11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PNSA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSA11A")
PNSS11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PNSS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSS11A")
PNRM11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PNRM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRM11A")
PNRA11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PNRA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRA11A")
PNRS11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PNRS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRS11A")
PNLM11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PNLM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLM11A")
PNLA11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PNLA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLA11A")
PNLS11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PNLS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLS11A")
PNZS11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PNZS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNZS11A")
PNPP11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PNPP11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNPP11A")
PNMM11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PNMM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMM11A")
PNME11A <- funcao11A(Progenesis_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PNME11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNME11A")

# prog log
PLDM11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PLDM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDM11A")
PLDA11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PLDA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDA11A")
PLDS11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PLDS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDS11A")
PLSM11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PLSM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSM11A")
PLSA11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PLSA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSA11A")
PLSS11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PLSS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSS11A")
PLRM11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PLRM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRM11A")
PLRA11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PLRA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRA11A")
PLRS11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PLRS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRS11A")
PLLM11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PLLM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLM11A")
PLLA11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PLLA11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLA11A")
PLLS11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PLLS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLS11A")
PLZS11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PLZS11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLZS11A")
PLPP11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PLPP11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLPP11A")
PLMM11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PLMM11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMM11A")
PLME11A <- funcao11A(Progenesis_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PLME11A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLME11A")


#alguns tem warnings pois tem algumas proteinas que todas as observacoes sao 0
#isso faz com que algumas contas sejam divisoes por 0, que causa resultados NaN (not a number)
#ou log(0) que da -Inf





#######################################
# 1AA
#######################################

# raw no log
RNDM1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, div, median, NULL, NULL) 
NormStat <- funcaoCompar(RNDM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDM1AA") 
RNDA1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RNDA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDA1AA")
RNDS1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RNDS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDS1AA")
RNSM1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RNSM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSM1AA")
RNSA1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RNSA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSA1AA")
RNSS1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RNSS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSS1AA")
RNRM1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RNRM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRM1AA")
RNRA1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RNRA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRA1AA")
RNRS1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RNRS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRS1AA")
RNLM1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RNLM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLM1AA")
RNLA1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RNLA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLA1AA")
RNLS1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RNLS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLS1AA")
RNZS1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RNZS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNZS1AA")
RNPP1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RNPP1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNPP1AA")
RNMM1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RNMM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMM1AA")
RNME1AA <- funcao1AA(raw_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RNME1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNME1AA")

# raw log
RLDM1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(RLDM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDM1AA")
RLDA1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RLDA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDA1AA")
RLDS1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RLDS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDS1AA")
RLSM1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RLSM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSM1AA")
RLSA1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RLSA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSA1AA")
RLSS1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RLSS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSS1AA")
RLRM1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RLRM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRM1AA")
RLRA1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RLRA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRA1AA")
RLRS1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RLRS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRS1AA")
RLLM1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RLLM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLM1AA")
RLLA1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RLLA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLA1AA")
RLLS1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RLLS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLS1AA")
RLZS1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RLZS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLZS1AA")
RLPP1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RLPP1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLPP1AA")
RLMM1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RLMM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMM1AA")
RLME1AA <- funcao1AA(raw_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RLME1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLME1AA")

# prog no log
PNDM1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PNDM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDM1AA")
PNDA1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PNDA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDA1AA")
PNDS1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PNDS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDS1AA")
PNSM1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PNSM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSM1AA")
PNSA1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PNSA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSA1AA")
PNSS1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PNSS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSS1AA")
PNRM1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PNRM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRM1AA")
PNRA1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PNRA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRA1AA")
PNRS1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PNRS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRS1AA")
PNLM1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PNLM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLM1AA")
PNLA1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PNLA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLA1AA")
PNLS1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PNLS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLS1AA")
PNZS1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PNZS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNZS1AA")
PNPP1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PNPP1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNPP1AA")
PNMM1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PNMM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMM1AA")
PNME1AA <- funcao1AA(Progenesis_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PNME1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNME1AA")

# prog log
PLDM1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PLDM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDM1AA")
PLDA1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PLDA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDA1AA")
PLDS1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PLDS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDS1AA")
PLSM1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PLSM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSM1AA")
PLSA1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PLSA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSA1AA")
PLSS1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PLSS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSS1AA")
PLRM1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PLRM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRM1AA")
PLRA1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PLRA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRA1AA")
PLRS1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PLRS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRS1AA")
PLLM1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PLLM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLM1AA")
PLLA1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PLLA1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLA1AA")
PLLS1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PLLS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLS1AA")
PLZS1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PLZS1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLZS1AA")
PLPP1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PLPP1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLPP1AA")
PLMM1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PLMM1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMM1AA")
PLME1AA <- funcao1AA(Progenesis_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PLME1AA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLME1AA")


#alguns tem warnings pois tem algumas proteinas que todas as observacoes sao 0
#isso faz com que algumas contas sejam divisoes por 0, que causa resultados NaN (not a number)
#ou log(0) que da -Inf






#######################################
# A1A
#######################################

# raw no log
RNDMA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, div, median, NULL, NULL) 
NormStat <- funcaoCompar(RNDMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDMA1A") 
RNDAA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RNDAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDAA1A")
RNDSA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RNDSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDSA1A")
RNSMA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RNSMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSMA1A")
RNSAA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RNSAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSAA1A")
RNSSA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RNSSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSSA1A")
RNRMA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RNRMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRMA1A")
RNRAA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RNRAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRAA1A")
RNRSA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RNRSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRSA1A")
RNLMA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RNLMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLMA1A")
RNLAA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RNLAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLAA1A")
RNLSA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RNLSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLSA1A")
RNZSA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RNZSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNZSA1A")
RNPPA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RNPPA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNPPA1A")
RNMMA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RNMMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMMA1A")
RNMEA1A <- funcaoA1A(raw_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RNMEA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMEA1A")

# raw log
RLDMA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(RLDMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDMA1A")
RLDAA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RLDAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDAA1A")
RLDSA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RLDSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDSA1A")
RLSMA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RLSMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSMA1A")
RLSAA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RLSAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSAA1A")
RLSSA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RLSSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSSA1A")
RLRMA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RLRMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRMA1A")
RLRAA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RLRAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRAA1A")
RLRSA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RLRSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRSA1A")
RLLMA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RLLMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLMA1A")
RLLAA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RLLAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLAA1A")
RLLSA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RLLSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLSA1A")
RLZSA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RLZSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLZSA1A")
RLPPA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RLPPA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLPPA1A")
RLMMA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RLMMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMMA1A")
RLMEA1A <- funcaoA1A(raw_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RLMEA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMEA1A")

# prog no log
PNDMA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PNDMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDMA1A")
PNDAA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PNDAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDAA1A")
PNDSA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PNDSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDSA1A")
PNSMA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PNSMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSMA1A")
PNSAA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PNSAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSAA1A")
PNSSA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PNSSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSSA1A")
PNRMA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PNRMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRMA1A")
PNRAA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PNRAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRAA1A")
PNRSA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PNRSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRSA1A")
PNLMA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PNLMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLMA1A")
PNLAA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PNLAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLAA1A")
PNLSA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PNLSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLSA1A")
PNZSA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PNZSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNZSA1A")
PNPPA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PNPPA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNPPA1A")
PNMMA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PNMMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMMA1A")
PNMEA1A <- funcaoA1A(Progenesis_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PNMEA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMEA1A")

# prog log
PLDMA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PLDMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDMA1A")
PLDAA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PLDAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDAA1A")
PLDSA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PLDSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDSA1A")
PLSMA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PLSMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSMA1A")
PLSAA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PLSAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSAA1A")
PLSSA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PLSSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSSA1A")
PLRMA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PLRMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRMA1A")
PLRAA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PLRAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRAA1A")
PLRSA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PLRSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRSA1A")
PLLMA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PLLMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLMA1A")
PLLAA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PLLAA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLAA1A")
PLLSA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PLLSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLSA1A")
PLZSA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PLZSA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLZSA1A")
PLPPA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PLPPA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLPPA1A")
PLMMA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PLMMA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMMA1A")
PLMEA1A <- funcaoA1A(Progenesis_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PLMEA1A, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMEA1A")




#######################################
# AA1
#######################################

# raw no log
RNDMAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, div, median, NULL, NULL) #W #################### W -> CORRIGIDA para funcionar com 1, 2, 3 ou mais fatores e trazendo o NULL para cá, em vez de deixar na função
NormStat <- funcaoCompar(RNDMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDMAA1") # chama a função para incluir este método de normalização na tabela de estatística comparativa entre os métodos.
RNDAAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RNDAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDAAA1")
RNDSAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RNDSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDSAA1")
RNSMAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RNSMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSMAA1")
RNSAAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RNSAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSAAA1")
RNSSAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RNSSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSSAA1")
RNRMAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RNRMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRMAA1")
RNRAAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RNRAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRAAA1")
RNRSAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RNRSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRSAA1")
RNLMAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RNLMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLMAA1")
RNLAAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RNLAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLAAA1")
RNLSAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RNLSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLSAA1")
RNZSAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RNZSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNZSAA1")
RNPPAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RNPPAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNPPAA1")
RNMMAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RNMMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMMAA1")
RNMEAA1 <- funcaoAA1(raw_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RNMEAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMEAA1")

# raw log
RLDMAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(RLDMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDMAA1")
RLDAAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RLDAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDAAA1")
RLDSAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RLDSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDSAA1")
RLSMAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RLSMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSMAA1")
RLSAAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RLSAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSAAA1")
RLSSAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RLSSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSSAA1")
RLRMAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(RLRMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRMAA1")
RLRAAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RLRAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRAAA1")
RLRSAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RLRSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRSAA1")
RLLMAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RLLMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLMAA1")
RLLAAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RLLAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLAAA1")
RLLSAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RLLSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLSAA1")
RLZSAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RLZSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLZSAA1")
RLPPAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RLPPAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLPPAA1")
RLMMAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(RLMMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMMAA1")
RLMEAA1 <- funcaoAA1(raw_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(RLMEAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMEAA1")

# prog no log
PNDMAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PNDMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDMAA1")
PNDAAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PNDAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDAAA1")
PNDSAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PNDSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDSAA1")
PNSMAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PNSMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSMAA1")
PNSAAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PNSAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSAAA1")
PNSSAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PNSSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSSAA1")
PNRMAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PNRMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRMAA1")
PNRAAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PNRAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRAAA1")
PNRSAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PNRSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRSAA1")
PNLMAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PNLMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLMAA1")
PNLAAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PNLAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLAAA1")
PNLSAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PNLSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLSAA1")
PNZSAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PNZSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNZSAA1")
PNPPAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PNPPAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNPPAA1")
PNMMAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PNMMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMMAA1")
PNMEAA1 <- funcaoAA1(Progenesis_data, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PNMEAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMEAA1")

# prog log
PLDMAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, div, median, NULL, NULL)
NormStat <- funcaoCompar(PLDMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDMAA1")
PLDAAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PLDAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDAAA1")
PLDSAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PLDSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDSAA1")
PLSMAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PLSMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSMAA1")
PLSAAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PLSAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSAAA1")
PLSSAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PLSSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSSAA1")
PLRMAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, root,median, NULL, NULL)
NormStat <- funcaoCompar(PLRMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRMAA1")
PLRAAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PLRAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRAAA1")
PLRSAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PLRSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRSAA1")
PLLMAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PLLMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLMAA1")
PLLAAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PLLAAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLAAA1")
PLLSAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PLLSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLSAA1")
PLZSAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PLZSAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLZSAA1")
PLPPAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PLPPAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLPPAA1")
PLMMAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, mm,min,max, NULL)
NormStat <- funcaoCompar(PLMMAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMMAA1")
PLMEAA1 <- funcaoAA1(Progenesis_data_log, num.factors, num.levels, num.replic, me,mean,max,min)
NormStat <- funcaoCompar(PLMEAA1, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMEAA1")


#alguns tem warnings pois tem algumas proteinas que todas as observacoes sao 0
#isso faz com que algumas contas sejam divisoes por 0, que causa resultados NaN (not a number)
#ou log(0) que da -Inf






#RNDMAAA <- funcaoAAA(raw_data, num.factors, div, median, NULL, NULL) 


#######################################
# AAA
#######################################

# raw no log
RNDMAAA <- funcaoAAA(raw_data, num.factors, div, median, NULL, NULL) 

NormStat <- funcaoCompar(RNDMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDMAAA") # chama a função para incluir este método de normalização na tabela de estatística comparativa entre os métodos.
RNDAAAA <- funcaoAAA(raw_data, num.factors, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RNDAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDAAAA")
RNDSAAA <- funcaoAAA(raw_data, num.factors, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RNDSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNDSAAA")
RNSMAAA <- funcaoAAA(raw_data, num.factors, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RNSMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSMAAA")
RNSAAAA <- funcaoAAA(raw_data, num.factors, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RNSAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSAAAA")
RNSSAAA <- funcaoAAA(raw_data, num.factors, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RNSSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNSSAAA")
RNRMAAA <- funcaoAAA(raw_data, num.factors, root,median, NULL, NULL)
NormStat <- funcaoCompar(RNRMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRMAAA")
RNRAAAA <- funcaoAAA(raw_data, num.factors, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RNRAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRAAAA")
RNRSAAA <- funcaoAAA(raw_data, num.factors, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RNRSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNRSAAA")
RNLMAAA <- funcaoAAA(raw_data, num.factors, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RNLMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLMAAA")
RNLAAAA <- funcaoAAA(raw_data, num.factors, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RNLAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLAAAA")
RNLSAAA <- funcaoAAA(raw_data, num.factors, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RNLSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNLSAAA")
RNZSAAA <- funcaoAAA(raw_data, num.factors, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RNZSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNZSAAA")
RNPPAAA <- funcaoAAA(raw_data, num.factors, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RNPPAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNPPAAA")
RNMMAAA <- funcaoAAA(raw_data, num.factors, mm,min,max, NULL)
NormStat <- funcaoCompar(RNMMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMMAAA")
RNMEAAA <- funcaoAAA(raw_data, num.factors, me,mean,max,min)
NormStat <- funcaoCompar(RNMEAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RNMEAAA")

# raw log
RLDMAAA <- funcaoAAA(raw_data_log, num.factors, div, median, NULL, NULL)
NormStat <- funcaoCompar(RLDMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDMAAA")
RLDAAAA <- funcaoAAA(raw_data_log, num.factors, div, mean, NULL, NULL)
NormStat <- funcaoCompar(RLDAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDAAAA")
RLDSAAA <- funcaoAAA(raw_data_log, num.factors, div,sum, NULL, NULL)
NormStat <- funcaoCompar(RLDSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLDSAAA")
RLSMAAA <- funcaoAAA(raw_data_log, num.factors, sub,median, NULL, NULL)
NormStat <- funcaoCompar(RLSMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSMAAA")
RLSAAAA <- funcaoAAA(raw_data_log, num.factors, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(RLSAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSAAAA")
RLSSAAA <- funcaoAAA(raw_data_log, num.factors, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(RLSSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLSSAAA")
RLRMAAA <- funcaoAAA(raw_data_log, num.factors, root,median, NULL, NULL)
NormStat <- funcaoCompar(RLRMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRMAAA")
RLRAAAA <- funcaoAAA(raw_data_log, num.factors, root,mean, NULL, NULL)
NormStat <- funcaoCompar(RLRAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRAAAA")
RLRSAAA <- funcaoAAA(raw_data_log, num.factors, root,sum, NULL, NULL)
NormStat <- funcaoCompar(RLRSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLRSAAA")
RLLMAAA <- funcaoAAA(raw_data_log, num.factors, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(RLLMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLMAAA")
RLLAAAA <- funcaoAAA(raw_data_log, num.factors, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(RLLAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLAAAA")
RLLSAAA <- funcaoAAA(raw_data_log, num.factors, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(RLLSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLLSAAA")
RLZSAAA <- funcaoAAA(raw_data_log, num.factors, zs,mean,sd, NULL)
NormStat <- funcaoCompar(RLZSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLZSAAA")
RLPPAAA <- funcaoAAA(raw_data_log, num.factors, pp,mean,sd, NULL)
NormStat <- funcaoCompar(RLPPAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLPPAAA")
RLMMAAA <- funcaoAAA(raw_data_log, num.factors, mm,min,max, NULL)
NormStat <- funcaoCompar(RLMMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMMAAA")
RLMEAAA <- funcaoAAA(raw_data_log, num.factors, me,mean,max,min)
NormStat <- funcaoCompar(RLMEAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "RLMEAAA")

# prog no log
PNDMAAA <- funcaoAAA(Progenesis_data, num.factors, div, median, NULL, NULL)
NormStat <- funcaoCompar(PNDMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDMAAA")
PNDAAAA <- funcaoAAA(Progenesis_data, num.factors, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PNDAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDAAAA")
PNDSAAA <- funcaoAAA(Progenesis_data, num.factors, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PNDSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNDSAAA")
PNSMAAA <- funcaoAAA(Progenesis_data, num.factors, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PNSMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSMAAA")
PNSAAAA <- funcaoAAA(Progenesis_data, num.factors, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PNSAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSAAAA")
PNSSAAA <- funcaoAAA(Progenesis_data, num.factors, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PNSSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNSSAAA")
PNRMAAA <- funcaoAAA(Progenesis_data, num.factors, root,median, NULL, NULL)
NormStat <- funcaoCompar(PNRMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRMAAA")
PNRAAAA <- funcaoAAA(Progenesis_data, num.factors, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PNRAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRAAAA")
PNRSAAA <- funcaoAAA(Progenesis_data, num.factors, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PNRSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNRSAAA")
PNLMAAA <- funcaoAAA(Progenesis_data, num.factors, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PNLMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLMAAA")
PNLAAAA <- funcaoAAA(Progenesis_data, num.factors, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PNLAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLAAAA")
PNLSAAA <- funcaoAAA(Progenesis_data, num.factors, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PNLSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNLSAAA")
PNZSAAA <- funcaoAAA(Progenesis_data, num.factors, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PNZSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNZSAAA")
PNPPAAA <- funcaoAAA(Progenesis_data, num.factors, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PNPPAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNPPAAA")
PNMMAAA <- funcaoAAA(Progenesis_data, num.factors, mm,min,max, NULL)
NormStat <- funcaoCompar(PNMMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMMAAA")
PNMEAAA <- funcaoAAA(Progenesis_data, num.factors, me,mean,max,min)
NormStat <- funcaoCompar(PNMEAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PNMEAAA")

# prog log
PLDMAAA <- funcaoAAA(Progenesis_data_log, num.factors, div, median, NULL, NULL)
NormStat <- funcaoCompar(PLDMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDMAAA")
PLDAAAA <- funcaoAAA(Progenesis_data_log, num.factors, div, mean, NULL, NULL)
NormStat <- funcaoCompar(PLDAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDAAAA")
PLDSAAA <- funcaoAAA(Progenesis_data_log, num.factors, div,sum, NULL, NULL)
NormStat <- funcaoCompar(PLDSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLDSAAA")
PLSMAAA <- funcaoAAA(Progenesis_data_log, num.factors, sub,median, NULL, NULL)
NormStat <- funcaoCompar(PLSMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSMAAA")
PLSAAAA <- funcaoAAA(Progenesis_data_log, num.factors, sub,mean, NULL, NULL)
NormStat <- funcaoCompar(PLSAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSAAAA")
PLSSAAA <- funcaoAAA(Progenesis_data_log, num.factors, sub,sum, NULL, NULL)
NormStat <- funcaoCompar(PLSSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLSSAAA")
PLRMAAA <- funcaoAAA(Progenesis_data_log, num.factors, root,median, NULL, NULL)
NormStat <- funcaoCompar(PLRMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRMAAA")
PLRAAAA <- funcaoAAA(Progenesis_data_log, num.factors, root,mean, NULL, NULL)
NormStat <- funcaoCompar(PLRAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRAAAA")
PLRSAAA <- funcaoAAA(Progenesis_data_log, num.factors, root,sum, NULL, NULL)
NormStat <- funcaoCompar(PLRSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLRSAAA")
PLLMAAA <- funcaoAAA(Progenesis_data_log, num.factors, logarit,median, NULL, NULL)
NormStat <- funcaoCompar(PLLMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLMAAA")
PLLAAAA <- funcaoAAA(Progenesis_data_log, num.factors, logarit,mean, NULL, NULL)
NormStat <- funcaoCompar(PLLAAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLAAAA")
PLLSAAA <- funcaoAAA(Progenesis_data_log, num.factors, logarit,sum, NULL, NULL)
NormStat <- funcaoCompar(PLLSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLLSAAA")
PLZSAAA <- funcaoAAA(Progenesis_data_log, num.factors, zs,mean,sd, NULL)
NormStat <- funcaoCompar(PLZSAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLZSAAA")
PLPPAAA <- funcaoAAA(Progenesis_data_log, num.factors, pp,mean,sd, NULL)
NormStat <- funcaoCompar(PLPPAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLPPAAA")
PLMMAAA <- funcaoAAA(Progenesis_data_log, num.factors, mm,min,max, NULL)
NormStat <- funcaoCompar(PLMMAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMMAAA")
PLMEAAA <- funcaoAAA(Progenesis_data_log, num.factors, me,mean,max,min)
NormStat <- funcaoCompar(PLMEAAA, num.factors, num.levels, num.replic, tst_thrsh, NormStat, empty.NormStat, bal.design, rep.meas, "PLMEAAA")


#alguns tem warnings pois tem algumas proteinas que todas as observacoes sao 0
#isso faz com que algumas contas sejam divisoes por 0, que causa resultados NaN (not a number)
#ou log(0) que da -Inf







######
#salvando tudo em csv

ctr <- c("RN.....", "RL.....", "PN.....", "PL.....") #W tabelas controle
l1 <- rep(rep(c("R","P"), each=32),7) #W O último número (atualmente 2) representa o número de condições de normalização testadas. A cada condição acrescentada em l5, aumentar uma unidade nesse número.
l2 <- rep(rep(rep(c("N","L"), each=16),2),7) #W O último número (atualmente 2) representa o número de condições de normalização testadas. A cada condição acrescentada em l5, aumentar uma unidade nesse número.
l3 <- rep(rep(c(rep(c("D","S","R","L"),each=3),"Z","P","M","M"),4),7) #W O último número (atualmente 2) representa o número de condições de normalização testadas. A cada condição acrescentada em l5, aumentar uma unidade nesse número.
l4 <- rep(rep(c(rep(c("M","A","S"),4),"S","P","M","E"),4),7) #W O último número (atualmente 2) representa o número de condições de normalização testadas. A cada condição acrescentada em l5, aumentar uma unidade nesse número.
l5 <- rep(c("1A1","A11","11A","1AA","A1A","AA1","AAA"), each=64) #W acrescentar 1A1, 11A, 1AA
tabelas <- paste(paste(paste(paste(l1,l2,sep=""),l3,sep=""),l4,sep=""),l5,sep="")
tabelas <- c(ctr, tabelas)

#escrever o caminho pra onde salvar no pc
caminho  <- paste(dlgDir(default = getwd(), title='Choose a folder to save the RESULT files:')$res, "/", sep="") 
for(i in 1:length(tabelas)){
  nome <- paste(caminho,tabelas[i],sep="")
  local <- paste(nome, ".csv",sep="")
  write.csv(get(tabelas[i]),local)
}
write.csv(NormStat, paste(caminho, "NormStat.csv"), sep="")



#W #######################################################################
#W  CHANGELOG
#W #######################################################################
#W  w1 - Criar novo dataframe (chamado NormStat) contendo:
#W       1a col = método de normalização (ex: PLMEAA1)
#W       2a col = teste estatístico
#W       3a col = subdiv. do teste (ex. p/ ANOVA 2-way: fator 1, fator 2, interação)
#W       4a col = limite do teste
#W       5a col = quantas ptn passam no teste
#W       a partir da 6a col - 1 coluna p/ cada PTN
#W       linhas (definidas em funçào da 2a e 3a col): 1-way ANOVA, 2-way ANOVA (fator1, fator2, interação), Tukey HSD (para as várias combinações pareadas de condições), teste t (várias comb), Benjamini-Hochberg (várias comb), Limma ?, Rank Products?
#W       Esse dataframe deve ser criado aos poucos, em uma função, chamada depois de cada função para criar as tabelas de normalização, inserindo o conjunto de linhas (aprox.21?) referente a cada método de normalização - ver planilha MOCK.
#W       tentativa de implementar tudo no mesmo script (3-way (ou fatorial), 2-way, 1-way, t-test), dependendo do número de fatores apresentado como input do usuário
#W       O código começou a ficar muito grande e confuso, pois foi surgindo uma série de "if" nested, tanto no código geral quanto na função para calcular normalização.
#W       Em 19/5 pensei na possibilidade de fazer um script para cada número de fatores, sendo o w2_1way, w2_2way, w2_3way e w2_fatorial. 
#W       Vou primeiro tentar fazer o w2 (completo, sem dividir), mas se começar a dar problemas, mudo para esses tipos diferentes de arquivos.
#W       Ajustada uma das funções de normalização (A11).
#W       Concluído o W1 com ANOVA de 1 e 2 vias + Shapiro-Wilk + Levene funcionando para A11, mas com o código ainda não colocado em função.
#W 
#W  w2 - enviar informações para a função de normalização: número de fatores, número de níveis (p/ cada fator) e número de replicatas. O nome de cada fator vai ser obtido na própria função, pois a tabela contém isso.
#W       FUNCIONANDO com:
#W         A11 e AA1, 64 variações para cada.
#W         testes estat = 1-way ou 2-way ANOVA, Levene, Shapiro-Wilk
#W         escolha automática do teste a partir do número de fatores
#W       PNDMAA1 e PNDMAA1 comparados com Excel - OK.
#W 
#W  w3 - incluido Tukey HSD (post-hoc p/ 1-way e p/ 2-way se um dos fatores tiver mais que 2 níveis)
#W       incluidos scores (score 2-way = soma núm. de significativas p/ cada fator + signif p/ interação + signif p/ Levene + Signif p/ Shapiro-Wilk) (score 1-way = soma núm. de significativas ANOVA + signif p/ Levene + Signif p/ Shapiro-Wilk) (score teste t, ou B-H = soma núm. de significativas + signif p/ Levene + Signif p/ Shapiro-Wilk)
#W       FUNCIONANDO
#W 
#W  w4 - incluido teste t (se tiver 1 fator e 2 níveis) 
#W       com Shapiro-Wilk e com Levene
#W       pareado e não pareado
#W       FUNCIONANDO
#W
#W  w5 - incluido AAA
#W       
#W  w6 - incluido A1A
#W       
#W  w7 - incluido 11A
#W       
#W  w8 - incluido 1A1
#W       
#W  w9 - incluido 1AA
#W       
#W  w10 - incluida imputação
#W       
#W  w11 - incluir ART p/ 2-way ANOVA sem distribuição normal
#W       
#W       
#W
#W #######################################################################
#W  TO DO
#W #######################################################################
#W
#W  DONE - 1- Acrescentar as outras estratégias de agrupamento para normalização 1AA, 1A1, A11, 11A, A1A, AA1 e AAA já estão prontas, 111 não faz sentido (corresponde aos controles sem normalização)
#W  DONE - 2- Converter a planilha original exportada do Progenesis para os dados que este script usa
#W  3- Implementar MANOVA? São diversas variáveis dependentes (cada proteína é uma) e mShapiro para fazer o Shapiro de multivariada. MANOVA vai dezer se há alguma associação entre as condições testadas e o relacionamento interno do conjunto de proteínas (ver o segundo link - statistics by Jim). Depois é necessário fazer ANOVA para cada uma (o que este script já faz).
#W             https://statisticsbyjim.com/anova/multivariate-anova-manova-benefits-use/
#W             http://www.sthda.com/english/wiki/manova-test-in-r-multivariate-analysis-of-variance
#W  4- Achar jeitos eficazes de mostrar os resultados de estatística, para facilitar a definição do melhor método de normalização em questão (gráficos? score?)
#W  5- Ver se convém, e como implementar outras formas e/ou combinações de normalização, scaling e standardization - independentes ou combinados - estudar diferenças (standardization: média=0, stdev=1) (scaling: todos na mesma escala - ex. entre 0 e 1) (normalization:fit para uma distribuição, frequentemente normal, mas pode ser outra - Poisson ou outra), BoxCox?
#W      https://stats.stackexchange.com/questions/35591/normalization-vs-scaling
#W      https://medium.com/@nsethi610/data-cleaning-scale-and-normalize-data-4a7c781dd628
#W  6- Implementar Welch ANOVA para variâncias não homogêneas e Games-Howell multiple comparisons method para fazer o post-hoc desse caso. https://statisticsbyjim.com/anova/welchs-anova-compared-to-classic-one-way-anova/#more-708
#W  DONE - 7- Implementar imputation p/ substituir zeros - conta o número de zeros na tabela, determina o mínimo da tabela toda (exceto zeros), define um vetor de aleatórios (runif) entre o mínimo e mínimo/1000, com o número de componentes igual ao número de zeros, substitui cada um dos zeros por um dos componentes do vetor
#W  8- Ajustar o score para os métodos. Por enquanto soma as proteínas significativas para cada fator, interação , Levene e Shapiro. Outro score poderá considerar também o número de significativas para cada par de comparações no post-hoc
#W  9- Repeatd measures ANOVA p/ 1-way, 2-way e 3-way (condição similar ao teste t pareado = medidas no mesmo indivíduo), Mauchly's test - elipticidade (variáncia conservada entre os pares) (https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/)
#W 10- Bayesian ANOVA (https://psyarxiv.com/spreb/download) (https://www.flutterbys.com.au/stats/tut/tut7.4b.html)
#W 11- Inverter a ordem: calcular primeiro teste de normalidade, depois de homogeneidade de variâncias, depois os testes ANOVA ou T. O tipo destes deve ser definido em função dos resultados dos outros.
#W 12- Incluir teste de outliers. Para ANOVA de medidas repetidas precisa - identify_outliers() [rstatix package]. (https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/) ou https://www.r-bloggers.com/identify-describe-plot-and-remove-the-outliers-from-the-dataset/#:~:text=To%20detect%20the%20outliers%20I,data%20with%20and%20without%20outliers.
#W 13- transformar teste t em função
#W 14- usar teste t como post-hoc para 1-way, 2-way e 3-way ANOVA se algum dos fatores tiver mais que 3 níveis, mas corrigir o p-valor p/ repetições Benjamini-Hockberg e/ou FDR (p.adjust)
#W 15- extrair as informações de condições, replicatas, desenho balanceado direto da planilha do Progenesis. Acho que a única informação que precise de input é o número de fatores (talvez os níveis de cada fator, se tiver mais que 1)
#W 16- Ver se convém implementar nested ANOVA (http://www.biostathandbook.com/nestedanova.html)
#W 17- opção de criar uma planilha de resultados com todos os p-valores corrigidos (corrige todos os p-valores de NormStat), com opção de método de correção e de número de repetições (número de pares testados -> correção válida só para os post-hoc) (número de proteínas; número de proteínas significativas) - comando = p.adjust
#W 18- opção p/ 2-way não paramétrico (não normal) = ARTool (http://depts.washington.edu/acelab/proj/art/index.html), c/ post-hoc emmeans, c/ opção p/ medidas repetidas ou não - ou ver outros exemplos de código abaixo dos use cases (https://www.researchgate.net/post/Is_there_a_non-parametric_equivalent_of_a_2-way_ANOVA)
#W 19- bootstrapping?
#W 20- opções para 2-way heteroscedático (variâncias não homogêneas): permutation ANOVA (https://rdrr.io/cran/RVAideMemoire/man/perm.anova.html); Robust ANOVA (https://rcompanion.org/rcompanion/d_08a.html)
#W 21- organizar testes (ver esquema abaixo)
#W 22- Incluir cálculo da potência do teste (https://www.statmethods.net/stats/power.html) (https://advstats.psychstat.org/book/power/index.php)
#W 23- Three-way ANOVA (inclusive opções para não normal, variâncias não homogêneas, não balanceado e repeated measures)
#W 24- Full factorial ANOVA???; n-way ANOVA???
#W 25- Implementar post-hoc para comparações inter-fatores (atualmente só está implementado intra-fatores) ver se a função do Tukey já fz isso automático, sem especificar o which. Para o não paramétrico (ART), usar o testInteractions (library phia). "is the difference between a and b in condition d significantly different from the difference between a and b in condition e?"
#W 26- Implementar: Kruskal-Wallis, Welch, Friedman, Benjamini-Hochberg, automação do uso de paramétricos ou não p/ cada ptn, Limma?, Rank Products?
#W
#W
#W
#W
#W #######################################################################
#W  USE CASES
#W #######################################################################
#W
#W  1- Análises com 2 fatores, 2 níveis em cada fator (ex: projeto Isabelle M x F)
#W  2- Análises com 2 fatores, mais que 2 níveis em um dos fatores (ex: projeto Wendy manose - manose/glicose; ctrl/PMA/fMLP)
#W  3- Análises com 1 fator, 2 níveis (ex: projetos Phillippe IL-8; Adriano TNF)
#W  4- Análises com 2 fatores, mais que 2 níveis em cada fator (ex: projeto Hylane ctrl/PMA s_inib/inib_primeiro/inib_depois)


#W #######################################################################
#W  REFERENCES (to cite in our articles that use this script)
#W #######################################################################
#W
#W - ART - Wobbrock, J.O., Findlater, L., Gergle, D. and Higgins, J.J. (2011). The Aligned Rank Transform for nonparametric factorial analyses using only ANOVA procedures. Proceedings of the ACM Conference on Human Factors in Computing Systems (CHI '11). Vancouver, British Columbia (May 7-12, 2011). New York: ACM Press, pp. 143-146.
#W - Site com algumas referências úteis para citar http://depts.washington.edu/acelab/proj/art/index.html






# if(!require(rcompanion)){install.packages("rcompanion")}
# if(!require(WRS2)){install.packages("WRS2")}
# if(!require(ARTool)){install.packages("ARTool")}
# if(!require(car)){install.packages("car")}
# ### Data from Sokal and Rohlf
# Value = c(709,679,699,657,594,677,592,538,476,508,505,539)
# Sex = c(rep("Male",3), rep("Female",3), rep("Male",3), rep("Female",3))
# Fat = c(rep("Fresh", 6), rep("Rancid", 6))
# Sokal = data.frame(Value, Sex, Fat)
# ### Scheirer-Ray-Hare test
# library(rcompanion)
# scheirerRayHare(Value ~ Sex * Fat, data=Sokal)
# ### Aligned ranks transformation anova
# library(ARTool)
# model = art(Value ~ Sex * Fat, data=Sokal)
# anova(model)
# ### 2-way median test
# library(WRS2)
# med2way(Value ~ Sex * Fat, data=Sokal)
# ### General linear model
# modelm = lm(Value ~ Sex * Fat, data=Sokal)
# library(car)
# Anova(modelm)
# hist(residuals(modelm), col="darkgray")
# ### Normal scores transformation
# Elfving = qnorm((rank(Value)-pi/8)/(length(Value)-pi/4+1))
# modelint = lm(Elfving ~ Sex * Fat, data=Sokal)
# library(car)
# Anova(modelint)
# hist(residuals(modelint), col="darkgray")
# plot(predict(modelint), residuals(modelint))
# 



# 
#                S                                             N
# 2ou3 fatores  -----distribuição normal (Shapiro-Wilk > 0.5) ---> não paramétrico usar ART package; post-hoc = emmeans (no ART) (http://depts.washington.edu/acelab/proj/art/index.html)
#               | S                                                                  N
#               |----homogeneidade e variâncias (homoscedasticidade) (Levene > 0.05)---> permutation ANOVA (https://rdrr.io/cran/RVAideMemoire/man/perm.anova.html) or Robust ANOVA (https://rcompanion.org/rcompanion/d_08a.html)
#               || S                                                                      N
#               ||---estudo balanceado (mesmo número de replicatas em todas as condições)---> car package, type III sum of squares (http://www.sthda.com/english/wiki/two-way-anova-test-in-r#tukey-multiple-pairwise-comparisons)
#               ||| N                                                                                               S
#               |||--medidas repetidas (within subject; as diversas condições foram medidas para o mesmo indivíduo)---> type III package rstatix (https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/)
#               ||||
#               2-way (ou 3-way) ANOVA, post-hoc = Tukey
#                              
#               S                                             N
# 1 fator       -----distribuição normal (Shapiro-Wilk > 0.5) ---> não paramétrico usar Kruskal-Wallis; post-hoc = Wilcoxon pairwise (http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r)
# >2 níveis     | S                                                                  N
#               |----homogeneidade e variâncias (homoscedasticidade) (Levene > 0.05)---> Welch ANOVA ; post-hoc = Games-Howell package = userfriendlyscience (https://statisticsbyjim.com/anova/welchs-anova-compared-to-classic-one-way-anova/#more-708) (https://rpubs.com/aaronsc32/games-howell-test) (https://cran.r-project.org/web/packages/userfriendlyscience/index.html) (http://www.sthda.com/english/wiki/one-way-anova-test-in-r)
#               || S ou N                                                                      
#               ||---estudo balanceado (mesmo número de replicatas em todas as condições)---> o fato de ser balanceado ou não, não muda o teste. estudo balanceado tem maior potência.
#               ||| N                                                                                               S
#               |||--medidas repetidas (within subject; as diversas condições foram medidas para o mesmo indivíduo)---> type III package rstatix (https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/) ou package car (https://www.r-bloggers.com/r-tutorial-series-one-way-repeated-measures-anova/) ou package userfriendlyscience (https://cran.r-project.org/web/packages/userfriendlyscience/userfriendlyscience.pdf)
#               ||||
#               1-way ANOVA, post-hoc = Tukey
#                                                                                                                                   
#                                                                                                                                   
#                    
#               S                                             N
# 1 fator       -----distribuição normal (Shapiro-Wilk > 0.5) ---> não paramétrico usar Wilcoxon (=Mann-Whitney) (http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r)
# 2 níveis      | S                                                                  N
#               |----homogeneidade e variâncias (homoscedasticidade) (Levene > 0.05)---> Welch test usa o mesmo t.test, mas com a opção var.equal=FALSE (http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r)
#               || S ou N                                                                      
#               ||---estudo balanceado (mesmo número de replicatas em todas as condições)---> o fato de ser balanceado ou não, não muda o teste. estudo balanceado tem maior potência.
#               ||| N                                                                                                                         S
#               |||--dependente (medidas repetidas ou pareado, ou within subject; as diversas condições foram medidas para o mesmo indivíduo)---> usa o mesmo t.test, mas com a opção paired=TRUE (http://www.sthda.com/english/wiki/t-test)
#               ||||
#               t-test
#                                                                                                                                                             
#                    
#                    