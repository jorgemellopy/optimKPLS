#' @title Optimizacion de KPLS mediante un algoritmmo memetico
#' @description Ajuste metaheurístico de parámetros de la regresión KPLS
#' @param datos1 conjunto de datos de variables independientes.
#' @param datos2 conjunto de datos de variables dependientes.
#' @param numC numero de componentes
#' @param Li cota inferior
#' @param Ls cota superior
#' @param Pob numero de poblacion
#' @param maxE numero de evaluacion maxima
#' @return sigma y numero de componente
#' @examples
#' numCompoments <- 7
#' path<-"parkinson.xlsx"
#' dataset <- read_excel(path)
#' database <- data.frame(dataset)
#' data1<-database[2:20]
#' data2<-database[,21:22]
#' lo<-0
#' up<-100
#' po<-40
#' mE<-100
#' Resultado<-optim_kpls(data1,data2,numCompoments,lo,up,po,mE)
#' print(Resultado)
#' @export
optimkpls <- function(datos1,datos2,numC,Li,Ls,Pob,maxE){

  resultado <- matrix(0,1,4)

  tinicial <- proc.time() # Inicia el cronómetro
  ## Calculo de la solucion optima por medio utilizando el algoritmo memetico
  res.fun <- malschains(function(x) {-funcion_objetivo(x,datos1,datos2,numC)}, lower=c(Li), upper=c(Ls), maxEvals=maxE,
                        control=malschains.control(popsize=Pob, istep=100, ls="cmaes", optimum=-5))
  result <- res.fun$sol

  tfinal <- proc.time()-tinicial    # Detiene el cronómetro
  ##
  optimum.value <- funcion_objetivo_2(result,datos1,datos2,numC)

  ##########################################################################
  #Guardar resultados
  ##########################################################################
  resultado[1,1]<-optimum.value[1] #Q2cum maximo funcion 2
  resultado[1,2]<-optimum.value[2] #Num componente
  resultado[1,3]<-optimum.value[3] #degree
  resultado[1,4]<-tfinal[3]        #Tiempo de computo
  Res <- data.frame(resultado)
  return(Res)
}

#' Funcion Objetivo
#' @param x valor de sigma.
#' @param datos1 conjunto de datos de variables independientes.
#' @param datos2 conjunto de datos de variables dependientes.
#' @param numC numero de componentes.
#' @export
funcion_objetivo <- function(x,datos1,datos2,numC){
  vDegree <- x[1]
  if(vDegree==0){
    vDegree <- 0.00001
  }
  kergauss<-rbfdot(sigma = vDegree)
  datos1<-data.matrix(datos1)
  K<-kernelMatrix(kergauss,datos1)
  jd_kpls = plsreg2(K[,],datos2,comps = numC)
  Q2 <- jd_kpls$Q2cum
  Res <- data.frame(Q2)
  v <- c() #Declaramos un vector vacio
  q <- length(Res) #Cantidad de columnas

  for (i in 1:numC){
    for (j in 1:q){
      if(is.null(Q2[i,j])){
        Q2[i,j] <- 0
      }else {
        Q2[i,j] <- Q2[i,j]
      }
    }
  }

  for (i in 1:numC){
    v[i] <- Q2[i,q]
  }

  s <- c() #Declaramos un vector vacio
  s[1] <- max(v)
  for (i in 1:numC){
    if(s[1]==v[i]){
      s[2] <- i
    }
  }

  return(s[1])  # Q2cum -> retorna el valor maximo

}

#' Funcion Objetivo de Evaluacion Final
#' @param x valor de sigma.
#' @param datos1 conjunto de datos de variables independientes.
#' @param datos2 conjunto de datos de variables dependientes.
#' @param numC numero de componentes.
#' @export
funcion_objetivo_2 <- function(x,datos1,datos2,numC){ #Para evaluacion final
  vDegree <- x[1]
  if(vDegree==0){
    vDegree <- 0.00001
  }
  kergauss<-rbfdot(sigma = vDegree)
  datos1<-data.matrix(datos1)
  K<-kernelMatrix(kergauss,datos1)
  jd_kpls = plsreg2(K[,],datos2,comps = numC)
  Q2 <- jd_kpls$Q2cum
  Res <- data.frame(Q2)
  v <- c() #Declaramos un vector vacio
  q <- length(Res) #Cantidad de columnas

  for (i in 1:numC){
    for (j in 1:q){
      if(is.null(Q2[i,j])){
        Q2[i,j] <- 0
      }else {
        Q2[i,j] <- Q2[i,j]
      }
    }
  }

  for (i in 1:numC){
    v[i] <- Q2[i,q]
  }

  s <- c() #Declaramos un vector vacio
  s[1] <- max(v)
  for (i in 1:numC){
    if(s[1]==v[i]){
      s[2] <- i
      s[3] <- vDegree
    }
  }

  return(s)  # retorna vector de valores

}
