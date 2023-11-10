# Cargar paquetes
library(vegan)
library(RColorBrewer)
# Cargar datos
data("BCI")
data("BCI.env")
# Funcion para barras de error
error.bar <- function(x,y,upper,lower=upper,length=0.1,type=c("X","Y"),...){if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
  stop("vectors must be same length")
  if(missing(type)){type="Y"}
   if(type=="Y"){
    arrows(x,upper, x,lower, angle=90, code=3, length=length, ...)
  }
  else if(type=="X"){
    arrows(upper,y,lower,y, angle=90, code=3, length=length, ...)
  }
}

# Estimar la diversidad utilizando los números de Hill
# Gamma
# Sumar todos los individuos por especie
BCI.tot<-apply(BCI,2,sum)
q0<-specnumber(BCI.tot)
q1<-exp(diversity(BCI.tot))
q2<-diversity(BCI.tot,"invsimpson")

# Gamma por habitat
# Sumar los individuos por habitat
BCI.hab<-apply(BCI,2,function(x)tapply(x,BCI.env$Habitat,sum))
q0.hab<-specnumber(BCI.hab)
q1.hab<-exp(diversity(BCI.hab))
q2.hab<-diversity(BCI.hab,"invsimpson")

# Diversidad alpha por parcela
q0.parcela<-specnumber(BCI)
q1.parcela<-exp(diversity(BCI))
q2.parcela<-diversity(BCI,"invsimpson")

# Promedio de la diversidad por parcela
q0.prom.hab<-tapply(q0.parcela,BCI.env$Habitat,mean)
q1.prom.hab<-tapply(q1.parcela,BCI.env$Habitat,mean)
q2.prom.hab<-tapply(q2.parcela,BCI.env$Habitat,mean)

# Intervalo de confianza de la diversidad por parcela
q0.UP.hab<-q0.prom.hab+tapply(q0.parcela,BCI.env$Habitat,sd)*1.96
q0.LO.hab<-q0.prom.hab-tapply(q0.parcela,BCI.env$Habitat,sd)*1.96
q1.UP.hab<-q1.prom.hab+tapply(q1.parcela,BCI.env$Habitat,sd)*1.96
q1.LO.hab<-q1.prom.hab-tapply(q1.parcela,BCI.env$Habitat,sd)*1.96
q2.UP.hab<-q2.prom.hab+tapply(q2.parcela,BCI.env$Habitat,sd)*1.96
q2.LO.hab<-q2.prom.hab-tapply(q2.parcela,BCI.env$Habitat,sd)*1.96

# Ver las diferencias por habitat
habitats<-levels(BCI.env$Habitat)
cols<-brewer.pal(5,"Accent")
# q0
plot(1:5,q0.prom.hab,pch=19,xlab="",ylab="Especies",ylim=c(65,110),xaxt="n")
axis(1,1:5,labels = levels(BCI.env$Habitat))
error.bar(x=1:5,y=q0.prom.hab,upper=q0.UP.hab,lower=q0.LO.hab)

#q1
plot(1:5,q1.prom.hab,pch=19,xlab="",ylab="exp(Shannon)",ylim=c(0,80),xaxt="n")
axis(1,1:5,labels = levels(BCI.env$Habitat))
error.bar(x=1:5,y=q1.prom.hab,upper=q1.UP.hab,lower=q1.LO.hab)

#q2
plot(1:5,q2.prom.hab,pch=19,xlab="",ylab="1/D",ylim=c(0,50),xaxt="n")
axis(1,1:5,labels = levels(BCI.env$Habitat))
error.bar(x=1:5,y=q2.prom.hab,upper=q2.UP.hab,lower=q2.LO.hab)

# Diversidad Beta
beta<-vegdist(BCI)
ordenamiento<-cmdscale(beta)

# Veamos los resultados
par(mar=c(4.5,4.5,2,2))
plot(ordenamiento,type="n",xlab="Eje 1",ylab="Eje 2")
for(i in 1:5){
  points(ordenamiento[BCI.env$Habitat==habitats[i],],pch=19,col=cols[i])
}
ordiellipse(ordenamiento,BCI.env$Habitat,col=cols)
par(mar=c(0,0,0,0),fig=c(0,1,0,1),new=TRUE)
plot(0,0,type="n",axes=FALSE)
legend("top",habitats,pch=19,col=cols,bty="n",horiz = TRUE)

# Curvas de acumlación
mods<-by(BCI,BCI.env$Habitat,specaccum,method="rarefaction",simplify=FALSE)

par(mar=c(4.5,4.5,2,2))
plot(0,0,xlim=c(300,11100),ylim=c(85,220),type="n",xlab="Individuos"
     ,ylab="Especies")
for(i in 1:5){
  plot(mods[[i]],xvar="individuals",ci=1.96,add=TRUE,col=cols[i])
}


#Estimadores de riqueza
S.est<-by(BCI,BCI.env$Habitat,specpool,simplify = FALSE)
S.est<-matrix(unlist(S.est),ncol=9,nrow=5,byrow=TRUE,dimnames=list(habitats,names(S.est[[1]])))

plot(1:5,S.est[,"chao"],pch=19,xlab="",ylab="Especies",ylim=c(100,250),xaxt="n")
axis(1,1:5,labels = levels(BCI.env$Habitat))
error.bar(x=1:5,y=S.est[,"chao"]
          ,upper=S.est[,"chao"]+S.est[,"chao.se"]*1.96
            ,lower=S.est[,"chao"]-S.est[,"chao.se"]*1.96)
