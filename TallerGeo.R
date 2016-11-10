####TALLER GEOESTADISTICA####

####PAQUETES PARA LA INSTALACIÓN####

install.packages('sp')		#Paquete específico para la representación espacial
install.packages('gstat')	#Paquete más usado para análisis geoestadístico en R

library(sp)
library(gstat)

#Explorando los datos
data(meuse)		#Datos río Meuse

head(meuse)		#Encabezado de los datos
summary(meuse)	#Resumen de los datos

par(mfrow=c(2,2))
plot(density(meuse$zinc),main='Densidad Zn',col="blue")
plot(density(meuse$cadmium),main='Densidad Cd',col="orange")
plot(density(meuse$copper),main='Densidad Cu',col="red")
plot(density(meuse$lead),main='Densidad Pb',col="grey")

#Nuestro principal objetivo es ver el comportamiento de las variables Zinc, Cadmio, Cobre y Plomo.
AnData <- meuse[,c(1:6,8)] 
cor(AnData[,3:6])	#Matriz de correlacion

#Estandarizamos los datos
medias	<-	apply(AnData[,3:6],2,mean)
desv	<-	apply(AnData[,3:6],2,sd)

for (i in 1:4)
{AnData[,i+2] <- (AnData[,i+2]-medias[i])/desv[i]}

boxplot(AnData$cadmium,AnData$copper,AnData$lead,AnData$zinc,names=c('Cd','Cu','Pb','Zn'))

#Transformación espacial de los datos
coordinates(AnData) <- ~x+y
coordinates(meuse) <- ~x+y

#Gráficos Bubble
b1 <- bubble(meuse, "cadmium",col="orange",alpha=0.4, main = "Ca concentrations (ppm)")
b2 <- bubble(meuse, "copper",col="red",alpha=0.4, main = "Cu concentrations (ppm)")
b3 <- bubble(meuse, "lead",col="grey",alpha=0.4, main = "Pb concentrations (ppm)")
b4 <- bubble(meuse, "zinc",col="blue",alpha=0.4, main = "Zn concentrations (ppm)")

print(b1, split = c(1, 1, 2, 2), more = TRUE)
print(b2, split = c(2, 1, 2, 2), more = TRUE)
print(b3, split = c(1, 2, 2, 2), more = TRUE)
print(b4, split = c(2, 2, 2, 2), more = FALSE)

##plotGoogleMaps##
#Visualización sobre mapa 
install.packages('plotGoogleMaps')
library(plotGoogleMaps)

#Añadimos las referencias geográficas a nuestros datos
proj4string(meuse) <- CRS('+init=epsg:28992') 

#Creamos nuestros mapas web
m1	<-	plotGoogleMaps(meuse,filename='MiMapa.html')
m2	<-	bubbleGoogleMaps(meuse,zcol='zinc', max.radius = 80,filename='MiMapa2.html')
m3	<-	plotGoogleMaps(meuse,filename='MiMapa3.htm',iconMarker='http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png',mapTypeId='ROADMAP',layerName = 'MEUSE POINTS') 


##KRIGING##
#Variograma#
vg.zn <- variogram(zinc~1, AnData)
plot(vg.zn)
vg.zn

#Ajustamos una función a nuestro variograma
vg.zn.fit <- fit.variogram(vg.zn, model = vgm(1, "Sph", 900, 1))
vg.zn.fit

#Visualizamos el ajuste
plot(vg.zn,vg.zn.fit)

#Mallado específico para nuestros datos
data(meuse.grid)
head(meuse.grid)
coordinates(meuse.grid) <- ~x+y
plot(meuse.grid)

##Ordinary Kriging##
zn.krig <- krige(zinc~1, AnData, meuse.grid, model = vg.zn.fit)
names(zn.krig)
p1 <- spplot(zn.krig["var1.pred"])
p2 <- spplot(zn.krig["var1.var"])
print(p1, split=c(1,1,2,1), more=TRUE)
print(p2, split=c(2,1,2,1), more=TRUE)

#Crossvalidation#
cv.ok <- krige.cv(zinc~1, AnData, model=vg.zn.fit)
summary(cv.ok)
resid.ok <- as.data.frame(cv.ok)$residual
sqrt(mean(resid.ok^2))
mean(resid.ok)
mean(resid.ok^2/as.data.frame(cv.ok)$var1.var)

##Universal Kriging##
vg.zn.dist <- variogram(zinc~sqrt(dist), AnData)
vg.zn.dist.fit <- fit.variogram(vg.zn.dist, model = vgm(1, "Sph", 900, 1))
zn.krig.dist <- krige(zinc~sqrt(dist), AnData, meuse.grid, model = vg.zn.dist.fit)

#Crossvalidation#
cv.uk <- krige.cv(zinc~sqrt(dist), AnData, model=vg.zn.dist.fit)
summary(cv.uk)
resid.uk <- as.data.frame(cv.uk)$residual
sqrt(mean(resid.uk^2))
mean(resid.uk)
mean(resid.uk^2/as.data.frame(cv.uk)$var1.var)

##COKRIGING##
#Ordinary Cokriging#
gock <- gstat(NULL, "zn", zinc~1, AnData)
gock <- gstat(gock, "cd", cadmium~1, AnData)
gock <- gstat(gock, "cu", copper~1, AnData)
gock <- gstat(gock, "pb", lead~1, AnData)

#Se puede realizar un Cokriging Universal de manera parecida a cómo se realizó en el caso del kriging#
#Universal Cokriging#
guck <- gstat(NULL, "zn", zinc~sqrt(dist), AnData)
guck <- gstat(guck, "cd", cadmium~sqrt(dist), AnData)
guck <- gstat(guck, "cu", copper~sqrt(dist), AnData)
guck <- gstat(guck, "pb", lead~sqrt(dist), AnData)

#Variogramas Cruzados#
#Ordinary Cokriging#
cvock <- variogram(gock)
#Universal Cokriging#
cvuck <- variogram(guck)

#Ajuste del modelo inicial#
#Ordinary Cokriging#
gock <- gstat(gock, model = vgm(1, "Sph", 900, 1), fill.all = TRUE)
#Universal Cokriging#
guck <- gstat(guck, model = vgm(1, "Sph", 900, 1), fill.all = TRUE)

#Modelo Lineal de Corregionalización (LMC)
#Ordinary Cokriging#
gock.fit <- fit.lmc(cvock, gock)
plot(cvock, gock.fit)
#Universal Cokriging#
guck.fit <- fit.lmc(cvuck, guck)
plot(cvuck, guck.fit)

#Cokriging#
#Ordinary Cokriging#
cko <- predict(gock.fit, meuse.grid)
#Universal Cokriging#
cku <- predict(guck.fit, meuse.grid)

#Cross-validation
#Ordinary Cokriging#
crva.cko <- gstat.cv(gock.fit,all.residuals = T)
summary(crva.cko)
plot(density(crva.cko$zn))

#Universal Cokriging#
crva.cku <- gstat.cv(guck.fit,all.residuals = T)
summary(crva.cku)
plot(density(cv.cku$cu))

#Visualizamos los resultados
pl1 <- spplot(cku["zn.pred"], main="Zinc")
pl2 <- spplot(cku["cu.pred"], main="Cobre")
pl3 <- spplot(cku["cd.pred"], main="Cadmio")
pl4 <- spplot(cku["pb.pred"], main="Plomo")

print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))