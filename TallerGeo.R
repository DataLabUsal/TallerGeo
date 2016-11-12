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
#Veamos qué sucede aplicando una transformación logarítmica
plot(density(log(meuse$zinc)),main='Densidad Zn',col="blue")
plot(density(log(meuse$cadmium)),main='Densidad Cd',col="orange")
plot(density(log(meuse$copper)),main='Densidad Cu',col="red")
plot(density(log(meuse$lead)),main='Densidad Pb',col="grey")

#Nuestro principal objetivo es ver el comportamiento de las variables Zinc, Cadmio, Cobre y Plomo.
AnData <- meuse[,c(1:6,8)] 
cor(AnData[,3:6])	#Matriz de correlacion

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
vg.zn <- variogram(log(zinc)~1, AnData)
dev.off() #Borramos los gráficos anteriores
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
zn.krig <- krige(log(zinc)~1, AnData, meuse.grid, model = vg.zn.fit)
names(zn.krig)
p1 <- spplot(zn.krig["var1.pred"])
p2 <- spplot(zn.krig["var1.var"])
print(p1, split=c(1,1,2,1), more=TRUE)
print(p2, split=c(2,1,2,1), more=TRUE)

#Crossvalidation#
cv.ok <- krige.cv(log(zinc)~1, AnData, model=vg.zn.fit)
summary(cv.ok)
resid.ok <- as.data.frame(cv.ok)$residual
sqrt(mean(resid.ok^2))
mean(resid.ok)
mean(resid.ok^2/as.data.frame(cv.ok)$var1.var)

##Universal Kriging##
vg.zn.dist <- variogram(log(zinc)~sqrt(dist), AnData)
vg.zn.dist.fit <- fit.variogram(vg.zn.dist, model = vgm(1, "Sph", 900, 1))
zn.krig.dist <- krige(log(zinc)~sqrt(dist), AnData, meuse.grid, model = vg.zn.dist.fit)

#Crossvalidation#
cv.uk <- krige.cv(log(zinc)~sqrt(dist), AnData, model=vg.zn.dist.fit)
summary(cv.uk)
resid.uk <- as.data.frame(cv.uk)$residual
sqrt(mean(resid.uk^2))
mean(resid.uk)
mean(resid.uk^2/as.data.frame(cv.uk)$var1.var)

##COKRIGING##
#Ordinary Cokriging#
gock <- gstat(NULL, "zn", log(zinc)~1, AnData)
gock <- gstat(gock, "cd", log(cadmium)~1, AnData)
gock <- gstat(gock, "cu", log(copper)~1, AnData)
gock <- gstat(gock, "pb", log(lead)~1, AnData)

#Se puede realizar un Cokriging Universal de manera parecida a cómo se realizó en el caso del kriging#
#Universal Cokriging#
guck <- gstat(NULL, "zn", log(zinc)~sqrt(dist), AnData)
guck <- gstat(guck, "cd", log(cadmium)~sqrt(dist), AnData)
guck <- gstat(guck, "cu", log(copper)~sqrt(dist), AnData)
guck <- gstat(guck, "pb", log(lead)~sqrt(dist), AnData)

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

#Visualizamos los resultados
pl1 <- spplot(cku["zn.pred"], main="Zinc Universal")
pl1b <- spplot(cko["zn.pred"], main="Zinc Ordinary")
pl2 <- spplot(cku["cu.pred"], main="Cobre Universal")
pl2b <- spplot(cko["cu.pred"], main="Cobre Ordinary")
pl3 <- spplot(cku["cd.pred"], main="Cadmio Universal")
pl3b <- spplot(cko["cd.pred"], main="Cadmio Ordinary")
pl4 <- spplot(cku["pb.pred"], main="Plomo Universal")
pl4b <- spplot(cku["pb.pred"], main="Plomo Ordinary")

print(pl1, split = c(1,1,4,2), more=TRUE)
print(pl1b, split = c(2,1,4,2), more=TRUE)
print(pl2, split = c(3,1,4,2), more=TRUE)
print(pl2b, split = c(4,1,4,2), more=TRUE)
print(pl3, split = c(1,2,4,2), more=TRUE)
print(pl3b, split = c(2,2,4,2), more=TRUE)
print(pl4, split = c(3,2,4,2), more=TRUE)
print(pl4b, split = c(4,2,4,2))

#Crossvalidation#
#Ordinary#
cv.ck.ok <- gstat.cv(gock.fit)
summary(cv.ck.ok)
resid.cko <- as.data.frame(cv.ck.ok)$residual
sqrt(mean(resid.cko^2))
mean(resid.cko)
mean(resid.cko^2/as.data.frame(cv.ck.ok)$zn.var)
##Si quieres calcular todos los residuales añade la opción "all.residuals = TRUE" en gstat.cv#
cv.ck.ok.all <- gstat.cv(gock.fit,all.residuals=TRUE)
summary(cv.ck.ok.all)


#Universal#
cv.ck.uk <- gstat.cv(guck.fit)
summary(cv.ck.uk)
resid.cku <- as.data.frame(cv.ck.uk)$residual
sqrt(mean(resid.cku^2))
mean(resid.cku)
mean(resid.cku^2/as.data.frame(cv.ck.uk)$zn.var)
#Todos los residuales#
cv.ck.uk.all <- gstat.cv(gock.fit,all.residuals=TRUE)
summary(cv.ck.uk.all)
