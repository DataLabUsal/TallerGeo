###Taller Geoestadística (Ejemplo 2)###

##Instalar paquetes en caso de no haberlo hecho anteriormente#
install.packages('sp')		#Paquete específico para la representación espacial
install.packages('gstat')	#Paquete más usado para análisis geoestadístico en R
install.packages('RCurl')	#Paquete HTTP requests

#Librerías#
library(sp)
library(gstat)
library(Rcurl)

#Descargar archivo SunshineMine.csv que puede recuperarse desde la url > https://github.com/DataLabUsal/TallerGeo#
url1 <- getURL('https://raw.githubusercontent.com/DataLabUsal/TallerGeo/master/SunshineMine.csv')
AnData <- read.csv(text = url1)
#Recuerda poner el camino de tu archivo. Si no, a través de RStudio puedes descargarlo desde la pestaña "Import Dataset"#


#Exploración de los datos#
head(AnData)	#Encabezado de los datos
summary(AnData)	#Resumen de los datos

par(mfrow=c(1,2))
plot(density(AnData$Au),main='Densidad Au',col="red")
plot(density(AnData$Ag),main='Densidad Ag',col="blue")
#Veamos lo que sucede aplicando una transformación logarítmica a nuestros datos
plot(density(log(AnData$Au)),main='Densidad Au',col="red")
plot(density(log(AnData$Ag)),main='Densidad Ag',col="blue")

#Transformamos a SpatialPointsDataFrame (librería 'sp')
coordinates(AnData) <- ~Easting+Elevation

#Gráficos de concentración
b1 <- bubble(AnData, "Au",col="red",alpha=0.4, main = "Au concentrations")
b2 <- bubble(AnData, "Ag",col="blue",alpha=0.4, main = "Ag concentrations")
print(b1, split = c(1, 1, 1, 2), more = TRUE)
print(b2, split = c(1, 2, 1, 2), more = TRUE)

##KRIGING##
#Comenzamos con un Ordinary Kriging#
#Variograma#
vg.au <- variogram(log(Au)~1, AnData)
plot(vg.au)
vg.au

#Ajustamos una función a nuestro variograma
vg.au.fit <- fit.variogram(vg.au, model = vgm(1, "Sph", 400, 1))
vg.au.fit

#Visualizamos el ajuste
plot(vg.au,vg.au.fit)

#Mallado específico para nuestros datos (se puede descargar desde la página https://github.com/DataLabUsal/TallerGeo)#
url2 <- getURL('https://raw.githubusercontent.com/DataLabUsal/TallerGeo/master/SunshineMine.grid.csv')
grid <- read.csv(text = url2)
coordinates(grid) <- ~Easting+Elevation

##Ordinary Kriging##
au.krig <- krige(log(Au)~1, AnData, grid, model = vg.au.fit)
names(au.krig)
spplot(au.krig["var1.pred"])

#Crossvalidation#
cv.ok <- krige.cv(log(Au)~1, AnData, grid, model = vg.au.fit)
summary(cv.ok)
resid.ok <- as.data.frame(cv.ok)$residual
sqrt(mean(resid.ok^2))
mean(resid.ok)
mean(resid.ok^2/as.data.frame(cv.ok)$var1.var)

##Universal Kriging##
vg.au.width <- variogram(log(Au)~sqrt(Vein_Width), AnData)
vg.au.width.fit <- fit.variogram(vg.au.width, model = vgm(1, "Sph", 400, 1))
plot(vg.au.width,vg.au.width.fit)
au.krig.width <- krige(log(Au)~sqrt(Vein_Width), AnData, grid, model = vg.au.width.fit)
spplot(au.krig.width["var1.pred"])

#Crossvalidation#
cv.uk <- krige.cv(log(Au)~sqrt(Vein_Width), AnData, model=vg.au.width.fit)
summary(cv.uk)
resid.uk <- as.data.frame(cv.uk)$residual
sqrt(mean(resid.uk^2))
mean(resid.uk)
mean(resid.uk^2/as.data.frame(cv.uk)$var1.var)

##COKRIGING##
#Ordinary Cokriging#
gock <- gstat(NULL, "au", log(Au)~1, AnData)
gock <- gstat(gock, "ag", log(Ag)~1, AnData)


#Se puede realizar un Cokriging Universal de manera parecida a cómo se realizó en el caso del kriging#
#Universal Cokriging#
guck <- gstat(NULL, "au", log(Au)~sqrt(Vein_Width), AnData)
guck <- gstat(guck, "ag", log(Ag)~sqrt(Vein_Width), AnData)

#Variogramas Cruzados#
#Ordinary Cokriging#
cvock <- variogram(gock)
#Universal Cokriging#
cvuck <- variogram(guck)

#Ajuste del modelo inicial#
#Ordinary Cokriging#
gock <- gstat(gock, model = vgm(1, "Sph", 630, 1), fill.all = TRUE)
#Universal Cokriging#
guck <- gstat(guck, model = vgm(1, "Sph", 630, 1), fill.all = TRUE)

#Modelo Lineal de Corregionalización (LMC)
#Ordinary Cokriging#
gock.fit <- fit.lmc(cvock, gock)
plot(cvock, gock.fit)
#Universal Cokriging#
guck.fit <- fit.lmc(cvuck, guck)
plot(cvuck, guck.fit)

#Cokriging#
#Ordinary Cokriging#
cko <- predict(gock.fit, grid)
#Universal Cokriging#
cku <- predict(guck.fit, grid)

#Visualizamos los resultados
pl1 <- spplot(cko["au.pred"], main="Au Ordinary")
pl2 <- spplot(cko["ag.pred"], main="Ag Ordinary")
pl3 <- spplot(cku["au.pred"], main="Au Universal")
pl4 <- spplot(cku["ag.pred"], main="Ag Universal")

print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))

#Visualizamos la predicción del error de la varianza y covarianzas
spplot.vcov(cko)
spplot.vcov(cku)

#Crossvalidation#
#Ordinary#
cv.ck.ok <- gstat.cv(gock.fit)
summary(cv.ck.ok)
resid.cko <- as.data.frame(cv.ck.ok)$residual
sqrt(mean(resid.cko^2))
mean(resid.cko)
mean(resid.cko^2/as.data.frame(cv.ck.ok)$au.var)
##Si quieres calcular todos los residuales añade la opción "all.residuals = TRUE" en gstat.cv#
cv.ck.ok.all <- gstat.cv(gock.fit,all.residuals=TRUE)
summary(cv.ck.ok.all)


#Universal#
cv.ck.uk <- gstat.cv(guck.fit)
summary(cv.ck.uk)
resid.cku <- as.data.frame(cv.ck.uk)$residual
sqrt(mean(resid.cku^2))
mean(resid.cku)
mean(resid.cku^2/as.data.frame(cv.ck.uk)$au.var)
#Todos los residuales#
cv.ck.uk.all <- gstat.cv(guck.fit,all.residuals=TRUE)
summary(cv.ck.uk.all)