####TALLER GEOESTADISTICA

####PAQUETES PARA LA INSTALACIÓN

install.packages('sp')		#Paquete específico para la representación espacial
install.packages('gstat')	#Paquete más usado para análisis geoestadístico en R

library(sp)
library(gstat)

#Explorando los datos
data(meuse)

head(meuse)		#Encabezado de los datos
summary(meuse)	#Resumen de los datos

#Nuestro principal objetivo es ver el comportamiento de las variables Zinc, Cadmio, Cobre y Plomo.
AnData <- meuse[,1:6] 
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

##Kriging##
zn.krig <- krige(zinc~1, AnData, meuse.grid, model = vg.zn.fit)
names(zn.krig)
p1 <- spplot(zn.krig["var1.pred"])
p2 <- spplot(zn.krig["var1.var"])
print(p1, split=c(1,1,2,1), more=TRUE)
print(p2, split=c(2,1,2,1), more=TRUE)

##Cokriging##
g <- gstat(NULL, "zn", zinc~1, AnData)
g <- gstat(g, "cd", cadmium~1, AnData)
g <- gstat(g, "cu", copper~1, AnData)
g <- gstat(g, "pb", lead~1, AnData)

#Variogramas Cruzados
cv <- variogram(g)

#Ajuste del modelo inicial
g <- gstat(g, model = vgm(1, "Sph", 900, 1), fill.all = TRUE)

#Modelo Lineal de Corregionalización (LMC)
g.fit <- fit.lmc(cv, g)
plot(cv, g.fit)

#Cokriging
ck <- predict(g.fit, meuse.grid)

#Visualizamos los resultados
pl1 <- spplot(ck["zn.pred"], main="Zinc")
pl2 <- spplot(ck["cu.pred"], main="Cobre")
pl3 <- spplot(ck["cd.pred"], main="Cadmio")
pl4 <- spplot(ck["pb.pred"], main="Plomo")

print(pl1, split = c(1,1,2,2), more=TRUE)
print(pl2, split = c(1,2,2,2), more=TRUE)
print(pl3, split = c(2,1,2,2), more=TRUE)
print(pl4, split = c(2,2,2,2))
