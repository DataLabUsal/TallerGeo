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

#A partir de aquí toca trabajar. Utiliza el archivo TallerGeo.R para ver cómo puedes continuar.