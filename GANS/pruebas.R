library(tuneR) # leer wave
library(dplyr) # para el map
library(tidyverse) # para el map
library(signal, warn.conflicts = F, quietly = T) # signal processing functions
library(oce, warn.conflicts = F, quietly = T) # image plotting functions and nice color maps


dir_inicial = getwd()
# Lectura de controles
# leemos todos los archivos de la carpeta de controles en la subcarpeta Grabaciones


controles = map(fs::dir_ls ("Controles/Grabaciones"), ~ readWave(.x)) 
# removemos el nombre asignado por repeticion de la carpeta 
names(controles) = str_remove(names(controles) , "Controles/Grabaciones/") 

# Lectura de parkinson
# Leemos archivos de la carpeta Parkinson de la subcarpeta Grabaciones
parkinson = map(fs::dir_ls ("Parkinson/Grabaciones"), ~ readWave(.x))
#Cambiamos el nombre asignado por repeticion de la carpeta
names(parkinson) = str_remove(names(parkinson), "Parkinson/Grabaciones/")


# PLOTS 
# Numero de puntos a usar para la transformada rápida de Fourier
#nfft=1024

# Tamaño de la ventana (en puntos)
#window=256

# Superposicion  (en puntos)
#overlap=128


# PRUEBA PARA LA PRIMER GRABACION DE PACIENTE DE CONTROL RECORD1 
#spec = specgram(x = controles$`CONT2014-0002T1R1.wav`@left,
 #               n = nfft,
  #              Fs = controles$`CONT2014-0002T1R1.wav`@samp.rate,
   #             window = window,
    #            overlap = overlap)

# Descartar la informacion de la fase
#P = abs(spec$S) 
#P = P/max(P)  
#P = 10*log10(P)

# config time axis
#t = spec$t




# plot spectrogram
# Aqui cambiamos de directorio para guardar las imagenes
# LA carpeta es la asignada en varchar
#setwd(paste(getwd(),"/Imagenes_pruebas", sep = ""))

#Aqui ya se genera el plot
#png("Espectrograma1.png")
#imagep(x = t, # Eje x
#       y = spec$f, # Eje y
#       z = t(P), # Plotea el espectro, en intensidad
#       col = oce.colorsViridis, # Color
#       drawPalette = F, # Dice si la paleta se muestra o no
#       decimate = F 
#)
#dev.off()

#AQUI TERMINA LA PRUEBA, SALIO BIEN!!



# ASIGNAMOS LOS NOMBRES QUE CONTENDRAN LOS PLOTS 
nombres = paste("Espectrograma", 1:length(controles), ".png", sep = "")
#CReamos listas para ir guardando las iteraciones
tasa_muestras = espectrogramas = seniales = p = t = f = list(list())
#Reasignamos el directorio de controles donde queremos que se guarden los espectros
dir_cont_espec = paste(getwd(),"/Controles/Espectrogramas", sep = "")

#Aqui se crea la funcion
f_creacion_espectros = function(list_grabs,nfft,window,overlap,nombres){
        for(i in 1:length(list_grabs)){
                seniales[[i]] = list_grabs[[i]]@left
                tasa_muestras[[i]] = list_grabs[[i]]@samp.rate
                espectrogramas[[i]]  = specgram(x = seniales[[i]],
                                                n = nfft,
                                                Fs = tasa_muestras[[i]],
                                                window = window,
                                                overlap = overlap)
                p[[i]] = abs(espectrogramas[[i]]$S)  
                p[[i]] = p[[i]]/max(p[[i]]) 
                p[[i]] = 10*log10(p[[i]])
                t[[i]] = espectrogramas[[i]]$t
                f[[i]] = espectrogramas[[i]]$f
                png(nombres[i])
                imagep(x = t[[i]],
                       y = f[[i]],
                       z = t(p[[i]]),
                       col = oce.colorsViridis,
                       ylab = 'Frequency [Hz]',
                       xlab = 'Time [s]',
                       drawPalette = F,
                       decimate = F)
                dev.off()
                }
}

#Aqui se manda llamar la funcion para los controles
f_creacion_espectros(controles,1024,256,128,nombres)

#Ahora aplicaremos la funcion para los pacientes con parkinson

#Cambiamos o volvemos al directorio inicial

setwd(dir_inicial)
#Ahora nos movemos a la carpeta donde queremos que se guarden los espectros de parkinson
dir_park_espec = paste(getwd(),"/Parkinson/Espectrogramas", sep = "")
#MAndamos llamar la funcion
f_creacion_espectros(parkinson,1024,256,128,nombres)
#Regresamos al directorio inicial 
setwd(dir_inicial)


# CREACION DEL ALGORITMO 
# GANS



#Preparacion de los datos
library(keras) # Mnist
library(EBImage) #mnist
library(png) # para usar la funcion readPNG (hasta ahora no empleada)
library(imager) # para la funcion load.image 
library(caret)
library(caTools)
###########


##GRUPO DE CONTROL
#Listamos los files de la carpeta donde estan guardados los espectros del grupo de control
imag_control = list.files(path = dir_cont_espec, pattern="*.png", full.names = TRUE)
#Aqui leemos las imagenes y las guardamos en una lista 
all_im_control <- lapply(imag_control, load.image )

##GRUPO DE PARKINSON
#Listamos los files de la carpeta donde estan guardados los espectros del grupo de parkinson
imag_parkinson = list.files(path = dir_park_espec, pattern="*.png", full.names = TRUE)
#Aqui leemos las imagenes y las guardamos en una lista 
all_im_parkinson <- lapply(imag_parkinson, load.image )

#PREUBA PLOT"
# Agrupamos las imagenes y creamos el plot
#par(mfrow = c(4,4), mar = rep(0,4))
#for (i in 1:16) plot(as.raster(all_im_control[[i]]))
#par(mfrow = c(1,1))
## FIN PRUEBA PLOT

# Creamos las muestras de test y de entrenamiento para el modelo de grupo de control
sample_control<-sample.split(all_im_control, SplitRatio = 0.75)
train_control = subset(all_im_control, sample == TRUE)
test_control  = subset(all_im_control, sample == FALSE)


# Creamos las muestras de test y de entrenamiento para el modelo de grupo de parkinson
sample_parkinson<-sample.split(all_im_parkinson, SplitRatio = 0.75)
train_parkinson = subset(all_im_parkinson, sample == TRUE)
test_parkinson  = subset(all_im_parkinson, sample == FALSE)



conjunta = list(train = list(x = train_control, y = train_parkinson) , 
                test = list(x = test_control , y = test_parkinson)) 


str(conjunta_control)

c(c(traincont, trainpark), c(testcont,testpark)) %<-% conjunta

for (i in length(traincont)){
  print(dim(traincont[[i]]))
 reshape_c = array_reshape(traincont[[i]] , c(length(traincont), 28,28,1))
  print(dim(reshape_c))
}

str(traincont[[1]])


vector = array()
for(i in length(traincont)){
  vector = unlist(traincont[[i]])
  vector = append(vector)
  return(vector)
}

unlist_traincont = unlist(traincont[[1]])
str(unlist_traincont)

getwd()

setwd(dir_cont_espec)

im1 = readPNG("Espectrograma1.png")
m = as.matrix(im1)

nrow(unlist_traincont)
retrainx 




#PRUEBA 
mnist = dataset_mnist()
str(mnist)

str(all_im_control)

c(c(trainx,trainy), c(testx,testy)) %<-% mnist

trainx = array_reshape(trainx,
                       c(nrow(trainx), 28,28,1))

trainx = trainx / 255

summary(trainx)
#Red de generadores



#Red de discriminadores



















