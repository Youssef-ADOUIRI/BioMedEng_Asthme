# BioMedEng_Asthme

Ce dépôt contient les codes MATLAB utilises dans le Projet Biomedical Engineering pour la simulation de l’écoulement de l’air dans les bronches asthmatiques.

# Prérequis



MATLAB

# Fichiers



## Bronche saine

### [Tracer_sans_obstacle.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/Tracer_sans_obstacle.m)

Permet de tracer les profils de vitesse et de pression pour un tube droit, et affiche l’erreur calculée dans les tests de validation. Nécessite les fichiers “laplace2d_General_v0.m” et “Validate_StokesEq.m”.

### [laplace2d_General_v0.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/laplace2d_General_v0.m)

Calculer l’inverse de la matrice A dans le cas d’un tube droit.

### [Validate_StokesEq.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/Validate_StokesEq.m)

Tests de validation de la vitesse et la pression.

## Bronche asthmatique

### [Tracer_avec_obstacle.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/Tracer_avec_obstacle.m)

Permet de tracer les profils de vitesse et de pression pour un tube déformé avec une seule déformation, ainsi que le débit suivant la longueur du tube, et affiche l’erreur calculée dans les tests de validation. Nécessite les fichiers “laplace2d_General_v1.m” et “Validate_StokesEq_obs.m”.

### [laplace2d_General_v1.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/laplace2d_General_v1.m)

Calculer l’inverse de la matrice A dans le cas d’une seule déformation.

### [Validate_StokesEq_obs.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/Validate_StokesEq_obs.m)

Tests de validation de la vitesse et la pression, ne prenant pas en compte les points a l’intérieur de l’obstacle.

### [Tracer_obs2.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/Tracer_obs2.m)

Permet de tracer les profiles de vitesse et de pression pour un tube déformé avec deux déformations. Nécessite le fichier “laplace2d_General_v4.m”.

### [laplace2d_General_v4.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/laplace2d_General_v4.m)

Calculer l’inverse de la matrice A dans le cas de deux déformations.

### [Reistance_long.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/Reistance_long.m)

Tracer la résistance en fonction de la longueur de la déformation.

### [Reistance_ord.m](https://github.com/Youssef-ADOUIRI/BioMedEng_Asthme/blob/master/Reistance_ord.m)

Tracer la résistance en fonction de la hauteur de la déformation.
