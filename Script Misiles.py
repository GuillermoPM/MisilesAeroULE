"""
Script para los cálculos del primer ejercicio entregable de la asignatura de vehículos lanzadores y misiles.

@Autores: Guillermo Peña Martínez, Alejandro Paz Rodríguez, Raúl Ordás Collado
@Fecha: 23/11/2022
"""
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols

#%% Variables entrada
# Variables geométricas
D = 0.177				# Diámetro del misil.
r0 = D/2				# Radio del misil.
l = 0.5					# Longitud de la ojiva.
tr = 1                  # Ancho máximo de las alas.
cr = 1                  # Cuerda superficies de control.
cdc = 1                 # Cf viscous crossflow.
m = 1                   # Distancia borde de ataque a punto espesor máximo.
Bc = 0.525              # Wingspan
deflx_w = 0.6           # Deflexión de la estela, término (1-dew/dalpha)
Sm = 10                 # Superficie proyectada del misil
Scontrol = 10           # Superficie alar de los controles (m^2)
Sw = 10                 # Superficie alar (m^2)
Bw = 0.56               # Wingspan ala (m)
Xcg = 0.2               # Distancia al centro de gravedad
Xcp = 0.1               # Distancia al centro de presiones
Xw = 0.3                # Distancia adimensional al centro de presiones del ala
Xb = 0.4                # Distancia adimensional al centro de presiones del fuselaje
Xc = 0.5                # Distancia adimensional al centro de presiones del control
Iy = 1                  # Momento de inercia del misil respecto a un eje perpendicular al plano del movimiento

# Características del misil
masa = 1                # Masa del misil

# Variables de vuelo
AOA = 4 				# Ángulo de ataque (deg).
AOA = AOA*np.pi/180		# Ángulo de ataque (rad).
M = 2.3                 # Mach de vuelo
h = 10000               # Altura de vuelo (m)
#%% Variables globales
B = (M**2-1)**0.5
dens = 1.225*(1-22.558*10 **(-6)*h)**4.2559
Temp = 288.15-6.5*h/1000
Vsound = (1.4*Temp*287)**0.5
Vinf = M*Vsound
q = 0.5*Vinf**2*dens

#%% Cálculo de la resistencia
# El drag debido al lift en un cuerpo esbelto de base cilíndrica
# es la mitad del que se produce en una placa plana (Nielsen).
L_q = 2*np.pi*AOA*r0**2
D_q = 0.5*L_q*AOA                                   # Resitencia debida a la sustentación / presión dinámica

# Viscous Crossflow drag:
Sc = np.pi*r0**2                                    # Superficie transversal
Dc_q = cdc*AOA**3*Sc                                # Viscous crossflow drag / presión dinámica

# Drag frontal en las alas para AOA = 0
Cd0 = 1/(4*m*(1-m))*4*(tr/cr)**2/B

# Drag de fricción






#%% Fuerzas en el giro configuración "clásica"
# Fuerzas debidas a alpha
CNic = 1                                            # Pendiente del coeficiente de fuerza normal del control aislado
CNiw = 1                                            # Pendiente del coeficiente de fuerza normal del ala aislada
Cnalpha_b = 2                                       # Coeficiente de fuerza normal debido a alpha del fuselaje
Cnalpha_c = (1+D/Bc)*deflx_w*Scontrol/Sm*CNic       # Coeficiente de fuerza normal debido a alpha del control
Cnalpha_w = CNiw*Sw/Sm*(1+D/Bw)                     # Coeficiente de fuerza normal debido a alpha de las alas
CN_alpha = Cnalpha_b + Cnalpha_c + Cnalpha_w        # Coeficiente de fuerza normal debido a alpha

# Fuerzas debidas a delta
Cndelta_c = CNic*Scontrol/Sm*(1+D/Bc)               # Coeficiente de fuerza normal debido a delta
Cn = CN_alpha + Cndelta_c                           # Coeficiente de fuerza normal total

#%% Momentos en el giro
# Momentos debidos a alpha y delta
Margen_estatico = (Cnalpha_b*(Xcg-Xb)+Cnalpha_c*(Xcg-Xc)+Cnalpha_w*(Xcg-Xw))/(Cnalpha_b+Cnalpha_c+Cnalpha_w)
Cm_alpha = CN_alpha * Margen_estatico               # Coeficiente de momentos debido al ángulo de ataque
Cm_delta = Cndelta_c * (Xcg-Xc)                     # Coeficiente de momentos debido al ángulo del control

# Momentos debidos a la velocidad angular y la variación del AOA
Cm_ang = -2*(Cnalpha_b*(Xcg-Xb)**2+Cnalpha_c*(Xcg-Xc)**2+Cnalpha_w*(Xcg-Xw)**2)
# Cmf_ang = gasto_m*(Iy/masa-re**2)/(0.5/Vinf*q*S*d**2)
Cm_variacion_AOA = -2*CNic*(Sc/Sm*r0)*(1+D/Bc)*deflx_w*(Xw-Xc)*(Xw-Xcg)             # Coeficiente de momentos debido a la variación del AOA
