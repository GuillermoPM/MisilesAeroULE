"""
Script para los cálculos del primer ejercicio entregable de la asignatura de vehículos lanzadores y misiles.

@Autores: Guillermo Peña Martínez, Alejandro Paz Rodríguez, Raúl Ordás Collado
@Fecha: 23/11/2022
@Referencias: Missile Aerodynamics (Nielsen)
"""
#%% Módulos
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
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
B = (M**2-1)**0.5										# Parámetro de corrección de compresibilidad
dens = 1.225*(1-22.558*10 **(-6)*h)**4.2559				# Densidad a la altura de vuelo
T = 288.15-6.5*h/1000									# Temperatura a la altura de vuelo
Vsound = (1.4*T*287)**0.5								# Velocidad del sonido a la altura de vuelo
Vinf = M*Vsound											# Velocidad de la corriente libre
q = 0.5*Vinf**2*dens									# Presión dinámica
gamma0 = 1.403											# Relación de calores específicos
theta = 5500											# Constante en grados Rankine
R = 1718												# Constante gases ideales en ft^2/sec^2*R
cp0 = gamma0/(gamma0-1)*R								# Cp para gas calóricamente perfecto
cv0 = cp0/gamma0											# Cv para gas calóricamente perfecto
T0 = 491.7												# Temperatura de referencia para la viscosidad en grados Rankine
mu0 = 3.58*10**(-7)										# Viscosidad a la temperatura de referencia em slug/sec*ft
#%% Datos para gráfica variaciones con la temperatura. Valores extraídos del Nielsen
Temp_Pr = [200,300,350,400,460,500,600,700,800,850,900,1000,1100,1200,1300,1400,1500,1580]
Pr_datos = [0.768,0.75,0.74,0.73,0.72,0.714,0.7,0.69,0.684,0.682,0.68,0.679,0.68,0.682,0.685,0.689,0.692,0.695]
Temp_gamma = [200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1350,1400,1500]
gamma_datos = [0.76,0.758,0.755,0.752,0.751,0.75,0.749,0.747,0.744,0.74,0.738,0.734,0.732,0.72,0.728,0.725]
Pr = interpolate.CubicSpline(Temp_Pr,Pr_datos)
plt.plot(Temp_Pr,Pr(Temp_Pr))
plt.plot(Temp_Pr,Pr_datos)
plt.show()


#%% Gráfica variaciones con la temperatura
Temp = np.linspace(0,2000,100)
Temp_gm = np.linspace(200,2000,100)
mu_rel = lambda Tst: (Tst/T0)**(3/2)*(T0+198.72)/(Tst+198.72)       # Ley de sutherland para temperatura en grados Rankine
gamma = lambda Tst: 1+(gamma0-1)/(1+(gamma0-1)*((theta/Tst)**2*np.exp(theta/Tst)/(np.exp(theta/Tst)-1)**2))
cp = lambda Tst: cp0*(1+(gamma0-1)/gamma0*((theta/Tst)**2*np.exp(theta/Tst)/(np.exp(theta/Tst)-1)**2))/R
cv = lambda Tst: cv0*(1+(gamma0-1)*((theta/Tst)**2*np.exp(theta/Tst)/(np.exp(theta/Tst)-1)**2))/R
k = lambda Tst: 5.75*10**(-5)*(1+0.00317*Tst-0.0000021*Tst**2)*0.5781759824/460.67
# k = lambda Tst: 0.25*(9*gamma(Tst)-5)*mu_rel(Tst)*mu0*cv(Tst)*R
Pr = lambda Tst: cp(Tst)*R*mu_rel(Tst)*mu0/k(Tst)


fig, ax = plt.subplots()
cp_pl = ax.twinx()
muRel_pl = ax.twinx()
Pr_pl = ax.twinx()
cp_pl.plot(Temp,cp(Temp_gm),'r-',label="Cp/R")
cp_pl.set_ylim(3,4.5)
cp_pl.set_ylabel("Cp/R")

muRel_pl.plot(Temp,mu_rel(Temp),'g-',label="Sutherland")
muRel_pl.spines.right.set_position(("axes",1.2))
muRel_pl.set_ylim(0,3)
muRel_pl.set_ylabel("Sutherland")

Pr_pl.plot(Temp, Pr(Temp_gm), 'b-', label="Prandtl")
Pr_pl.set_ylim(0.65,0.8)
Pr_pl.spines.right.set_position(("axes",1.4))
Pr_pl.set_ylabel("Pr")

ax.plot(Temp,gamma(Temp_gm),"k-",label="Gamma")
ax.set_ylim(1.3,1.45)
ax.set_ylabel("Gamma")
ax.grid(True)

plt.show()
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
