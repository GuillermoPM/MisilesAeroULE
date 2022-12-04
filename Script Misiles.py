"""
Script para los cálculos del primer ejercicio entregable de la asignatura de vehículos lanzadores y misiles.

@Autores: Guillermo Peña Martínez, Alejandro Paz Rodríguez, Raúl Ordás Collado
@Fecha: 23/11/2022
@Referencias: Missile Aerodynamics (Nielsen)
"""
#%% Módulos
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate

#%% Variables entrada

# Variables geométricas
D = 0.177				# Diámetro del misil. (m)
r0 = D/2				# Radio del misil. (m)
l = 0.5					# Longitud de la ojiva. (m)
l_misil = 3.65*3.281	# Longitud del misil ft
tr = 1					# Ancho máximo de las alas.
cr = 1                  # Cuerda superficies de control.
cdc = 1                 # Cf viscous crossflow.
m = 0.3                   # Distancia borde de ataque a punto espesor máximo.
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

# Variables de vuelo
AOA = 4 				# Ángulo de ataque (deg).
AOA = AOA*np.pi/180		# Ángulo de ataque (rad).
M = 3                 	# Mach de vuelo
h = 10000               # Altura de vuelo (m)
#%% Variables globales
B = (M**2-1)**0.5											# Parámetro de corrección de compresibilidad
dens = (1.225*(1-22.558*10 ** (-6)*h)**4.2559)*0.00194032	# Densidad a la altura de vuelo
R = 1718													# Constante gases ideales en ft^2/sec^2*R
T = 1.8*(288.15-6.5*h/1000)  								# Temperatura a la altura de vuelo en Rankine
Patm = dens*T*R
Vsound = (1.4*T*287)**0.5									# Velocidad del sonido a la altura de vuelo
Vinf = M*Vsound												# Velocidad de la corriente libre
T0 = 491.7													# Temperatura de referencia para la viscosidad en grados Rankine
mu0 = 3.58*10**(-7)											# Viscosidad a la temperatura de referencia em slug/sec*ft
dl = 50														# Divisiones del misil para las funciones de la capa límite.
#%% Datos para gráfica variaciones con la temperatura. Valores extraídos del Nielsen
Temp_Pr = [200,300,350,400,460,500,600,700,800,850,900,1000,1100,1200,1300,1400,1500,1580]
Pr_datos = [0.768,0.75,0.74,0.73,0.72,0.714,0.7,0.69,0.684,0.682,0.68,0.679,0.68,0.682,0.685,0.689,0.692,0.695]
Temp_gamma = [200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1350,1400,1500,1600,1700,1800,1900,2000]
gamma_datos = [1.41,1.408,1.405,1.402,1.401,1.4,1.398,1.397,1.392,1.39,1.388,1.385,1.382,1.38,1.378,1.375,1.372,1.37,1.368,1.366,1.365]
Temp_cp = [200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900]
cp_datos = [3.59,3.52,3.51,3.5,3.51,3.52,3.55,3.58,3.62,3.68,3.72,3.78,3.82,3.87,3.9,3.94,3.98,4.02]
Temp_mu = [200,300,400,500,550,600,650,700,750,800,900,1000,1050,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]
mu_datos = [0.2,0.63,0.82,1,1.1,1.19,1.24,1.31,1.38,1.41,1.58,1.69,1.74,1.79,1.89,2,2.1,2.18,2.24,2.34,2.4,2.51,2.58]
Pr = interpolate.CubicSpline(Temp_Pr,Pr_datos)
gamma = interpolate.CubicSpline(Temp_gamma,gamma_datos)
cp = interpolate.CubicSpline(Temp_cp,cp_datos)
mu = interpolate.CubicSpline(Temp_mu,mu_datos)

fig, ax = plt.subplots()
cp_pl = ax.twinx()
muRel_pl = ax.twinx()
Pr_pl = ax.twinx()
cp_pl.plot(Temp_cp, cp(Temp_cp), 'r-', label="Cp/R")
cp_pl.set_ylim(3, 4.5)
cp_pl.spines.right.set_position(("axes",1.2))
cp_pl.set_ylabel("Cp/R")

muRel_pl.plot(Temp_mu, mu(Temp_mu), 'g-', label="Sutherland")
muRel_pl.spines.right.set_position(("axes", 1.4))
muRel_pl.set_ylim(0, 3)
muRel_pl.set_ylabel("Sutherland")

Pr_pl.plot(Temp_Pr, Pr(Temp_Pr), 'b-', label="Prandtl")
Pr_pl.set_ylim(0.65, 0.8)
Pr_pl.spines.right.set_position(("axes", 1.6))
Pr_pl.set_ylabel("Pr")

ax.plot(Temp_gamma, gamma(Temp_gamma), "k-", label="Gamma")
ax.set_ylim(1.3, 1.45)
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
def Laminar(T,Patm):
	Titer = 1800 # Temperatura de inicio de la iteración
	Tref = Titer 
	while True:
		r = Pr(Titer)**0.5
		Ts = T*(1+(gamma(T)-1)/2*M**2)
		Trec = T + r*(Ts-T)
		Tref = T + 0.5*(Trec-T)+0.22*r*(Ts-T)
		if abs(Titer-Tref) <= 0.00001:
			break
		Titer = Tref
	dens_ref = Patm/(R*Tref)
	V0 = (gamma(T)*R*T)**0.5*M
	mu_ref = mu(Tref)*mu0
	Re_local = lambda x: V0*dens_ref/mu_ref*x
	Re_medio = V0*dens_ref/mu_ref*l_misil
	cf_medio = 0.664*2/(Re_medio**0.5)
	cf = lambda x: 0.664*Re_local(x)**(-0.5)
	x_misil = np.linspace(0.6,l_misil,dl)
	cf_local = []
	i = 0
	while Re_local(x_misil[i]) < 10**7:
		cf_local.append(cf(x_misil[i])) 
		i += 1
	delta_Blayer = 5*x_misil/Re_local(x_misil)**0.5
	displ_thickness = 1.72*x_misil/(Re_local(x_misil)**0.5)
	return cf_local,x_misil[i],delta_Blayer,dens_ref,cf_medio,displ_thickness


xT = Laminar(T, Patm)[1]				# Punto de transición

def Turbulenta(T,Patm,x):
	Titer = 1800
	Tref = Titer
	while True:
		r = Pr(Titer)**(3/2)
		Ts = T*(1+(gamma(T)-1)/2*M**2)
		Trec = T + r*(Ts-T)
		Tref = T + 0.5*(Trec-T)+0.22*r*(Ts-T)
		if abs(Titer-Tref) <= 0.00001:
			break
		Titer = Tref
	dens_ref = Patm/(R*Tref)
	V0 = (gamma(T)*R*T)**0.5*M
	mu_ref = mu(Tref)*mu0
	Re_medio = V0*dens_ref/mu_ref*l_misil
	cf_medio = 0.455/(np.log10(Re_medio)**2.58)
	x_misil = np.linspace(x,l_misil,dl)
	
	Re_local = lambda x: V0*dens_ref/mu_ref*x
	cf = lambda x: 0.370*np.log10(Re_local(x))**(-2.584)
	delta_BLayer = 0.37*x_misil*Re_local(x_misil)**(-1/5)
	cf_local = cf(x_misil)
	return cf_local,delta_BLayer,V0,dens_ref


dens_lam, dens_turb, displ_thickness = Laminar(T, Patm)[3], Turbulenta(T, Patm, xT)[3], Laminar(T, Patm)[5]/(D*3.281)
V0 = Turbulenta(T, Patm, xT)[2]
BL_laminar = np.array(Laminar(T, Patm)[0])
cf_local = np.append(BL_laminar,Turbulenta(T,Patm,xT)[0])
cf_local_func = interpolate.CubicSpline(np.linspace(0,l_misil,len(cf_local)),cf_local)
# Relación entre el espesor de la capa límite laminar y el diámetro del misil
delta_rel = Laminar(T, Patm)[2]/(D*3.281)
# Relación entre el espesor de la capa límite turbulenta y el diámetro del misil
delta_rel_turb = Turbulenta(T, Patm,xT)[1]/(D*3.281)
def Plots():
	plt.figure(1)
	plt.plot(np.linspace(0,xT,len(delta_rel)),delta_rel,label = 'Laminar')
	plt.plot(np.linspace(xT, l_misil, len(delta_rel_turb)), delta_rel_turb,label = 'Turbulenta')
	plt.plot([xT,xT],[delta_rel[-1],delta_rel_turb[0]],':')
	plt.legend()
	plt.ylabel("δ/D")
	plt.xlabel("Distancia desde el extremo de la ojiva (ft)")

	plt.figure(2)
	plt.plot(np.linspace(0,xT,len(displ_thickness)),displ_thickness)
	plt.ylabel("δ*/D")
	plt.xlabel("Distancia desde el extremo de la ojiva (ft)")

	plt.figure(3)
	plt.plot(np.linspace(0,l_misil,len(cf_local)),cf_local)
	plt.ylabel("Coeficiente de fricción local")
	plt.xlabel("Distancia desde el extremo de la ojiva (ft)")
	plt.show()
Plots()

# Coeficiente de fricción medio. Se integra la función generada como suma de los coeficientes de fricción locales de la
# parte laminar y de la parte turbulenta

cf_medio_lam = integrate.quad(cf_local_func,0,l)[0]*3**0.5 + integrate.quad(cf_local_func,l,xT)[0]
cf_medio_turb = integrate.quad(cf_local_func, xT, l_misil)[0]

Dfricc_lam, Dfricc_turb = cf_medio_lam * V0**2 / \
	2*dens_lam, cf_medio_turb * V0**2/2*dens_turb

Dfricc = Dfricc_lam + Dfricc_turb				# En unidades del sistema imperial
r_ojiva = lambda x: r0/l*x
# Integral de revolución para obtener el drag total
Dfricc_ojiva, Dfricc_fus = Dfricc_lam*2*np.pi*integrate.quad(r_ojiva, 0, l)[0], (xT- l)*Dfricc_lam*2*np.pi*r0*3.281 + 2*np.pi*r0*3.281*(l_misil-xT)*Dfricc_turb
Dfricc_misil = Dfricc_ojiva + Dfricc_fus
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

