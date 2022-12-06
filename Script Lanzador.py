"""
Script para los cálculos del segundo ejercicio entregable de la asignatura de Vehículos Lanzadores y Misiles.

@Autores: Guillermo Peña Martínez, Alejandro Paz Rodríguez, Raúl Ordás Collado
@Fecha: 2/12/2022
@Referencias: 
"""

import numpy as np
import pandas as pd

datos_etapas = {
	"impulso especifico": [315, 330, 345],
	"fraccion estructural": [0.10, 0.15, 0.20]
}

etapas = pd.DataFrame(datos_etapas, index=["etapa 1", "etapa 2", "etapa 3"])

g0 = 0.00981	# Aceleración gravitatoria (m/s^2)
mf = 2000		# Carga útil (kg)
vbo = 7.8		# Velocidad en el apagado (km/s)

Is = etapas["impulso especifico"].mean()
ε = etapas["fraccion estructural"].mean()


""" Primera parte """

πPL = (np.exp((-vbo)/(Is*g0)) - ε)/(1-ε)
print(πPL)
