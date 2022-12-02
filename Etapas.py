import numpy as np

#-------------EJERCICIO 2----------------------

mpl = 2000
v = 7800
Is = 330
epsilon = 0.15
g = 9.81

n = np.exp(v/(3*Is*g))
pi_pl = ((1/n-epsilon)/(1-epsilon))**3

# Masas totales

m0 = mpl/pi_pl
m02 = mpl/pi_pl**(2/3)
m03 = mpl/pi_pl**(1/3)
print("Masas de las etapas")
print(m0,m02,m03)

# Masas en vacio

mE1 = epsilon*(m0-m02)
mE2 = epsilon*(m02-m03)
mE3 = epsilon*(m03-mpl)
print("Masas en vacio")
print(mE1,mE2,mE3)

# Masas de propelente

mP1 = m0-mE1-m02
mP2 = m02-mE2-m03
mP3 = m03-mE3-mpl
print("Masas de propelente")
print(mP1,mP2,mP3)