"""
Script para los cálculos del segundo ejercicio entregable de la asignatura de Vehículos Lanzadores y Misiles.

@Autores: Guillermo Peña Martínez, Alejandro Paz Rodríguez, Raúl Ordás Collado
@Fecha: 2/12/2022
@Referencias: 
"""

from numpy import exp, mean, linspace
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator

etapas = {
	"Is": [315, 330, 345],
	"ε": [0.10, 0.15, 0.20]
}

def figure_features(tex=True, font="serif", dpi=180):
    """Customize figure settings.
    Args:
        tex (bool, optional): use LaTeX. Defaults to True.
        font (str, optional): font type. Defaults to "serif".
        dpi (int, optional): dots per inch. Defaults to 180.
    """
    plt.rcParams.update(
        {
            "font.size": 16,
            "font.family": font,
			"font.serif": ['Palatino'],
            "text.usetex": tex,
			"text.latex.preamble": r'\usepackage{amsmath}',
            "figure.subplot.top": 0.9,
            "figure.subplot.right": 0.9,
            "figure.subplot.left": 0.15,
            "figure.subplot.bottom": 0.12,
            "figure.subplot.hspace": 0.4,
            "savefig.dpi": dpi,
            "savefig.format": "png",
            "axes.titlesize": 16,
            "axes.labelsize": 18,
            "axes.axisbelow": True,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 5,
            "xtick.minor.size": 2.25,
            "xtick.major.pad": 7.5,
            "xtick.minor.pad": 7.5,
            "ytick.major.pad": 7.5,
            "ytick.minor.pad": 7.5,
            "ytick.major.size": 5,
            "ytick.minor.size": 2.25,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 16,
            "legend.framealpha": 1,
            "figure.titlesize": 16,
            "lines.linewidth": 2,
        }
    )

def add_grid(ax, lines=True, locations=None):
    """Add a grid to the current plot.
    Args:
        ax (Axis): axis object in which to draw the grid.
        lines (bool, optional): add lines to the grid. Defaults to True.
        locations (tuple, optional):
            (xminor, xmajor, yminor, ymajor). Defaults to None.
    """

    if lines:
        ax.grid(lines, alpha=0.5, which="minor", ls=":")
        ax.grid(lines, alpha=0.7, which="major")

    if locations is not None:

        assert (
            len(locations) == 4
        ), "Invalid entry for the locations of the markers"

        xmin, xmaj, ymin, ymaj = locations

        ax.xaxis.set_minor_locator(MultipleLocator(xmin))
        ax.xaxis.set_major_locator(MultipleLocator(xmaj))
        ax.yaxis.set_minor_locator(MultipleLocator(ymin))
        ax.yaxis.set_major_locator(MultipleLocator(ymaj))

g0 = 0.00981	# Aceleración gravitatoria (m/s^2)
mf = 2000		# Carga útil (kg)
vbo = 7.8		# Velocidad en el apagado (km/s)

πPL = lambda Is, ε : (exp((-vbo)/(Is*g0)) - ε)/(1-ε)

# Cálculo de la fracción de carga con los datos medioss
print(πPL(mean(etapas["Is"]), mean(etapas["ε"])))

Is = linspace(315, 550, 100)
ε = etapas["ε"]

figure_features()

fig = plt.figure(figsize=(8, 6))
ax = plt.axes(xlim=(315, 550))

ax.set_xlabel("$I_{sp}$")
ax.set_ylabel("$\\pi_{\\text{PL}}$", labelpad=12.0)

add_grid(ax, locations=(25, 50, 0.05, 0.1))

ax.axvline(x=345.31, ymax=0.5, color='grey', linestyle='--')
ax.axvline(x=419.11, ymax=0.5, color='grey', linestyle='--')
ax.axvline(x=494.027, ymax=0.5, color='grey', linestyle='--')
ax.plot(Is, πPL(Is, ε[0]), color="coral", label="$\\varepsilon = 0.10$")
ax.plot(Is, πPL(Is, ε[1]), color="mediumseagreen", label="$\\varepsilon = 0.15$")
ax.plot(Is, πPL(Is, ε[2]), color="royalblue", label="$\\varepsilon = 0.20$")
ax.plot(345.31, 0, 'o', color='red')
ax.plot(419.11, 0, 'o', color='red')
ax.plot(494.027, 0, 'o', color='red')


ax.legend(loc='lower right')
plt.show()