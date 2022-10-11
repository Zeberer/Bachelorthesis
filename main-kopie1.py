import matplotlib.pyplot as plt
import numpy as np

daten = np.genfromtxt(r"Gerade.csv", float, delimiter=" ")


# ax = plt.axes(projection='3d')

# Kreis mit Radius 2

def Kreis(t):
    x = 2 * np.cos(2 * np.pi * t)
    y = 2 * np.sin(2 * np.pi * t)
    return x, y


x, y = Kreis(np.linspace(0, 1, 100))


f = open("kartesische_Koordinaten2.csv", "a")

# Data for a three-dimensional line

xline = daten[:, 0] * np.sin(daten[:, 1]) * np.cos(daten[:, 2])  # daten[:, 0] = r, daten[:,1] = theta, daten[:,2] = phi
yline = daten[:, 0] * np.sin(daten[:, 1]) * np.sin(daten[:, 2])
zline = daten[:, 0] * np.cos(daten[:, 1])

for i in range(0, len(xline)):
    f.write(str(xline[i]) + ' ' + str(yline[i]) + ' ' + str(zline[i]) + '\n')

f.close()

"""
ax.plot3D(xline, yline, zline, 'gray', label=r'dr=-0.5, $d\phi$=0.0583, $d\theta$=0.03')
ax.plot_surface(xK, yK, zK)
ax.set_xlabel('$x[solar masses]$')
ax.set_ylabel('$y[solar masses]$')
ax.set_zlabel('$z[solar masses]$')
ax.set_title("3D plot of a light curve")
#plt.legend(loc='upper left')
plt.savefig(r"C:\\Users\\SoZeb\\Documents\\Uni\\Bachelorarbeit\\Plots\\3D_3d.pdf")
plt.show()
"""


r_min = min(daten[:, 0])
print(r_min)
r_min = r_min * 1480


G = 6.67 * 10**(-11)
M = 2 * 10**30
c = 3 * 10**8

alpha_b = (4 * G * M) / (r_min * (c**2))
print('alpha_b = ', alpha_b)

alpha_bd = alpha_b * 180 / np.pi
print('alpha_bd = ', alpha_bd)

alpha_bs = alpha_b * (3600 * 180) / np.pi
print('alpha_bs = ', alpha_bs)





AV = [xline[1] - xline[0], yline[1] - yline[0], zline[1] - zline[0]]
EV = [xline[-1] - xline[-2], yline[-1] - yline[-2],
      zline[-1] - zline[-2]]  # -1 ist letzter Eintrag in Array, -2 sind letzte beide Einträge

print(AV)
print(EV)

kreuz = np.cross(AV, EV)

alpha = np.arctan(np.cross(AV, EV) / np.dot(AV, EV))

print('alpha = ', alpha)

alpha = np.sqrt(alpha[0] ** 2 + alpha[1] ** 2 + alpha[2] ** 2)

# alpha = np.arctan((np.sqrt(kreuz[0]**2 + kreuz[1]**2 + kreuz[2]**2)) / np.dot(AV, EV))

print('alpha=', alpha)

grad = alpha * 180 / np.pi

print('grad=', grad)

winkelsekunden = alpha * (3600 * 180) / np.pi

print('winkelsekunden=', winkelsekunden)

"""
plt.plot(xline, yline, label="x=10, y=0, z=0, dx=1.3, dy=0, dz=0")
plt.plot(x, y, color='black', linestyle='--')
plt.legend(loc="lower right")
plt.xlabel(r'x in $M_\odot$', fontsize=13)
plt.ylabel(r'y in $M_\odot$', fontsize=13)
plt.title('radial component', fontsize=18)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.savefig(r"C:\\Users\\SoZeb\\Documents\\Uni\\Bachelor\\Bachelorarbeit\\Plots\\Gerade.pdf")
plt.show()
"""

plt.plot(xline, yline)
plt.plot(x, y, color='black', linestyle='--')
plt.xlabel(r'x in $M_\odot$', fontsize=13)
plt.ylabel(r'y in $M_\odot$', fontsize=13)
plt.title('ISCO', fontsize=18)
plt.savefig(r"C:\\Users\\SoZeb\\Documents\\Uni\\Bachelor\\Bachelorarbeit\\Plots\\ISCO1.pdf")
plt.show()

