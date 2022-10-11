import matplotlib.pyplot as plt
import numpy as np


daten_RK4_006 = np.genfromtxt(r"y1RK4_0.06.csv", float, delimiter=" ")
daten_ana_006 = np.genfromtxt(r"y1_0.06.csv", float, delimiter=" ")
daten_t_006 = np.genfromtxt(r"t1RK4_0.06.csv", float, delimiter=" ")

daten_RK4_003 = np.genfromtxt(r"y1RK4_0.03.csv", float, delimiter=" ")
daten_ana_003 = np.genfromtxt(r"y1_0.03.csv", float, delimiter=" ")
daten_t_003 = np.genfromtxt(r"t1RK4_0.03.csv", float, delimiter=" ")



delta_y_006 = np.subtract(daten_ana_006, daten_RK4_006)
print('delta_y_006 = ', delta_y_006)


delta_y_003 = np.subtract(daten_ana_003, daten_RK4_003)
print('delta_y_003 = ', delta_y_003)

delta_ana = np.subtract(daten_ana_003, daten_ana_003)


print('delta_y_006_t10 = ', delta_y_006[168])
print('delta_y_003_t10 = ', delta_y_003[335])


delta_y_003_richtig = delta_y_006[168] * 1/16
print('delta_y_003_richtig = ', delta_y_003_richtig)


delta_y_skaliert = delta_y_003 * 16
print('delta_y_003_skaliert = ', delta_y_skaliert)


plt.plot(daten_t_006, delta_y_006)
plt.xlabel('t')
plt.ylabel('delta_y')
plt.title('Fehler/differenz_0.06')
plt.savefig(r"C:\Users\SoZeb\PycharmProjects\pythonProject1\Fehler_Differenz_0.06.pdf")
plt.show()


plt.plot(daten_t_003, delta_y_003)
plt.xlabel('t')
plt.ylabel('Delta_y')
plt.title('Fehler/differenz_0.03')
plt.savefig(r"C:\Users\SoZeb\PycharmProjects\pythonProject1\Fehler_Differenz_0.03.pdf")
plt.show()



plt.plot(daten_t_006, delta_y_006, label="h = 0.06")
plt.plot(daten_t_003, delta_y_003, label="h = 0.03")
plt.plot(daten_t_003, delta_ana, label="analytic solution")
plt.legend(loc="upper left")
plt.xlabel('t', fontsize=13)
plt.ylabel(r'$\Delta y$', fontsize=13)
plt.title('convergence test', fontsize=18)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.savefig(r"C:\Users\SoZeb\PycharmProjects\pythonProject1\convergence test.png")
plt.show()


plt.plot(daten_t_006, delta_y_006, label="h = 0.06")
plt.plot(daten_t_003, delta_y_skaliert, color='orange', linewidth=2, linestyle=':', label="h = 0.03 scaled")
plt.plot(daten_t_003, delta_ana, color='green', label="analytic solution")
plt.legend(loc="upper left")
plt.xlabel('t', fontsize=13)
plt.ylabel(r'$\Delta y$', fontsize=13)
plt.title('convergence test with scaling', fontsize=18)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.savefig(r"C:\Users\SoZeb\PycharmProjects\pythonProject1\scaling.png")
plt.show()
