import matplotlib.pyplot as plt
import numpy as np



t = np.linspace(0, 100.02, 1668)
y = np.sin(t)

print(np.array_repr(t).replace("\n", ""))
print(np.array_repr(y).replace("\n", ""))
np.savetxt('t1_0.06.csv', t, delimiter='')
np.savetxt('y1_0.06.csv', y, delimiter='')

plt.plot(t, y)
plt.xlabel('t')
plt.ylabel('sin(t)')
plt.title('Sinus analytisch')
plt.savefig(r"C:\Users\SoZeb\PycharmProjects\pythonProject1\Sinusana_0.06.pdf")
plt.show()






