import matplotlib.pyplot as plt
import numpy as np


def f(y, t):
    return np.cos(t)


y = 0
t = 0
h = 0.06
dtend = 100

y1werte = []
t1werte = []

y1werte.append(y)
t1werte.append(t)

while(t<dtend):

    k1 = f(y, t)
    k2 = f(y + h*k1/2, t + h/2)
    k3 = f(y + h*k2/2, t + h/2)
    k4 = f(y + h*k3, t + h)

    y = y + h/6*(k1 + 2*k2 + 2*k3 + k4)
    t = t + h

    y1werte.append(y)
    t1werte.append(t)


print(y1werte)
print(t1werte)

#print(np.array_repr(t).replace("\n", ""))
#print(np.array_repr(y).replace("\n", ""))
np.savetxt('t1RK4_0.06.csv', t1werte, delimiter='')
np.savetxt('y1RK4_0.06.csv', y1werte, delimiter='')

plt.plot(t1werte, y1werte)
plt.xlabel('t')
plt.ylabel('sin(t)')
plt.title('Sinus RK4 0.06')
plt.savefig(r"C:\Users\SoZeb\PycharmProjects\pythonProject1\RK4_0.06.pdf")
plt.show()






