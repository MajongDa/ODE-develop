import matplotlib.pyplot as plt
import numpy as np

d_path = '/home/majong/Documents/simfor/build/examples/ODUtest/'

data = np.loadtxt(d_path + 'odu_1_RK.txt')
x = data[0, :]
y = data[1, :]
y_appr = data[2, :]
err = data[3, :]

fig, (ax1, ax2) = plt.subplots(2, sharex=True)

ax1.plot(x, y, label='Аналит. решение')
ax1.plot(x, y_appr, 'x', label='Метод Рунге-Кутта')
ax1.set_ylabel('y', fontsize=16)

ax2.plot(x, err, label='Метод Рунге-Кутта')


data = np.loadtxt(d_path + 'odu_1_AdBash.txt')
x = data[0, :]
y = data[1, :]
y_appr = data[2, :]
err = data[3, :]
ax1.plot(x, y_appr, 'o', label='Метод Адамса-Башфорта')
ax1.set_ylabel('y', fontsize=16)

ax2.plot(x, err, label='Метод Адамса-Башфорта')

data = np.loadtxt(d_path + 'odu_1_AdMltn.txt')
x = data[0, :]
y = data[1, :]
y_appr = data[2, :]
err = data[3, :]


ax1.plot(x, y_appr, 'v', label='Метод Адамса-Мултона')
ax1.set_ylabel('y', fontsize=16)

ax2.plot(x, err, label='Метод Адамса-Мултона')
ax2.set_ylabel('Абс. ошибка', fontsize=16)



plt.xlabel('x', fontsize=16)
ax1.legend()
ax2.legend()
fig.tight_layout()
plt.show()
