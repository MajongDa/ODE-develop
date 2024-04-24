import numpy  as np
import matplotlib.pyplot as plt

d_path = '/home/majong/Documents/simfor/build/examples/ODUtest/Osc/'

mu0_1 = np.loadtxt ( d_path + 'VdP_osc_mu=0.1.txt' )
mu_1 = np.loadtxt ( d_path + 'VdP_osc_mu=1.txt' )
mu_5 = np.loadtxt ( d_path + 'VdP_osc_mu=5.txt')
mu_10 = np.loadtxt ( d_path + 'VdP_osc_mu=10.txt' )

plt.figure(figsize=(16,9))

# Lighten borders
plt.gca().spines["top"].set_alpha(0)
plt.gca().spines["bottom"].set_alpha(.3)
plt.gca().spines["right"].set_alpha(0)
plt.gca().spines["left"].set_alpha(.3)

plt.title('Изменение формы предельного цикла при изменении μ', fontsize=22)
plt.autoscale(enable=None, axis="x", tight=True)
plt.ylim(-2.5,2.5)
plt.xlim(-15,15)
plt.grid(alpha=0.3)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.ylabel('u', fontsize=22)
plt.xlabel('v', fontsize=22)

plt.plot(mu0_1[:,1], mu0_1[:,0], label='μ=0.5')
plt.plot(mu_1[:,1], mu_1[:,0], label='μ=1')
plt.plot(mu_5[:,1], mu_5[:,0], label='μ=5')
plt.plot(mu_10[:,1], mu_10[:,0], label='μ=10')

plt.autoscale(enable=False, axis='both', tight=None)
plt.legend(fontsize=22)
plt.show()

