from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def Langmuir(P, Qm, B):
    return Qm * B * P / (1 + B * P)


press = [5, 10, 25, 50, 100, 200]
conc = [1.73e-04, 2.51e-4, 4.7e-4, 1.16e-3, 2.01e-3, 3.36e-3]

params, cov = curve_fit(Langmuir, press, conc)

print(f'\nQm: {params[0]}\nB: {params[1]}')

conc_test = [Langmuir(p, params[0], params[1]) for p in press]

plt.plot(press, conc, marker='o')
plt.plot(press, conc_test, '--')
plt.show()