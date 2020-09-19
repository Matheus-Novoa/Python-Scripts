from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def Langmuir(P, Qm, B):
    return Qm * B * P / (1 + B * P)


press = [5, 10, 25, 50, 100, 200]
conc = [2.35e-04, 3.92e-04, 8.63e-04, 1.65e-03, 2.82e-03, 4.74e-03]

params, cov = curve_fit(Langmuir, press, conc)

print(f'\nQm: {params[0]}\nB: {params[1]}')

conc_test = [Langmuir(p, params[0], params[1]) for p in press]

plt.plot(press, conc, marker='o')
plt.plot(press, conc_test, marker='o')
plt.show()