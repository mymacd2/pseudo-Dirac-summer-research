import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib



# Norms for the PMNS matrix squared:

u_e = [0, 0.674743, 0.302844, 0.0224125]
u_m = [0, 0.0946105, 0.360415, 0.544974]
u_t = [0, 0.230646, 0.33674,  0.432613]


# Effective distance traveled in 1oeV:
l_eff = 4.4434e23 * 5.06773093741 * 1e6

# Mass-difference squared (under assumption that all flavors have the same mass difference) in eV^2:
dm2_1 = 0
dm2_2 = 10**(-17.3)
dm2_3 = 10**(-17.3)

# Energy dist:
a = np.logspace(11, 13, num=2000)


# cos squared term:
osc1 = (np.cos((dm2_1 * l_eff)/(4*a)))**2
osc2 = (np.cos((dm2_2 * l_eff)/(4*a)))**2
osc3 = (np.cos((dm2_3 * l_eff)/(4*a)))**2


# Define start and end flavor states:
ui = u_m
uf1 = u_e
uf2 = u_t
uf3 = u_m


prob1 = osc1*(ui[1]*uf1[1]) + osc2*(ui[2]*uf1[2]) + osc3*(ui[3]*uf1[3])
prob2 = osc1*(ui[1]*uf2[1]) + osc2*(ui[2]*uf2[2]) + osc3*(ui[3]*uf2[3])
prob3 = osc1*(ui[1]*uf3[1]) + osc2*(ui[2]*uf3[2]) + osc3*(ui[3]*uf3[3])
prob4 = 1 - prob1 - prob2 - prob3

fig1, ax1 = plt.subplots()
µe, = ax1.plot(a, prob1)
µτ, = ax1.plot(a, prob2)
µµ, = ax1.plot(a, prob3)
µs, = ax1.plot(a, prob4)
ax1.set_xscale('log')
ax1.set_xticks([10**11, 10**12, 10**13])
ax1.set_xlabel(r"$E_{\nu}$ (eV)")
ax1.set_ylabel(r"$P(\nu_{e} \to \nu_{\alpha})$")
ax1.set_title(r"$\delta m_1^2 = 0$ eV$^2$, $\delta m_2^2 = 10^{-17.8}$ eV$^2$, $\delta m_3^2 = 10^{-17.3}$ eV$^2$", loc="left")
ax1.legend([µe, µτ, µµ, µs], [r"$\nu_{e} \to \nu_{e}$", r"$\nu_{e} \to \nu_{\tau}$",r"$\nu_{e} \to \nu_{\mu}$",r"$\nu_{e} \to \nu_{s}$"])
plt.show()





