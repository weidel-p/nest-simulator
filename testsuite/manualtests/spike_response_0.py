import numpy as np
import pylab as plt
import nest
from nest import voltage_trace as vtr
from nest import raster_plot as rplt

eps0 = 20.
tau_mem = 20.
tau_syn = 5.

eps = lambda t: eps0 / (tau_mem - tau_syn) * (np.exp(-t/tau_mem) - np.exp(-t/tau_syn))
deps = lambda t: eps0 / (tau_mem - tau_syn) * ( 1/tau_syn * np.exp(-t/tau_syn) - 1/tau_mem * np.exp(-t/tau_mem))
ddeps = lambda t: eps0 / (tau_mem - tau_syn) * ( -1/(tau_syn * tau_syn) * np.exp(-t/tau_syn) + 1/(tau_mem * tau_mem) * np.exp(-t/tau_mem))


test = lambda t: (np.e/2)/(tau_syn) * np.exp(-t/(tau_syn)) - (np.e/2)/tau_mem * np.exp(-t/tau_mem)

I = [0.]
V_m = 0.

def update(dt):
    I.append(deps(dt))


tmp = [0]
for t in np.linspace(0, 20, 200):
    tmp.append(tmp[-1] + deps(t) * 0.1)


#plt.plot(eps(np.linspace(0, 20, 20)) )
#plt.plot(deps(np.linspace(0, 20, 20)) )
#plt.plot(ddeps(np.linspace(0, 20, 200)) )
#plt.plot(test(np.linspace(0, 20, 20)) )
#plt.plot(tmp)

plt.figure()


nest.SetKernelStatus({"resolution": 1.0, "print_time": True})

n = nest.Create("SRM0", params={"tau_syn": tau_syn, "tau_m": tau_mem, "xi": -5.0, "theta": 16., "du": 2., "rho0": 60.0})
mm = nest.Create("multimeter", params={"record_from": ["V_m"]})
pg = nest.Create("poisson_generator", params={"rate": 2000.})
sd = nest.Create("spike_detector")

nest.Connect(pg, n, syn_spec={"weight": 1.5})
nest.Connect(pg, n, syn_spec={"weight": -1.})
nest.Connect(mm, n)
nest.Connect(n, sd)

nest.Simulate(2000.)

vtr.from_device(mm)
#rplt.from_device(sd)

vtr.show()



#plt.show()
