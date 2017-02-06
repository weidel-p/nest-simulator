import nest
from nest import raster_plot as rplt
import pylab as plt

nest.ResetKernel()
nest.SetKernelStatus({"overwrite_files": True, "print_time": False, "local_num_threads": 1})
nest.set_verbosity("M_FATAL")

vt = nest.Create("volume_transmitter")
n = nest.Create("parrot_neuron", 2, {"tau_minus": 1000.})
pg = nest.Create("poisson_generator", 1, {"rate": 2.})
pg_dopa = nest.Create("poisson_generator", 1, {"rate": 10000.})
parrot = nest.Create("parrot_neuron", 1)

nest.Connect(pg, n)
nest.Connect(pg_dopa, parrot)
nest.Connect(parrot, vt)

nest.SetDefaults("state_readout_synapse", {"vt": vt[0], "A": 0.00001, "Kminus_threshold": 2.5, "tau_plus": 1000., "n_threshold": 10., "n": 10., "Wmax": 5., "Kplus": 2.5, "Kminus": 2.5, "weight": 2.5})
nest.Connect([n[0]], [n[1]], syn_spec={"model": "state_readout_synapse"})

sd = nest.Create("spike_detector")

nest.Connect(n, sd)

print nest.GetStatus(nest.GetConnections([n[0]], [n[1]]))

dopa = []
w = []
kplus = []
kminus = []

c = nest.GetConnections([n[0]], [n[1]])

simtime = 10000

for i in range(simtime):
    nest.Simulate(1)
    dopa.append(nest.GetStatus(c, "n"))
    w.append(nest.GetStatus(c, "weight"))
    kplus.append(nest.GetStatus(c, "Kplus"))
    kminus.append(nest.GetStatus(c, "Kminus"))
nest.SetStatus(pg_dopa, {"rate": 1000.})
for i in range(simtime):
    nest.Simulate(1)
    dopa.append(nest.GetStatus(c, "n"))
    w.append(nest.GetStatus(c, "weight"))
    kplus.append(nest.GetStatus(c, "Kplus"))
    kminus.append(nest.GetStatus(c, "Kminus"))
nest.SetStatus(pg_dopa, {"rate": 10000.})
for i in range(simtime):
    nest.Simulate(1)
    dopa.append(nest.GetStatus(c, "n"))
    w.append(nest.GetStatus(c, "weight"))
    kplus.append(nest.GetStatus(c, "Kplus"))
    kminus.append(nest.GetStatus(c, "Kminus"))
nest.SetStatus(pg_dopa, {"rate": 20000.})
for i in range(simtime):
    nest.Simulate(1)
    dopa.append(nest.GetStatus(c, "n"))
    w.append(nest.GetStatus(c, "weight"))
    kplus.append(nest.GetStatus(c, "Kplus"))
    kminus.append(nest.GetStatus(c, "Kminus"))


print nest.GetStatus(c)
print nest.GetStatus(sd)

plt.figure()
plt.plot(dopa)
plt.title("dopa")
plt.figure()
plt.plot(w)
plt.title("w")
plt.figure()
plt.plot(kplus)
plt.title("Kplus")
plt.figure()
plt.plot(kminus)
plt.title("Kminus")

#rplt.from_device(sd, hist=True)
plt.show()








