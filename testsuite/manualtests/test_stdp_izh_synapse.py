import nest
import numpy as np

update_interval = 1000

simtime = 1000.0
delay = 20.0

w0 = 6.0

nest.SetKernelStatus({"syn_update_interval": update_interval})

sg_pre = nest.Create("poisson_generator", params={"rate": 20.0})
sg_post = nest.Create("poisson_generator", params={"rate": 1000.0})

pre = nest.Create("parrot_neuron")
post = nest.Create("izhikevich")

sd_pre = nest.Create("spike_detector")
sd_post = nest.Create("spike_detector")

nest.Connect(sg_pre, pre)
nest.Connect(sg_post, post, syn_spec={"weight": 10.0, "delay": delay})

nest.Connect(pre, post, syn_spec={"model": "stdp_izh_synapse", "weight": w0, "delay": delay})

nest.Connect(pre, sd_pre)
nest.Connect(post, sd_post)

nest.Simulate(simtime)

pre_spikes = np.insert( nest.GetStatus(sd_pre, "events")[0]["times"] + delay, 0, float("-inf") )
post_spikes = np.insert( nest.GetStatus(sd_post, "events")[0]["times"], 0, float("-inf") )

n_updates = int( simtime / update_interval )

ref_weight = w0

sd = 0.0
for n in range( n_updates ):
    rel_pre_spikes = np.extract( np.all( [ pre_spikes > n * update_interval, pre_spikes <= ( n + 1 ) * update_interval ], axis = 0 ), pre_spikes )
    rel_post_spikes = np.extract( np.all( [ post_spikes > n * update_interval, post_spikes <= ( n + 1 ) * update_interval ], axis = 0 ), post_spikes )
    sd += np.sum( [ 0.1 * np.sum( np.exp( ( np.extract( pre_spikes <= post_spike, pre_spikes ) - post_spike ) / 20.0 ) ) for post_spike in rel_post_spikes ] )
    sd -= np.sum( [ 0.1 * 1.2 * np.sum( np.exp( ( np.extract( post_spikes <= pre_spike, post_spikes ) - pre_spike ) / 20.0 ) ) for pre_spike in rel_pre_spikes ] )
    ref_weight = w0 + sd
    #ref_weight = ref_weight + sd + 0.01
    #sd *= 0.9
    print ref_weight
