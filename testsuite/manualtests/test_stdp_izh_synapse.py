# -*- coding: utf-8 -*-
#
# test_stdp_izh_synapse.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

import nest
import numpy as np

update_interval = 1000

simtime = 1000.0 *100 
delay = 10.0

w0 = 6.0

nest.SetKernelStatus({"syn_update_interval": update_interval, "resolution": 1.0})

np.random.seed(0)

pre_spike_times = np.unique(np.sort(np.random.randint(1, simtime, int(simtime/1000 * 20)))) 
#pre_spike_times = np.array([988., 989., 990.])
post_spike_times = np.unique(np.sort(np.random.randint(1, simtime, int(simtime/1000 * 1000)))) 
#post_spike_times = np.array([986.])

sg_pre = nest.Create("spike_generator", params={"spike_times": pre_spike_times.astype("double")})
sg_post = nest.Create("spike_generator", params={"spike_times": post_spike_times.astype("double")})

pre = nest.Create("parrot_neuron")
post = nest.Create("izhikevich", 1, params={"consistent_integration": False})

sd_pre = nest.Create("spike_detector")
sd_post = nest.Create("spike_detector")

nest.Connect(sg_pre, pre)
nest.Connect(sg_post, post, syn_spec={"weight": 20.0, "delay": delay})

nest.SetDefaults("stdp_izh_synapse", {"consistent_integration": False})
nest.Connect(pre, post, syn_spec={"model": "stdp_izh_synapse", "weight": w0, "delay": delay})

nest.Connect(pre, sd_pre)
nest.Connect(post, sd_post)

nest.Simulate(simtime)

pre_spikes = nest.GetStatus(sd_pre, "events")[0]["times"].astype("int")
post_spikes = nest.GetStatus(sd_post, "events")[0]["times"].astype("int")

w_end = nest.GetStatus(nest.GetConnections(pre, post), "weight")
print w_end

np.savetxt("times_zero.txt", pre_spikes, fmt="%d")
np.savetxt("times_one.txt", post_spikes, fmt="%d")

#n_updates = int( simtime / update_interval )
#
#ref_weight = w0
#
#sd = 0.0
#for n in range( n_updates ):
#    rel_pre_spikes = np.extract( np.all( [ pre_spikes > n * update_interval, pre_spikes <= ( n + 1 ) * update_interval ], axis = 0 ), pre_spikes )
#    rel_post_spikes = np.extract( np.all( [ post_spikes > n * update_interval, post_spikes <= ( n + 1 ) * update_interval ], axis = 0 ), post_spikes )
#    sd += np.sum( [ 0.1 * np.sum( np.exp( ( np.extract( pre_spikes <= post_spike, pre_spikes ) - post_spike ) / 20.0 ) ) for post_spike in rel_post_spikes ] )
#    sd -= np.sum( [ 0.1 * 1.2 * np.sum( np.exp( ( np.extract( post_spikes <= pre_spike, post_spikes ) - pre_spike ) / 20.0 ) ) for pre_spike in rel_pre_spikes ] )
#    ref_weight = w0 + sd
#    #ref_weight = ref_weight + sd + 0.01
#    #sd *= 0.9
#print ref_weight
