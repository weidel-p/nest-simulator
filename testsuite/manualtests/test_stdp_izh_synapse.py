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

simtime = 5000.0
delay = 20.0

w0 = 6.0

nest.SetKernelStatus({"syn_update_interval": update_interval, "resolution": 1.0})
nest.set_verbosity("M_FATAL")

sg_pre = nest.Create("poisson_generator", params={"rate": 2.0})
sg_post = nest.Create("poisson_generator", params={"rate": 100.0})

pre = nest.Create("parrot_neuron")
post = nest.Create("izhikevich")

sd_pre = nest.Create("spike_detector")
sd_post = nest.Create("spike_detector")

nest.Connect(sg_pre, pre)
nest.Connect(sg_post, post, syn_spec={"weight": 10.0, "delay": delay})

nest.Connect(pre, post, syn_spec={"model": "stdp_izh_synapse", "weight": w0, "delay": delay, "alpha": 0., "lambda": 0.2})

nest.Connect(pre, sd_pre)
nest.Connect(post, sd_post)

n_updates = int( simtime / update_interval )
for _ in range(n_updates):
    nest.Simulate(1000.)
    print "nest", nest.GetStatus(nest.GetConnections(pre, post), "weight")[0]

pre_spikes = np.insert( nest.GetStatus(sd_pre, "events")[0]["times"] + delay, 0, float("-inf") )
post_spikes = np.insert( nest.GetStatus(sd_post, "events")[0]["times"], 0, float("-inf") )


ref_weight = w0

sd = 0.0
for n in range( n_updates ):
    rel_pre_spikes = np.extract( np.all( [ pre_spikes > n * update_interval, pre_spikes <= ( n + 1 ) * update_interval ], axis = 0 ), pre_spikes )
    rel_post_spikes = np.extract( np.all( [ post_spikes > n * update_interval, post_spikes <= ( n + 1 ) * update_interval ], axis = 0 ), post_spikes )
    sd += np.sum( [ 0.2 * np.sum( np.exp( ( np.extract( pre_spikes <= post_spike, pre_spikes ) - post_spike ) / 20.0 ) ) for post_spike in rel_post_spikes ] )
    sd -= np.sum( [ 0.2 * 0 *1.2 * np.sum( np.exp( ( np.extract( post_spikes <= pre_spike, post_spikes ) - pre_spike ) / 20.0 ) ) for pre_spike in rel_pre_spikes ] )
    ref_weight = w0 + sd
    #ref_weight = ref_weight + sd + 0.01
    #sd *= 0.9
    if ref_weight < 0:
        ref_weight = 0.
    if ref_weight > 10.:
        ref_weight = 10.
    print ref_weight


t_pre = nest.GetStatus(sd_pre, "events")[0]["times"]
t_post = nest.GetStatus(sd_post, "events")[0]["times"]

t_syn = []
for t_p in t_pre:
    for t_po in t_post:
        if (t_p + delay ) == t_po:
            t_syn.append(t_po)

print t_po
print t_pre, t_post


#print t_pre[np.where(np.logical_and(t_pre > 4000, t_pre < 5000))], t_post[np.where(np.logical_and(t_post > 4000, t_post < 5000))]


