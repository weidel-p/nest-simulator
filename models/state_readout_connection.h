/*
 *  state_readout_connection.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef STATE_READOUT_H
#define STATE_READOUT_H

/* BeginDocumentation

   Name: stdp_dopamine_synapse - Synapse type for dopamine-modulated
                                 spike-timing dependent plasticity.

   Description:
   stdp_dopamine_synapse is a connection to create synapses with
   dopamine-modulated spike-timing dependent plasticity (used as a
   benchmark model in [1], based on [2]). The dopaminergic signal is a
   low-pass filtered version of the spike rate of a user-specific pool
   of neurons. The spikes emitted by the pool of dopamine neurons are
   delivered to the synapse via the assigned volume transmitter. The
   dopaminergic dynamics is calculated in the synapse itself.

   Examples:
   /volume_transmitter Create /vol Set
   /iaf_neuron Create /pre_neuron Set
   /iaf_neuron Create /post_neuron Set
   /iaf_neuron Create /neuromod_neuron Set
   /stdp_dopamine_synapse  << /vt vol >>  SetDefaults
   neuromod_neuron vol Connect
   pre_neuron post_neuron /stdp_dopamine_synapse Connect

   Parameters:
     Common properties:
           vt        long   - ID of volume_transmitter collecting the spikes
                              from the pool of dopamine releasing neurons and
                              transmitting the spikes to the synapse. A value of
                              -1 indicates that no volume transmitter has been
                              assigned.
           A_plus    double - Amplitude of weight change for facilitation
           A_minus   double - Amplitude of weight change for depression
           tau_plus  double - STDP time constant for facilitation in ms
           tau_c     double - Time constant of eligibility trace in ms
           tau_n     double - Time constant of dopaminergic trace in ms
           b         double - Dopaminergic baseline concentration
           Wmin      double - Minimal synaptic weight
           Wmax      double - Maximal synaptic weight

     Individual properties:
           c         double - eligibility trace
           n         double - neuromodulator concentration

   Remarks:
     The common properties can only be set by SetDefaults and apply to all
     synapses of the model.

   References:
   [1] Potjans W, Morrison A and Diesmann M (2010). Enabling
       functional neural circuit simulations with distributed
       computing of neuromodulated plasticity.
       Front. Comput. Neurosci. 4:141. doi:10.3389/fncom.2010.00141
   [2] Izhikevich, E.M. (2007). Solving the distal reward problem
       through linkage of STDP and dopamine signaling. Cereb. Cortex,
       17(10), 2443-2452.

   Transmits: SpikeEvent

   Author: Susanne Kunkel
   Remarks:
   - based on an earlier version by Wiebke Potjans
   - major changes to code after code revision in Apr 2013

   SeeAlso: volume_transmitter
*/

// Includes from libnestutil:
#include "numerics.h"

// Includes from models:
#include "volume_transmitter2.h"

// Includes from nestkernel:
#include "connection.h"
#include "spikecounter.h"

#include <iostream> 

namespace nest
{

/**
 * Class containing the common properties for all synapses of type dopamine
 * connection.
 */
class StateReadoutCommonProperties : public CommonSynapseProperties
{
public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  StateReadoutCommonProperties();

  /**
   * Get all properties and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  Node* get_node();

  long get_vt_gid() const;

  volume_transmitter2* vt_;
  double Aplus_;
  double Aminus_;
  double tau_plus_;
  double tau_n_;
  double b_;
  double n_upper_threshold_;
  double n_lower_threshold_;
  double mean_firing_rate_;
  double Wmin_;
  double Wmax_;
};

inline long
StateReadoutCommonProperties::get_vt_gid() const
{
  if ( vt_ != 0 )
    return vt_->get_gid();
  else
    return -1;
}

/**
 * Class representing an StateReadoutConnection with homogeneous parameters,
 * i.e. parameters are the same for all synapses.
 */
template < typename targetidentifierT >
class StateReadoutConnection : public Connection< targetidentifierT >
{

public:
  typedef StateReadoutCommonProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  StateReadoutConnection();

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  StateReadoutConnection( const StateReadoutConnection& );

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay;
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   */
  void send( Event& e, thread t, double_t, const StateReadoutCommonProperties& cp );

  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
    }
  };

  /*
   * This function calls check_connection on the sender and checks if the
   * receiver accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   * We have to override the base class' implementation, since for STDP
   * connections we have to call register_stdp_pl_connection on the target
   * neuron to inform the Archiver to collect spikes for this connection.
   * Further, the STDP dopamine synapse requires a volume transmitter to be set
   * before any simulation is performed. Checking this satisfies ticket #926.
   *
   * \param s The source node
   * \param r The target node
   * \param receptor_type The ID of the requested receptor type
   * \param t_lastspike last spike produced by presynaptic neuron (in ms)
   */
  void
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    double t_lastspike,
    const CommonPropertiesType& cp )
  {
    if ( cp.vt_ == 0 )
      throw BadProperty(
        "No volume transmitter has been assigned to the synapse." );

    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike - get_delay() );
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:

  void update_weight_( double n, double t, const StateReadoutCommonProperties& cp );

  // data members of each connection
  double weight_;
  double Kplus_;
  double Kminus_;

  // time of last update, which is either time of last presyn. spike or
  // time-driven update
  double t_last_update_;
};

//
// Implementation of class StateReadoutConnection.
//

template < typename targetidentifierT >
StateReadoutConnection< targetidentifierT >::StateReadoutConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , Kplus_( 0.0 )
  , Kminus_( 0.0 )
  , t_last_update_( 0.0 )
{
}

template < typename targetidentifierT >
StateReadoutConnection< targetidentifierT >::StateReadoutConnection(
  const StateReadoutConnection& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , Kplus_( rhs.Kplus_ )
  , Kminus_( rhs.Kminus_)
  , t_last_update_( rhs.t_last_update_ )
{
}

template < typename targetidentifierT >
void
StateReadoutConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

  // base class properties, different for individual synapse
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );

  // own properties, different for individual synapse
  def< double >( d, "Kplus", Kplus_ );
  def< double >( d, "Kminus", Kminus_);
}

template < typename targetidentifierT >
void
StateReadoutConnection< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  // base class properties
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );

  updateValue< double >( d, "Kplus", Kplus_ );
  updateValue< double >( d, "Kminus", Kminus_);
}


template < typename targetidentifierT >
inline void
StateReadoutConnection< targetidentifierT >::update_weight_( double n, double t, const StateReadoutCommonProperties& cp )
{

    Kplus_ = Kplus_ * std::exp( ( t_last_update_ - t) / cp.tau_plus_ );
    t_last_update_ = t;

    double norm_w = (weight_ - cp.Wmin_) / (cp.Wmax_ - cp.Wmin_);
    double dw = Kplus_ * Kminus_;
    if (n > cp.n_upper_threshold_ ){
       //facilitate
       norm_w += cp.Aplus_ * (1 - norm_w) * (dw / (1 + std::abs(dw)));
    }
    else if (n < cp.n_lower_threshold_){
       //depress
       norm_w -= cp.Aminus_ * norm_w * (dw / (1 + std::abs(dw)));
    }
    weight_ = norm_w * (cp.Wmax_ - cp.Wmin_) + cp.Wmin_;

    if ( weight_ < cp.Wmin_ )
      weight_ = cp.Wmin_;
    if ( weight_ > cp.Wmax_ )
      weight_ = cp.Wmax_;

}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
template < typename targetidentifierT >
inline void
StateReadoutConnection< targetidentifierT >::send( Event& e,
  thread t,
  double_t,
  const StateReadoutCommonProperties& cp )
{

  Node* target = get_target( t );

  // purely dendritic delay
  double dendritic_delay = get_delay();

  double t_spike = e.get_stamp().get_ms();



  //// get dopa history in relevant range (t_last_update, t_spike] from
  //std::map< double, double >::iterator start;
  //std::map< double, double >::iterator finish;
  //cp.vt_->get_history( t_last_update_,
  //  t_spike,
  //  &start,
  //  &finish );
  //
  //for(; start != finish; start++) {
  //  double t = start->first;
  //  double n = start->second;
  //  Kminus_ = target->get_K_value( t );
  //
  //  //TODO when to update weight? only for presynaptic spikes? or for pre and post? or always?
  //  update_weight_( n, t, cp );
  //  
  //
  //}

  double n = cp.vt_->get_n( t_spike );
  Kminus_ = target->get_K_value( t_spike );
  Kplus_ = Kplus_ * std::exp( ( t_last_update_ - t_spike) / cp.tau_plus_ ) + 1.0;
  
  update_weight_( n, t_spike, cp );

  e.set_receiver( *target );
  e.set_weight( weight_ );
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e.set_Kplus( Kplus_ );
  e.set_Kminus( Kminus_ );
  e.set_dopa( n );
  e();

  t_last_update_ = t_spike;
}

} // of namespace nest

#endif // of #ifndef STATE_READOUT_H
