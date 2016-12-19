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

/* 
 * Needed quantities:
 * Presynaptic activity: Kplus_
 * Postsynaptic activity: Kminus_ (computed in target neuron)
 * dopamine concentration: n_
 * synaptic weight: weight_
 * last update: t_last_update_
 *
 * Assume:
 * Kplus_ and n_ and weight_ are already updated to t_last_update_ 
 * There was no presynaptic spike between t_last_update_ and the current update step
 * 
 * 
 * Pseudocode
 *
 * function send():
 *     t_now = time of current presynaptic spike
 *     dopa_spikes = get all dopamine spike from the volume transmitter
 *     process_dopa_spikes()
 *     e()
 *
 *
 * function trigger_update_weight(dopa_spikes):
 *     limit dopa spikes to first dopa spike later than t_last_update_
 *     process_dopa_spikes() 
 * 
 *
 * function process_dopa_spikes():
 *     for t_dopa in [dopa_spikes, t_now]:
 *         post_spikes = target->getHistory(t_last_update_, t_dopa)
 *         for t_s in [post_spikes, t_dopa]:
 *             Kminus_ = get Kminus_ at time t_last_update_
 *
 *             needed_t0 = check_if_update_weight_needed() 
 *             propagate Kplus_, Kminus_, n_ to t_s in temporary variables 
 *             needed_t1 = check_if_update_weight_needed() 
 *
 *             if needed_t0 == "facilitate":
 *                  t_update_until = calc_update_needed_until(t_s)
 *                  facilitate(t_last_update_, t_update_until)
 *             if needed_t0 == "no update needed" and needed_t1 == "facilitate":
 *                  t_update_from = calc_update_needed_from(t_s)
 *                  facilitate(t_update_from, t_s)
 *             
 *             if needed_t0 == "depress":
 *                  t_update_until = calc_update_needed_until(t_s)
 *                  depress(t_last_update_, t_update_until)
 *             if needed_t0 == "no update needed" and needed_t1 == "depress":
 *                  t_update_from = calc_update_needed_from(t_s)
 *                  depress(t_update_from, t_s)
 *             if needed_t0 == "facilitate" and needed_t1 == "depress":
 *                  t_update_until = calc_update_needed_until(t_s)
 *                  facilitate(t_last_update_, t_update_until)
 *                  t_update_from = calc_update_needed_from(t_s)
 *                  depress(t_update_from, t_s)
 *            
 *
 *             
 *             t_last_update_ = t_s
 *             update Kplus_, n_ according to temporary values
 *             n_ += 1
 *             Kplus += 1
 *
 *
 *
 * 
 * function check_if_update_weight_needed():
 *         if n_ > n_upper_threshold:
 *             if Kminus_ > Kminus_threshold and Kplus > Kplus_threshold:
 *                 return "facilitate"
 *             else:
 *                 return "no update needed"
 *         if n < n_lower_threshold:
 *             if Kminus_ > Kminus_threshold and Kplus > Kplus_threshold:
 *                 return "depress"
 *             else:
 *                 return "no update needed"
 *
 *
 * function facilitate(t0, t1):
 *     dt = t1 - t0
 *     // additive rule
 *     weight_ += A * Kminus_ * Kplus_ / -(1/tau_minus + 1/tau_plus) * (np.exp(-dt/tau_minus) * np.exp(-dt/tau_plus) - 1)
 *
 *
 * function depress(t0, t1):
 *     dt = t1 - t0
 *     // additive rule
 *     weight_ -= A *  Kminus_ * Kplus_ / -(1/tau_minus + 1/tau_plus) * (np.exp(-dt/tau_minus) * np.exp(-dt/tau_plus) - 1)
 *
 *
 *
 * function calc_update_needed_from(t1)
 *     t_threshold_cross_n = -tau_n * ln(n_lower_threshold / n_)
 *     return min(t_threshold_cross_n, t1)
 *
 *
 *
 * function calc_update_needed_until(t1)
 *     t_threshold_cross_K_minus = -tau_minus * ln(Kminus_threshold / Kminus)
 *     t_threshold_cross_K_plus = -tau_plus * ln(Kplus_threshold / Kplus)
 *     if n_ > n_upper_threshold:
 *         t_threshold_cross_n = -tau_n * ln(n_upper_threshold / n_)
 *     else:
 *         t_threshold_cross_n = t1 
 *
 *     return min(t_threshold_cross_K_minus, t_threshold_cross_K_plus, t_threshold_cross_n, t1)
 *     
 *         
 *
 *
 */




// Includes from libnestutil:
#include "numerics.h"

// Includes from models:
#include "volume_transmitter.h"

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

  volume_transmitter* vt_;
  double Aplus_;
  double Aminus_;
  double tau_plus_;
  double tau_n_;
  double b_;
  double Kplus_threshold_;
  double Kminus_threshold_;
  double n_upper_threshold_;
  double n_lower_threshold_;
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

  void trigger_update_weight( thread t,
    const std::vector< spikecounter >& dopa_spikes,
    double t_trig,
    const StateReadoutCommonProperties& cp );

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
  void process_dopa_spikes_(
    thread t,
    double t1,
    const std::vector< spikecounter >& dopa_spikes,
    const StateReadoutCommonProperties& cp);

  void facilitate_(
    thread t,
    double t0,
    double t1,
    const StateReadoutCommonProperties& cp);

  void depress_(
    thread t,
    double t0,
    double t1,
    const StateReadoutCommonProperties& cp);

  void process_next_(
    thread t,
    double t0,
    double t1,
    const StateReadoutCommonProperties& cp );

  double calc_update_needed_from(
    double t0,
    double t1,
    const StateReadoutCommonProperties& cp );

  double calc_update_needed_until( thread t,
    double t0,
    double t1,
    const StateReadoutCommonProperties& cp );

  int check_if_update_needed(
    double Kplus,
    double Kminus,
    double n,
    const StateReadoutCommonProperties& cp );




  // data members of each connection
  double weight_;
  double Kplus_;
  double Kminus_;
  double n_;


  // time of last update, which is either time of last presyn. spike or
  // time-driven update
  double t_last_update_;
  int dopa_spikes_idx_;
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
  , n_( 0.0 )
  , dopa_spikes_idx_( 0 )
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
  , n_( rhs.n_ )
  , dopa_spikes_idx_( rhs.dopa_spikes_idx_ )
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
  def< double >( d, "n", n_ );
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

  updateValue< double >( d, "n", n_ );
  updateValue< double >( d, "Kplus", Kplus_ );
  updateValue< double >( d, "Kminus", Kminus_);
}




template < typename targetidentifierT >
inline void
StateReadoutConnection< targetidentifierT >::facilitate_(
  thread t,
  double t0,
  double t1,
  const StateReadoutCommonProperties& cp){

  assert(t1 >= t0);
  double dt = t1 - t0;
  double tau_minus = get_target( t )->get_tau_minus();
  
  weight_ = cp.Wmax_ - (cp.Wmax_ - weight_)
            * std::exp( -1/cp.Wmax_ * cp.Aplus_ * Kminus_ * Kplus_ / -(1/tau_minus + 1/cp.tau_plus_)
            * (std::exp(-dt/tau_minus) * std::exp(-dt/cp.tau_plus_) - 1));

  if ( weight_ > cp.Wmax_ ){
    weight_ = cp.Wmax_;
  }
  //std::cout << "POTENTIATE " << weight_ << " n " << n_ << " kp " << Kplus_ << " km " << Kminus_ << " t0 " << t0 << " t1 " << t1 << std::endl;

}

template < typename targetidentifierT >
inline void
StateReadoutConnection< targetidentifierT >::depress_(
  thread t,
  double t0,
  double t1,
  const StateReadoutCommonProperties& cp){
  
  assert(t1 >= t0);
  double dt = t1 - t0;
  double tau_minus = get_target( t )->get_tau_minus();

  weight_ = cp.Wmin_ + (weight_ - cp.Wmin_)
            * std::exp( -1/cp.Wmin_ * cp.Aminus_ * Kminus_ * Kplus_ / -(1/tau_minus + 1/cp.tau_plus_)
            * (std::exp(-dt/tau_minus) * std::exp(-dt/cp.tau_plus_) - 1));
  
  if ( weight_ < cp.Wmin_ )
  {
    weight_ = cp.Wmin_;
  }
  //std::cout << "DEPRESS " << weight_ << " n " << n_ << " kp " << Kplus_ << " km " << Kminus_ << " t0 " << t0 << " t1 " << t1 << std::endl;
}



template < typename targetidentifierT >
inline void
StateReadoutConnection< targetidentifierT >::process_next_(
  thread t,
  double t0,
  double t1,
  const StateReadoutCommonProperties& cp )
{
  //std::cout << t0 << " " << t1 << std::endl;
 
  assert(t1 >= t0);
  Node* target = get_target( t );

  // update all traces to  t0
  Kminus_ = target->get_K_value( t0 );

  int update_needed_t0 = check_if_update_needed( Kplus_, Kminus_, n_, cp ); 

  // propagate Kplus_, Kminus_, n_ to t1 in temporary variables 
  double Kminus_t1 = Kminus_ * std::exp( ( t0 - t1 ) / target->get_tau_minus() );
  double Kplus_t1 = Kplus_ * std::exp( ( t0 - t1 ) / cp.tau_plus_ );
  double n_t1 = n_ * std::exp( ( t0 - t1) / cp.tau_n_ );

  int update_needed_t1 = check_if_update_needed( Kplus_t1, Kminus_t1, n_t1, cp ); 

  // if potentitaion at t0
  if (update_needed_t0 == 1 and update_needed_t1 != -1)
  {
      //std::cout << "FACI NEEDED until t0: " << t0 << " " << t1 << " " << n_ << " " << Kplus_ << " " << Kminus_ << std::endl;
      //std::cout << "FACI NEEDED until t1:" << n_t1 << " " << Kplus_t1 << " " << Kminus_t1 << std::endl;
       double t_until = calc_update_needed_until(t, t0, t1, cp);
       assert(t_until >= t0);
       facilitate_(t, t0, t_until, cp );
  }

  // if potentitaion at t1
  else if (update_needed_t0 == 0 and update_needed_t1 == 1)
  {
      //TODO this case should never happen! 
      std::cout << "FACI NEEDED from t0: " << t0 << " " << t1 << " " << n_ << " " << Kplus_ << " " << Kminus_ << std::endl;
      std::cout << "FACI NEEDED from t1:" << n_t1 << " " << Kplus_t1 << " " << Kminus_t1 << std::endl;
       double t_from = calc_update_needed_from(t0, t1, cp);
       assert(t1 >= t_from);
       facilitate_(t, t_from, t1, cp);
  }

  // if depression at t0
  else if (update_needed_t0 == -1)
  {
      //std::cout << "DEPRE NEEDED" << std::endl;
       double t_until = calc_update_needed_until(t, t0, t1, cp);
       assert(t_until >= t0);
       depress_(t, t0, t_until, cp);
  }

  // if depression at t1
  else if (update_needed_t0 == 0 and update_needed_t1 == -1)
  {
      //std::cout << "DEPRE NEEDED from" << std::endl;
       double t_from = calc_update_needed_from(t0, t1, cp);
       assert(t1 >= t_from);
       depress_(t, t_from, t1, cp);
  }

  // if potentiation at t0 and depression at t1
  else if (update_needed_t0 == 1 and update_needed_t1 == -1)
  {
      std::cout << "FACI AND DEPRE NEEDED" << std::endl;
       double t_until = calc_update_needed_until(t, t0, t1, cp);
       assert(t_until >= t0);
       facilitate_(t, t0, t_until, cp);
       double t_from = calc_update_needed_from(t0, t1, cp);
       assert(t1 >= t_from);
       depress_(t, t_from, t1, cp);
  }

  // update Kplus_, Kminus_, n_ according to temporary values
  Kplus_ = Kplus_t1;
  Kminus_ = Kminus_t1;
  n_ = n_t1;


}


template < typename targetidentifierT >
inline int 
StateReadoutConnection< targetidentifierT >::check_if_update_needed(
  double Kplus,
  double Kminus,
  double n,
  const StateReadoutCommonProperties& cp )
{

  if (Kminus > cp.Kminus_threshold_ and Kplus > cp.Kplus_threshold_)
  {
    if (n > cp.n_upper_threshold_)
    {
      //std::cout << "update needed facili " <<  " n " << n << " kp " << Kplus << " km " << Kminus  << std::endl;
      return 1; 
    }
    else if (n < cp.n_lower_threshold_)
    {
      //std::cout << "update needed depp " <<  " n " << n << " kp " << Kplus << " km " << Kminus  << std::endl;
      return -1; 
    }
  }
  return 0;
}

template < typename targetidentifierT >
inline double 
StateReadoutConnection< targetidentifierT >::calc_update_needed_from(
  double t0,
  double t1,
  const StateReadoutCommonProperties& cp )
{
  assert(t1 >= t0);
  double t_threshold_cross_n = -cp.tau_n_ * std::log(cp.n_lower_threshold_ / n_);
 // std::cout << t0 << " from " << std::min(t_threshold_cross_n, t1)<< std::endl;
  return t0 + std::min(t_threshold_cross_n, t1);
}


template < typename targetidentifierT >
inline double 
StateReadoutConnection< targetidentifierT >::calc_update_needed_until( thread t,
  double t0,
  double t1,
  const StateReadoutCommonProperties& cp )
{
  assert(t1 >= t0);

  Node* target = get_target( t );
  
  double t_threshold_cross_K_minus = -target->get_tau_minus() * std::log(cp.Kminus_threshold_ / Kminus_);
  if (t_threshold_cross_K_minus < 0){
	std::cout <<  " km " <<  Kminus_ << " " << cp.Kminus_threshold_ << std::endl;
  }
  double t_threshold_cross_K_plus = -cp.tau_plus_ * std::log(cp.Kplus_threshold_ / Kplus_);
  if (t_threshold_cross_K_plus < 0){
	std::cout << " kp " <<  Kplus_ << " " << cp.Kplus_threshold_ << std::endl;
  }


  double t_threshold_cross_n = t1; 
  if (n_ > cp.n_upper_threshold_)
      t_threshold_cross_n = -cp.tau_n_ * std::log(cp.n_upper_threshold_ / n_);
  
  if (t_threshold_cross_n < 0){
	std::cout << " n " <<  n_ << " " << cp.n_upper_threshold_ << std::endl;
  }

  //std::cout << t0 << " until " << std::min(t_threshold_cross_K_minus, std::min(t_threshold_cross_K_plus, std::min(t_threshold_cross_n, t1)))<< std::endl;
  return t0 + std::min(t_threshold_cross_K_minus, std::min(t_threshold_cross_K_plus, std::min(t_threshold_cross_n, t1)));
  
}



template < typename targetidentifierT >
inline void
StateReadoutConnection< targetidentifierT >::process_dopa_spikes_(
  thread t,
  double t1,
  const std::vector< spikecounter >& dopa_spikes,
  const StateReadoutCommonProperties& cp )
{

  Node* target = get_target( t );


  std::deque< histentry >::iterator start;
  std::deque< histentry >::iterator finish;

  for(std::vector< spikecounter >::const_iterator it = dopa_spikes.begin(); it != dopa_spikes.end(); ++it) {
    if (it->spike_time_ <= t_last_update_ or it->spike_time_ > t1)
        continue;
  
    // get history of postsynaptic neuron
    target->get_history( t_last_update_,
    it->spike_time_,
    &start,
    &finish );
    while ( start != finish )
    {
        //std::cout << "process post spike " << n_ << std::endl;
       process_next_(t, t_last_update_, start->t_, cp);

       t_last_update_ = start->t_;
       ++start;

    }
    
    process_next_(t, t_last_update_, it->spike_time_, cp);
    n_ += it->multiplicity_ / cp.tau_n_; 
    // std::cout << "process dopa spike " << n_ << std::endl;
    t_last_update_ = it->spike_time_;
  }



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
  // t_lastspike_ = 0 initially

  double t_spike = e.get_stamp().get_ms();


  // get history of dopamine spikes
  const std::vector< spikecounter >& dopa_spikes = cp.vt_->deliver_spikes();

  process_dopa_spikes_(t, t_spike,  dopa_spikes, cp);

  if (t_spike > t_last_update_){
  //    std::cout << "ERRER " << t_last_update_ << " " << t_spike << std::endl;
      process_next_(t, t_last_update_, t_spike, cp);
      t_last_update_ = t_spike;
  }


  Kminus_ = target->get_K_value( t_spike);

  e.set_receiver( *target );
  e.set_weight( weight_ );
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e.set_Kplus( Kplus_ );
  e.set_Kminus( Kminus_ );
  e.set_dopa( n_ );
  e();

  Kplus_ += 1.;
  t_last_update_ = t_spike;
}

template < typename targetidentifierT >
inline void
StateReadoutConnection< targetidentifierT >::trigger_update_weight( thread t,
  const std::vector< spikecounter >& dopa_spikes,
  const double t_trig,
  const StateReadoutCommonProperties& cp )
{

  process_dopa_spikes_(t, t_trig, dopa_spikes, cp);
  process_next_(t, t_last_update_, t_trig, cp);
  t_last_update_ = t_trig;


}





} // of namespace nest

#endif // of #ifndef STATE_READOUT_H
