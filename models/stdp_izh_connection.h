/*
 *  stdp_izh_connection.h
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

#ifndef STDP_IZH_CONNECTION_H
#define STDP_IZH_CONNECTION_H

/* BeginDocumentation
  Name: stdp_izh_synapse - Synapse type for spike-timing dependent
   plasticity.

  Description:
   stdp_izh_synapse is a connector to create synapses with spike time
   dependent plasticity (as defined in [1]). Here the weight dependence
   exponent can be set separately for potentiation and depression.

  Examples:
   multiplicative STDP [2]  mu_plus = mu_minus = 1.0
   additive STDP       [3]  mu_plus = mu_minus = 0.0
   Guetig STDP         [1]  mu_plus = mu_minus = [0.0,1.0]
   van Rossum STDP     [4]  mu_plus = 0.0 mu_minus = 1.0

  Parameters:
   tau_plus   double - Time constant of STDP window, potentiation in ms
                       (tau_minus defined in post-synaptic neuron)
   lambda     double - Step size
   alpha      double - Asymmetry parameter (scales depressing increments as
                       alpha*lambda)
   mu_plus    double - Weight dependence exponent, potentiation
   mu_minus   double - Weight dependence exponent, depression
   Wmax       double - Maximum allowed weight

  Transmits: SpikeEvent

  References:
   [1] Guetig et al. (2003) Learning Input Correlations through Nonlinear
       Temporally Asymmetric Hebbian Plasticity. Journal of Neuroscience

   [2] Rubin, J., Lee, D. and Sompolinsky, H. (2001). Equilibrium
       properties of temporally asymmetric Hebbian plasticity, PRL
       86,364-367

   [3] Song, S., Miller, K. D. and Abbott, L. F. (2000). Competitive
       Hebbian learning through spike-timing-dependent synaptic
       plasticity,Nature Neuroscience 3:9,919--926

   [4] van Rossum, M. C. W., Bi, G-Q and Turrigiano, G. G. (2000).
       Stable Hebbian learning from spike timing-dependent
       plasticity, Journal of Neuroscience, 20:23,8812--8821

  FirstVersion: March 2006
  Author: Moritz Helias, Abigail Morrison
  Adapted by: Philipp Weidel
  SeeAlso: synapsedict, tsodyks_synapse, static_synapse
*/

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"


namespace nest
{

// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template
template < typename targetidentifierT >
class STDPIzhConnection : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPIzhConnection();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPIzhConnection( const STDPIzhConnection& );

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
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
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses (empty).
   */
  void send( Event& e,
    thread t,
    double t_lastspike,
    const CommonSynapseProperties& cp );


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

  void
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    double t_lastspike,
    const CommonPropertiesType& )
  {
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

  // data members of each connection
  double weight_;
  double tau_plus_;
  double tau_minus_;
  double tau_deriv_;
  double lambda_;
  double alpha_;
  double wdev_;
  double dw0_;
  double Wmax_;
  double t_lastpostspike_;
  double next_update_;
  double update_interval_;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param t_lastspike Time point of last spike emitted
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline void
STDPIzhConnection< targetidentifierT >::send( Event& e,
  thread t,
  double t_lastspike,
  const CommonSynapseProperties& )
{
  // synapse STDP depressing/facilitation dynamics

  double_t t_spike = e.get_stamp().get_ms();
    
  // t_lastspike_ = 0 initially
  // assume axonal delay = dendritic delay
  double dendritic_delay = get_delay() / 2.;
  double_t axonal_delay = dendritic_delay;

  //update actual weight if we have gone past an update interval
  if (t_spike >= next_update_)
  {
    //if ( e.get_sender().get_gid() == 10 )
    //std::cout << "old weight: " << weight_ << std::endl;
    double_t w_new = weight_ + wdev_ + dw0_;
    if ( w_new > 0.0 )
    {
      if ( w_new < Wmax_ )
	weight_ = w_new;
      else
	weight_ = Wmax_;
    }
    else
      weight_ = 0.0;
    wdev_ *= std::exp(-update_interval_/tau_deriv_);
    //if ( e.get_sender().get_gid() == 10 )
    //std::cout << "new weight: " << weight_ << ", derivative updated to " << wdev_ << std::endl;
    next_update_ += update_interval_;
  }

  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;

  // For a new synapse, t_lastspike contains the point in time of the last spike.
  // So we initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of history[0, ..., t_last_spike - dendritic_delay] have been
  // incremented by Archiving_Node::register_stdp_izh_connection(). See bug #218 for details.
  Node* target = get_target( t );
  target->get_history(t_lastspike + axonal_delay - dendritic_delay,
		       t_spike + axonal_delay - dendritic_delay,
		       &start, &finish);

  //facilitation due to post-synaptic spikes since last pre-synaptic spike
  double minus_dt;
  if (start != finish)
  {
    while(start != finish)
    {
      t_lastpostspike_ = start->t_;
      start++;
    }
    minus_dt = t_lastspike + axonal_delay - (t_lastpostspike_ + dendritic_delay);
    wdev_ += lambda_ * std::exp(minus_dt / tau_plus_);
    
    //if ( e.get_sender().get_gid() == 10 )
    //std::cout << "facilitation for dt=" << -minus_dt << " is " << wdev_ << std::endl;
  }

  //depression due to new pre-synaptic spike
  minus_dt = t_lastpostspike_ + dendritic_delay - (t_spike + axonal_delay);
  wdev_ -= alpha_ * lambda_ * std::exp(minus_dt / tau_minus_);

  //if ( e.get_sender().get_gid() == 10 )
  //std::cout << "depression for dt=" << -minus_dt << " is " << wdev_ << std::endl;

  e.set_receiver(*target);
  e.set_weight(weight_);
  e.set_delay( get_delay() );
  e.set_rport( get_rport() );
  e();

}


template < typename targetidentifierT >
STDPIzhConnection< targetidentifierT >::STDPIzhConnection()
  : ConnectionBase()
  ,  weight_(1.)
  ,  tau_plus_(20.0)
  ,  tau_minus_(20.0)
  ,  tau_deriv_(10.0)
  ,  lambda_(0.1)
  ,  alpha_(1.2)
  ,  wdev_(0.0)
  ,  dw0_(0.01)
  ,  Wmax_(10.0)
  ,  t_lastpostspike_(-1000.0)
  ,  next_update_(1000.0)
  ,  update_interval_(1000.)
{
}

template < typename targetidentifierT >
STDPIzhConnection< targetidentifierT >::STDPIzhConnection(
  const STDPIzhConnection< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  ,  tau_plus_(rhs.tau_plus_)
  ,  tau_minus_(rhs. tau_minus_)
  ,  tau_deriv_(rhs.tau_deriv_)
  ,  lambda_(rhs.lambda_)
  ,  alpha_(rhs.alpha_)
  ,  wdev_(rhs.wdev_)
  ,  dw0_(rhs.dw0_)
  ,  Wmax_(rhs.Wmax_)
  ,  t_lastpostspike_(rhs.t_lastpostspike_)
  ,  next_update_(rhs.next_update_)
  ,  update_interval_(rhs.update_interval_)
{
}

template < typename targetidentifierT >
void
STDPIzhConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, "tau_plus", tau_plus_);
  def< double >( d, "tau_minus", tau_minus_);
  def< double >( d, "tau_deriv", tau_deriv_);
  def< double >( d, "lambda", lambda_);
  def< double >( d, "alpha", alpha_);
  def< double >( d, "wdev", wdev_);
  def< double >( d, "dw0", dw0_);
  def< double >( d, "Wmax", Wmax_);
  def< double >( d, "update_interval", update_interval_);
}

template < typename targetidentifierT >
void
STDPIzhConnection< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, "tau_plus", tau_plus_);
  updateValue< double >( d, "tau_minus", tau_minus_);
  updateValue< double >( d, "tau_deriv", tau_deriv_);
  updateValue< double >( d, "lambda", lambda_);
  updateValue< double >( d, "alpha", alpha_);
  updateValue< double >( d, "wdev", wdev_);
  updateValue< double >( d, "dw0", dw0_);
  updateValue< double >( d, "Wmax", Wmax_);
  updateValue< double >( d, "update_interval", update_interval_);

}

} // of namespace nest

#endif // of #ifndef STDP_IZH_CONNECTION_H
