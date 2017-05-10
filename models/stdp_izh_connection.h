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
   dependent plasticity (as defined in [1]).

  Parameters:
   tau_plus   double - Time constant of STDP window, potentiation in ms
   tau_minus  double - Time constant of STDP window, depression in ms
   lambda     double - Step size
   alpha      double - Asymmetry parameter (scales depressing increments as
                       alpha*lambda)
   Wmax       double - Maximum allowed weight

  Transmits: SpikeEvent

  References:
   [1] Izhikevich EM (2006) Polychronization: computation with spikes.
       Neural Comput. 18(2):245-82.

  Author: Abigail Morrison
  Adapted by: Philipp Weidel, Susanne Kunkel
  SeeAlso: synapsedict, stdp_synapse
*/

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"
#include "target_identifier.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"


namespace nest
{

class STDPIzhConnection : public Connection< TargetIdentifierPtrRport >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< TargetIdentifierPtrRport > ConnectionBase;

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
   * \param tid thread of this synapse
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses (empty).
   */
  void send( Event& e,
    thread tid,
    double t_lastspike,
    const CommonPropertiesType& cp );

  bool requires_time_driven_update() const
  {
    return true;
  }

  void time_driven_update( const thread tid, const double t_trig, const CommonPropertiesType& );

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
  check_connection( Node& source,
    Node& target,
    rport receptor_type,
    double t_lastspike,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, source, target, receptor_type );

    target.register_stdp_connection( t_lastspike - get_delay() );
  }

  double
  get_weight() const
  {
    return weight_;
  }

  void
  set_weight( const double weight )
  {
    weight_ = weight;
  }

private:

  // data members of each connection
  double weight_;
  double wdev_;
  double K_plus_;
  double K_minus_;
  double tau_plus_;
  double tau_minus_;
  double lambda_;
  double alpha_;
  double Wmax_;
  double t_last_update_;
  double t_last_post_spike_;
  bool   consistent_integration_;
  std::vector< double > pre_spikes_;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send.
 * \param tid The thread on which this connection is stored.
 */
inline void
STDPIzhConnection::send( Event& e,
  thread tid,
  double,
  const CommonPropertiesType& )
{
  // The update of the synaptic weight is implemented in
  // time_driven_update(), which is triggered in intervals
  // of one second (syn_update_interval can be adjusted).

  // keep track of the current presynaptic spike
  // add delay to spike time to account for purely axonal delay
  pre_spikes_.push_back( e.get_stamp().get_ms() + get_delay() );

  Node* target = get_target( tid );

  // send an event to the postsynaptic neuron
  e.set_receiver(*target);
  // This is a work-around in order to enable purely axonal delays:
  // As postsynaptic spikes that are relevant for the update
  // might not yet be available, do not send  the current synaptic
  // weight but a pointer to this synapse. The postsynaptic
  // neuron retrieves the correct synaptic weight later when the
  // spike is due at time t_pre_spike + delay.
  e.set_weight( static_cast< double >( reinterpret_cast< long >(this) ) );
  // Default multiplicity is 1
  // Use multiplicity -1 to signal to postsynaptic neuron that this
  // event is delivered through an STDPIzhConnection
  e.set_multiplicity( -1 );
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
}

} // of namespace nest

#endif // of #ifndef STDP_IZH_CONNECTION_H
