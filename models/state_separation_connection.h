/*
*  state_separation_connection.h
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

#ifndef STATE_SEPARATION_H
#define STATE_SEPARATION_H

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
#include <cmath>

namespace nest
{

/**
* Class containing the common properties for all synapses of type dopamine
* connection.
*/
class StateSeparationCommonProperties : public CommonSynapseProperties
{
public:
 /**
  * Default constructor.
  * Sets all property values to defaults.
  */
 StateSeparationCommonProperties();

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
 double A_;
 double tau_short_;
 double tau_long_;
 double tau_n_;
 double b_;
 double n_threshold_;
 double Wmin_;
 double Wmax_;
 double LTD_scaling_;
 double tau_decay_;
 double weight0_;
};

inline long
StateSeparationCommonProperties::get_vt_gid() const
{
 if ( vt_ != 0 )
   return vt_->get_gid();
 else
   return -1;
}

/**
* Class representing an StateSeparationConnection with homogeneous parameters,
* i.e. parameters are the same for all synapses.
*/
template < typename targetidentifierT >
class StateSeparationConnection : public Connection< targetidentifierT >
{

public:
 typedef StateSeparationCommonProperties CommonPropertiesType;
 typedef Connection< targetidentifierT > ConnectionBase;

 /**
  * Default Constructor.
  * Sets default values for all parameters. Needed by GenericConnectorModel.
  */
 StateSeparationConnection();

 /**
  * Copy constructor from a property object.
  * Needs to be defined properly in order for GenericConnector to work.
  */
 StateSeparationConnection( const StateSeparationConnection& );

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
 void send( Event& e, thread t, double_t, const StateSeparationCommonProperties& cp );

 void trigger_update_weight( thread t,
   const std::vector< spikecounter >& dopa_spikes,
   double t_trig,
   const StateSeparationCommonProperties& cp );

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
   const StateSeparationCommonProperties& cp);

 void process_next_(
   thread t,
   double t0,
   double t1,
   const StateSeparationCommonProperties& cp );


 // data members of each connection
 double weight_;
 double Kplus_short_;
 double Kplus_long_;
 double n_;


 // time of last update, which is either time of last presyn. spike or
 // time-driven update
 double t_last_update_;
 int dopa_spikes_idx_;
};

//
// Implementation of class StateSeparationConnection.
//

template < typename targetidentifierT >
StateSeparationConnection< targetidentifierT >::StateSeparationConnection()
 : ConnectionBase()
 , weight_( 1.0 )
 , Kplus_short_( 0.0 )
 , Kplus_long_( 0.0 )
 , n_( 0.0 )
 , dopa_spikes_idx_( 0 )
 , t_last_update_( 0.0 )
{
}

template < typename targetidentifierT >
StateSeparationConnection< targetidentifierT >::StateSeparationConnection(
 const StateSeparationConnection& rhs )
 : ConnectionBase( rhs )
 , weight_( rhs.weight_ )
 , Kplus_short_( rhs.Kplus_short_ )
 , Kplus_long_( rhs.Kplus_long_ )
 , n_( rhs.n_ )
 , dopa_spikes_idx_( rhs.dopa_spikes_idx_ )
 , t_last_update_( rhs.t_last_update_ )
{
}

template < typename targetidentifierT >
void
StateSeparationConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

 // base class properties, different for individual synapse
 ConnectionBase::get_status( d );
 def< double >( d, names::weight, weight_ );

 // own properties, different for individual synapse
 def< double >( d, "n", n_ );
 def< double >( d, "Kplus_short", Kplus_short_ );
 def< double >( d, "Kplus_long", Kplus_long_ );
}

template < typename targetidentifierT >
void
StateSeparationConnection< targetidentifierT >::set_status( const DictionaryDatum& d,
 ConnectorModel& cm )
{
 // base class properties
 ConnectionBase::set_status( d, cm );
 updateValue< double >( d, names::weight, weight_ );

 updateValue< double >( d, "n", n_ );
 updateValue< double >( d, "Kplus_short", Kplus_short_ );
 updateValue< double >( d, "Kplus_long", Kplus_long_ );
}




template < typename targetidentifierT >
inline void
StateSeparationConnection< targetidentifierT >::process_next_(
 thread t,
 double t0,
 double t1,
 const StateSeparationCommonProperties& cp )
{
  //std::cout << t0 << " " << t1 << std::endl;

  assert(t1 >= t0);
  Node* target = get_target( t );

  // update all traces to  t0
  double Kminus_short = target->get_firing_rate_short( t0 );
  //double Kminus_long = target->get_firing_rate_long( t0 );


  // calculate decay before weight update
  if (cp.tau_decay_ > 0.){
    if (weight_ > cp.weight0_) {
        weight_ -= (weight_ - cp.weight0_) * ( 1 - std::exp( (t0 - t1) / cp.tau_decay_ ) );

        if (weight_ < cp.weight0_)
            weight_ = cp.weight0_;
    }
    else{
        weight_ -= (weight_ - cp.weight0_) * ( 1 - std::exp( (t0 - t1) / cp.tau_decay_ ) );

        if (weight_ > cp.weight0_)
            weight_ = cp.weight0_;
    }
  }


  double n_diff = n_ - cp.n_threshold_;

  // update weight 
  //double dw = cp.A_ * (Kplus_short_ * Kminus_short - Kplus_long_ * Kminus_long) * 
  //                    std::pow(n_diff, 2) * (t1 - t0); 
  
  double dw = 0;
  if (Kminus_short > cp.b_){
      dw = cp.A_ * Kplus_short_* n_diff;
  }

  if (dw > 0){
      weight_ += dw;
  }
  else{
      weight_ += dw * cp.LTD_scaling_;
  }

  //if (dw > 0){
  //    weight_ += (1 - (weight_ - cp.Wmin_) / (cp.Wmax_ - cp.Wmin_) ) * dw;
  //}
  //else{
  //    weight_ += ((weight_ - cp.Wmin_) / (cp.Wmax_ - cp.Wmin_) ) * dw;
  //}


 if ( weight_ > cp.Wmax_ ){
   weight_ = cp.Wmax_;
 }
 if ( weight_ < cp.Wmin_){
   weight_ = cp.Wmin_;
 }

 // propagate Kplus_, Kminus_, n_ to t1 
 Kplus_short_ *= std::exp( ( t0 - t1 ) / cp.tau_short_ );
 //Kplus_long_ *= std::exp( ( t0 - t1 ) / cp.tau_long_ );


 n_ = n_ * std::exp( ( t0 - t1) / cp.tau_n_ );

}



template < typename targetidentifierT >
inline void
StateSeparationConnection< targetidentifierT >::process_dopa_spikes_(
 thread t,
 double t1,
 const std::vector< spikecounter >& dopa_spikes,
 const StateSeparationCommonProperties& cp )
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
StateSeparationConnection< targetidentifierT >::send( Event& e,
 thread t,
 double_t,
 const StateSeparationCommonProperties& cp )
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



 double Kminus_short = target->get_firing_rate_short( t_spike );
 //double Kminus_long = target->get_firing_rate_long( t_spike );

 e.set_receiver( *target );
 e.set_weight( weight_ );
 e.set_delay( get_delay_steps() );
 e.set_rport( get_rport() );
 e.set_Kplus_short( Kplus_short_ );
 e.set_Kplus_long( Kplus_short_ ); //TODO fix
 e.set_Kminus_short( Kminus_short );
 e.set_Kminus_long( Kminus_short ); //TODO fix
 e.set_dopa( n_ );
 e();

 Kplus_short_ += 1./cp.tau_short_;
 //Kplus_long_ += 1./cp.tau_long_;
 t_last_update_ = t_spike;
}

template < typename targetidentifierT >
inline void
StateSeparationConnection< targetidentifierT >::trigger_update_weight( thread t,
 const std::vector< spikecounter >& dopa_spikes,
 const double t_trig,
 const StateSeparationCommonProperties& cp )
{

 process_dopa_spikes_(t, t_trig, dopa_spikes, cp);
 process_next_(t, t_last_update_, t_trig, cp);
 t_last_update_ = t_trig;


}





} // of namespace nest

#endif // of #ifndef STATE_SEPARATION_H
