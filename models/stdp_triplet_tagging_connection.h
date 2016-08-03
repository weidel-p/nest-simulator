/*
*  stdp_triplet_tagging_connection.h
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

#ifndef STDP_TRIPLET_TAGGING_CONNECTION_H
#define STDP_TRIPLET_TAGGING_CONNECTION_H

/* BeginDocumentation
 Name: stdp_triplet_tagging_synapse - Synapse type with spike-timing dependent
                              plasticity (triplets).

 Description:
   stdp_triplet_synapse is a connection with spike time dependent
   plasticity accounting for spike triplet effects (as defined in [1]).

 STDP examples:
   pair-based   Aplus_triplet = Aminus_triplet = 0.0
   triplet      Aplus_triplet = Aminus_triplet = 1.0

 Parameters:
   tau_plus           double - time constant of short presynaptic trace
                             - (tau_plus of [1])
   tau_plus_triplet   double - time constant of long presynaptic trace
                             - (tau_x of [1])
   Aplus              double - weight of pair potentiation rule
                             - (A_plus_2 of [1])
   Aplus_triplet      double - weight of triplet potentiation rule
                             - (A_plus_3 of [1])
   Aminus             double - weight of pair depression rule
                               (A_minus_2 of [1])
   Aminus_triplet     double - weight of triplet depression rule
                             - (A_minus_3 of [1])
   Wmax               double - maximum allowed weight

 States:
   Kplus              double: pre-synaptic trace (r_1 of [1])
   Kplus_triplet      double: triplet pre-synaptic trace (r_2 of [1])

 Transmits: SpikeEvent

 References:
   [1] J.-P. Pfister & W. Gerstner (2006) Triplets of Spikes in a Model
       of Spike Timing-Dependent Plasticity.  The Journal of Neuroscience
       26(38):9673-9682; doi:10.1523/JNEUROSCI.1425-06.2006

 Notes:
   - Presynaptic traces r_1 and r_2 of [1] are stored in the connection as
     Kplus and Kplus_triplet and decay with time-constants tau_plus and
     tau_plus_triplet, respectively.
   - Postsynaptic traces o_1 and o_2 of [1] are acquired from the post-synaptic
     neuron states Kminus_ and triplet_Kminus_ which decay on time-constants
     tau_minus and tau_minus_triplet, respectively. These two time-constants
     can be set as properties of the postsynaptic neuron.
   - This version implements the 'all-to-all' spike interaction of [1]. The
     'nearest-spike' interaction of [1] can currently not be implemented
     without changing the postsynaptic archiving-node (clip the traces to a
     maximum of 1).

 FirstVersion: Nov 2007
 Author: Abigail Morrison, Eilif Muller, Alexander Seeholzer, Teo Stocco
 Adapted by: Philipp Weidel
 SeeAlso: stdp_triplet_synapse_hpc, synapsedict, stdp_synapse, static_synapse
*/

// C-header for math.h since copysign() is in C99 but not C++98
#include <math.h>

// Includes from models:
#include "volume_transmitter.h"

#include "connection.h"
#include "spikecounter.h"
#include <iostream>

namespace nest
{

class STDPTripletTaggingCommonProperties : public CommonSynapseProperties
{
public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  STDPTripletTaggingCommonProperties();

  /**
   * Get all properties and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  Node* get_node();

  long_t get_vt_gid() const;

  volume_transmitter* vt_;
  double_t A_plus_;
  double_t A_minus_;
  double_t tau_plus_;
  double_t tau_c_;
  double_t tau_n_;
  double_t b_;
  double_t Wmin_;
  double_t Wmax_;

  double_t tau_target_;
  double_t theta_n_;
  double_t theta_tag_;
};

inline long_t
STDPTripletTaggingCommonProperties::get_vt_gid() const
{
  if ( vt_ != 0 )
    return vt_->get_gid();
  else
    return -1;
}




// connections are templates of target identifier type
// (used for pointer / target index addressing)
// derived from generic connection template
template < typename targetidentifierT >
class STDPTripletTaggingConnection : public Connection< targetidentifierT >
{

public:
 typedef STDPTripletTaggingCommonProperties CommonPropertiesType;
 typedef Connection< targetidentifierT > ConnectionBase;

 /**
  * Default Constructor.
  * Sets default values for all parameters. Needed by GenericConnectorModel.
  */
 STDPTripletTaggingConnection();

 /**
  * Copy constructor.
  * Needs to be defined properly in order for GenericConnector to work.
  */
 STDPTripletTaggingConnection( const STDPTripletTaggingConnection& );

 /**
  * Default Destructor.
  */
 ~STDPTripletTaggingConnection()
 {
 }

 // Explicitly declare all methods inherited from the dependent base
 // ConnectionBase. This avoids explicit name prefixes in all places
 // these functions are used. Since ConnectionBase depends on the template
 // parameter, they are not automatically found in the base class.
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
   double_t t_lastspike,
   const STDPTripletTaggingCommonProperties& cp );

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
  * connections we have to call register_stdp_connection on the target neuron
  * to inform the Archiver to collect spikes for this connection.
  *
  * \param s The source node
  * \param r The target node
  * \param receptor_type The ID of the requested receptor type
  * \param t_lastspike last spike emitted by presynaptic neuron
  */
 void
 check_connection( Node& s,
   Node& t,
   rport receptor_type,
   double_t t_lastspike,
   const CommonPropertiesType& cp)
 {
   if ( cp.vt_ == 0 )
      throw BadProperty(
        "No volume transmitter has been assigned to the dopamine synapse." );

   ConnTestDummyNode dummy_target;

   ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

   t.register_stdp_connection( t_lastspike - get_delay() );
 }

 void
 set_weight( double_t w )
 {
   weight_ = w;
 }

 void trigger_update_weight( thread t,
   const std::vector< spikecounter >& dopa_spikes,
   double_t t_trig,
   const STDPTripletTaggingCommonProperties& cp );



private:

 // update dopamine trace from last to current dopamine spike and increment
 // index
 void update_dopamine_( const std::vector< spikecounter >& dopa_spikes,
   const STDPTripletTaggingCommonProperties& cp );


 void process_dopa_spikes_( const std::vector< spikecounter >& dopa_spikes,
   double_t t0,
   double_t t1,
   const STDPTripletTaggingCommonProperties& cp );

  void update_weight_( double_t c0,
    double_t n0,
    double_t minus_dt,
    const STDPTripletTaggingCommonProperties& cp );

 inline double_t
 facilitate_( double_t w, double_t kplus, double_t ky ,
 const STDPTripletTaggingCommonProperties& cp)
 {
   double_t new_w = std::abs( w ) + kplus * ( Aplus_ + Aplus_triplet_ * ky );
   if (new_w - w > cp.theta_tag_){
       tag_ = new_w;
   }
   return copysign( new_w < std::abs( cp.Wmax_ ) ? new_w : cp.Wmax_, cp.Wmax_ );
 }

 inline double_t
 depress_( double_t w, double_t kminus, double_t Kplus_triplet_ ,
 const STDPTripletTaggingCommonProperties& cp)
 {
   double_t new_w =
     std::abs( w ) - kminus * ( Aminus_ + Aminus_triplet_ * Kplus_triplet_ );
   if (w - new_w > cp.theta_tag_){
       tag_ = new_w;
    }
   return copysign( new_w > 0.0 ? new_w : 0.0, cp.Wmax_ );
 }

 inline void 
 set_target_( )
 {
     target_ = tag_;
 }

 inline void 
 reset_tag_( )
 {
     tag_ = target_;
 }

 // data members of each connection
 double_t weight_;
 double_t tau_plus_;
 double_t tau_plus_triplet_;
 double_t Aplus_;
 double_t Aminus_;
 double_t Aplus_triplet_;
 double_t Aminus_triplet_;
 double_t Kplus_;
 double_t Kplus_triplet_;

 double_t n_;
 double_t c_;

 double_t tag_;
 double_t target_;
 

 // dopa_spikes_idx_ refers to the dopamine spike that has just been processes
 // after trigger_update_weight a pseudo dopamine spike at t_trig is stored at
 // index 0 and dopa_spike_idx_ = 0
 index dopa_spikes_idx_;

 // time of last update, which is either time of last presyn. spike or
 // time-driven update
 double_t t_last_update_;
};

template < typename targetidentifierT >
inline void
STDPTripletTaggingConnection< targetidentifierT >::update_dopamine_(
  const std::vector< spikecounter >& dopa_spikes,
  const STDPTripletTaggingCommonProperties& cp )
{
  double_t minus_dt = dopa_spikes[ dopa_spikes_idx_ ].spike_time_
    - dopa_spikes[ dopa_spikes_idx_ + 1 ].spike_time_;
  ++dopa_spikes_idx_;
  n_ = n_ * std::exp( minus_dt / cp.tau_n_ )
    + dopa_spikes[ dopa_spikes_idx_ ].multiplicity_ / cp.tau_n_;
}

template < typename targetidentifierT >
inline void
STDPTripletTaggingConnection< targetidentifierT >::process_dopa_spikes_(
  const std::vector< spikecounter >& dopa_spikes,
  double_t t0,
  double_t t1,
  const STDPTripletTaggingCommonProperties& cp )
{
  // process dopa spikes in (t0, t1]
  // propagate weight from t0 to t1
  if ( ( dopa_spikes.size() > dopa_spikes_idx_ + 1 )
    && ( dopa_spikes[ dopa_spikes_idx_ + 1 ].spike_time_ <= t1 ) )
  {
    // there is at least 1 dopa spike in (t0, t1]
    // propagate weight up to first dopa spike and update dopamine trace
    // weight and eligibility c are at time t0 but dopamine trace n is at time
    // of last dopa spike
    double_t n0 =
      n_ * std::exp( ( dopa_spikes[ dopa_spikes_idx_ ].spike_time_ - t0 )
             / cp.tau_n_ ); // dopamine trace n at time t0
    update_weight_(
      c_, n0, t0 - dopa_spikes[ dopa_spikes_idx_ + 1 ].spike_time_, cp );
    update_dopamine_( dopa_spikes, cp );

    // process remaining dopa spikes in (t0, t1]
    double_t cd;
    while ( ( dopa_spikes.size() > dopa_spikes_idx_ + 1 )
      && ( dopa_spikes[ dopa_spikes_idx_ + 1 ].spike_time_ <= t1 ) )
    {
      // propagate weight up to next dopa spike and update dopamine trace
      // weight and dopamine trace n are at time of last dopa spike td but
      // eligibility c is at time
      // t0
      cd = c_ * std::exp( ( t0 - dopa_spikes[ dopa_spikes_idx_ ].spike_time_ )
                  / cp.tau_c_ ); // eligibility c at time of td
      update_weight_( cd,
        n_,
        dopa_spikes[ dopa_spikes_idx_ ].spike_time_
          - dopa_spikes[ dopa_spikes_idx_ + 1 ].spike_time_,
        cp );
      update_dopamine_( dopa_spikes, cp );
    }

    // propagate weight up to t1
    // weight and dopamine trace n are at time of last dopa spike td but
    // eligibility c is at time t0
    cd = c_ * std::exp( ( t0 - dopa_spikes[ dopa_spikes_idx_ ].spike_time_ )
                / cp.tau_c_ ); // eligibility c at time td
    update_weight_(
      cd, n_, dopa_spikes[ dopa_spikes_idx_ ].spike_time_ - t1, cp );
  }
  else
  {
    // no dopamine spikes in (t0, t1]
    // weight and eligibility c are at time t0 but dopamine trace n is at time
    // of last dopa spike

    double_t n0 =
      n_ * std::exp( ( dopa_spikes[ dopa_spikes_idx_ ].spike_time_ - t0 )
             / cp.tau_n_ ); // dopamine trace n at time t0
    update_weight_( c_, n0, t0 - t1, cp );
  }

  // update eligibility trace c for interval (t0, t1]
  c_ = c_ * std::exp( ( t0 - t1 ) / cp.tau_c_ );
}

template < typename targetidentifierT >
inline void
STDPTripletTaggingConnection< targetidentifierT >::trigger_update_weight( thread t,
  const std::vector< spikecounter >& dopa_spikes,
  const double_t t_trig,
  const STDPTripletTaggingCommonProperties& cp )
{
  // propagate all state variables to time t_trig
  // this does not include the depression trace K_minus, which is updated in the
  // postsyn. neuron

  // purely dendritic delay
  double_t dendritic_delay = get_delay();

  // get spike history in relevant range (t_last_update, t_trig] from postsyn.
  // neuron
  std::deque< histentry >::iterator start;
  std::deque< histentry >::iterator finish;
  get_target( t )->get_history( t_last_update_ - dendritic_delay,
    t_trig - dendritic_delay,
    &start,
    &finish );

  // facilitation due to postsyn. spikes since last update
  double_t t0 = t_last_update_;
  double_t minus_dt;
  while ( start != finish )
  {
    process_dopa_spikes_( dopa_spikes, t0, start->t_ + dendritic_delay, cp );
    t0 = start->t_ + dendritic_delay;
    minus_dt = t_last_update_ - t0;
    //facilitate_( Kplus_ * std::exp( minus_dt / cp.tau_plus_ ), cp );
    ++start;
  }

  // propagate weight, eligibility trace c, dopamine trace n and facilitation
  // trace K_plus to time t_trig but do not increment/decrement as there are no
  // spikes to be handled at t_trig
  process_dopa_spikes_( dopa_spikes, t0, t_trig, cp );
  n_ = n_ * std::exp( ( dopa_spikes[ dopa_spikes_idx_ ].spike_time_ - t_trig )
              / cp.tau_n_ );
  //Kplus_ = Kplus_ * std::exp( ( t_last_update_ - t_trig ) / cp.tau_plus_ );

  t_last_update_ = t_trig;
  dopa_spikes_idx_ = 0;
}



template < typename targetidentifierT >
inline void
STDPTripletTaggingConnection< targetidentifierT >::update_weight_( double_t c0,
  double_t n0,
  double_t minus_dt,
  const STDPTripletTaggingCommonProperties& cp )
{
  //weight_ = weight_
  //  - c0 * ( n0 / taus_ * numerics::expm1( taus_ * minus_dt )
  //           - cp.b_ * cp.tau_c_ * numerics::expm1( minus_dt / cp.tau_c_ ) );
  //
  
  if (n_ > cp.b_ + cp.theta_n_)
      set_target_();
  else if (n_ < cp.b_ - cp.theta_n_)
      reset_tag_();


  weight_ = weight_ + ((weight_ - target_) / cp.tau_target_) * minus_dt;
  tag_ = tag_ + ((tag_ - target_) / (cp.tau_target_ / 2.)) * minus_dt;

  if ( weight_ < cp.Wmin_ )
    weight_ = cp.Wmin_;
  if ( weight_ > cp.Wmax_ )
    weight_ = cp.Wmax_;
}




/**
* Send an event to the receiver of this connection.
* \param e The event to send
* \param t The thread on which this connection is stored.
* \param t_lastspike Time point of last spike emitted
* \param cp Common properties object, containing the stdp parameters.
*/
template < typename targetidentifierT >
inline void
STDPTripletTaggingConnection< targetidentifierT >::send( Event& e,
 thread t,
 double_t t_lastspike,
 const STDPTripletTaggingCommonProperties& cp)
{

 double_t t_spike = e.get_stamp().get_ms();
 double_t dendritic_delay = get_delay();
 Node* target = get_target( t );

 // get spike history in relevant range (t1, t2] from post-synaptic neuron
 std::deque< histentry >::iterator start;
 std::deque< histentry >::iterator finish;
 target->get_history(
   t_lastspike - dendritic_delay, t_spike - dendritic_delay, &start, &finish );

 // facilitation due to post-synaptic spikes since last pre-synaptic spike
 while ( start != finish )
 {
   // post-synaptic spike is delayed by dendritic_delay so that
   // it is effectively late by that much at the synapse.
   double_t minus_dt = t_lastspike - ( start->t_ + dendritic_delay );

   // subtract 1.0 yields the triplet_Kminus value just prior to
   // the post synaptic spike, implementing the t-epsilon in
   // Pfister et al, 2006
   double_t ky = start->triplet_Kminus_ - 1.0;
   ++start;
   if ( minus_dt == 0 )
   {
     continue;
   }

   weight_ =
     facilitate_( weight_, Kplus_ * std::exp( minus_dt / tau_plus_ ), ky, cp);
 }

 // depression due to new pre-synaptic spike
 Kplus_triplet_ *= std::exp( ( t_lastspike - t_spike ) / tau_plus_triplet_ );

 // dendritic delay means we must look back in time by that amount
 // for determining the K value, because the K value must propagate
 // out to the synapse
 weight_ = depress_(
   weight_, target->get_K_value( t_spike - dendritic_delay ), Kplus_triplet_ , cp);

 Kplus_triplet_ += 1.0;
 Kplus_ = Kplus_ * std::exp( ( t_lastspike - t_spike ) / tau_plus_ ) + 1.0;

 e.set_receiver( *target );
 e.set_weight( weight_ );
 e.set_delay( get_delay_steps() );
 e.set_rport( get_rport() );
 e();
}

// Defaults come from reference [1] data fitting and table 3.
template < typename targetidentifierT >
STDPTripletTaggingConnection< targetidentifierT >::STDPTripletTaggingConnection()
 : ConnectionBase()
 , weight_( 1.0 )
 , tau_plus_( 16.8 )
 , tau_plus_triplet_( 101.0 )
 , Aplus_( 5e-10 )
 , Aminus_( 7e-3 )
 , Aplus_triplet_( 6.2e-3 )
 , Aminus_triplet_( 2.3e-4 )
 , Kplus_( 0.0 )
 , Kplus_triplet_( 0.0 )
 , n_(0.)
 , c_(0.)
 , dopa_spikes_idx_( 0 )
 , t_last_update_( 0.0 )
 , tag_( 1.0 )
 , target_( 1.0 )
{
}

template < typename targetidentifierT >
STDPTripletTaggingConnection< targetidentifierT >::STDPTripletTaggingConnection(
 const STDPTripletTaggingConnection< targetidentifierT >& rhs )
 : ConnectionBase( rhs )
 , weight_( rhs.weight_ )
 , tau_plus_( rhs.tau_plus_ )
 , tau_plus_triplet_( rhs.tau_plus_triplet_ )
 , Aplus_( rhs.Aplus_ )
 , Aminus_( rhs.Aminus_ )
 , Aplus_triplet_( rhs.Aplus_triplet_ )
 , Aminus_triplet_( rhs.Aminus_triplet_ )
 , Kplus_( rhs.Kplus_ )
 , Kplus_triplet_( rhs.Kplus_triplet_ )
 , n_( rhs.n_ )
 , c_( rhs.c_ )
 , dopa_spikes_idx_( rhs.dopa_spikes_idx_ )
 , t_last_update_( rhs.t_last_update_ )
 , tag_( rhs.tag_)
 , target_( rhs.target_)
{
}

template < typename targetidentifierT >
void
STDPTripletTaggingConnection< targetidentifierT >::get_status(
 DictionaryDatum& d ) const
{
 ConnectionBase::get_status( d );
 def< double_t >( d, names::weight, weight_ );
 def< double_t >( d, "tau_plus", tau_plus_ );
 def< double_t >( d, "tau_plus_triplet", tau_plus_triplet_ );
 def< double_t >( d, "Aplus", Aplus_ );
 def< double_t >( d, "Aminus", Aminus_ );
 def< double_t >( d, "Aplus_triplet", Aplus_triplet_ );
 def< double_t >( d, "Aminus_triplet", Aminus_triplet_ );
 def< double_t >( d, "Kplus", Kplus_ );
 def< double_t >( d, "Kplus_triplet", Kplus_triplet_ );

 // own properties, different for individual synapse
 def< double_t >( d, "n", n_ );
 def< double_t >( d, "c", c_ );
 def< double_t >( d, "tag", tag_);
 def< double_t >( d, "target", target_ );
}

template < typename targetidentifierT >
void
STDPTripletTaggingConnection< targetidentifierT >::set_status(
 const DictionaryDatum& d,
 ConnectorModel& cm )
{
 ConnectionBase::set_status( d, cm );
 updateValue< double_t >( d, names::weight, weight_ );
 updateValue< double_t >( d, "tau_plus", tau_plus_ );
 updateValue< double_t >( d, "tau_plus_triplet", tau_plus_triplet_ );
 updateValue< double_t >( d, "Aplus", Aplus_ );
 updateValue< double_t >( d, "Aminus", Aminus_ );
 updateValue< double_t >( d, "Aplus_triplet", Aplus_triplet_ );
 updateValue< double_t >( d, "Aminus_triplet", Aminus_triplet_ );
 updateValue< double_t >( d, "Kplus", Kplus_ );
 updateValue< double_t >( d, "Kplus_triplet", Kplus_triplet_ );

 updateValue< double_t >( d, "n", n_ );
 updateValue< double_t >( d, "c", c_ );
 updateValue< double_t >( d, "tag", tag_);
 updateValue< double_t >( d, "target", target_ );

 if ( not( Kplus_ >= 0 ) )
 {
   throw BadProperty( "State Kplus must be positive." );
 }

 if ( not( Kplus_triplet_ >= 0 ) )
 {
   throw BadProperty( "State Kplus_triplet must be positive." );
 }
}

} // of namespace nest

#endif // of #ifndef STDP_TRIPLET_TAGGING_CONNECTION_H
