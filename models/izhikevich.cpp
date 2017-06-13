/*
 *  izhikevich.cpp
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

#include "izhikevich.h"

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from models:
#include "stdp_izh_connection.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"
#include "target_identifier.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::izhikevich > nest::izhikevich::recordablesMap_;

namespace nest
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< izhikevich >::create()
{
  // use standard names whereever you can for consistency!
  insert_( names::V_m, &izhikevich::get_V_m_ );
  insert_( names::U_m, &izhikevich::get_U_m_ );
  insert_( names::I, &izhikevich::get_I_ );
  insert_( "combined_current", &izhikevich::get_combined_current_ );


}
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::izhikevich::Parameters_::Parameters_()
  : a_( 0.02 )                                      // a
  , b_( 0.2 )                                       // b
  , c_( -65.0 )                                     // c without unit
  , d_( 8.0 )                                       // d
  , I_e_( 0.0 )                                     // pA
  , V_th_( 30.0 )                                   // mV
  , V_min_( -std::numeric_limits< double >::max() ) // mV
  , consistent_integration_( true )
{
}

nest::izhikevich::State_::State_()
  : v_( -65.0 ) // membrane potential
  , u_( 0.0 )   // membrane recovery variable
  , I_( 0.0 )   // input current
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::izhikevich::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::I_e, I_e_ );
  def< double >( d, names::V_th, V_th_ ); // threshold value
  def< double >( d, names::V_min, V_min_ );
  def< double >( d, names::a, a_ );
  def< double >( d, names::b, b_ );
  def< double >( d, names::c, c_ );
  def< double >( d, names::d, d_ );
  def< bool >( d, names::consistent_integration, consistent_integration_ );
}

void
nest::izhikevich::Parameters_::set( const DictionaryDatum& d )
{

  updateValue< double >( d, names::V_th, V_th_ );
  updateValue< double >( d, names::V_min, V_min_ );
  updateValue< double >( d, names::I_e, I_e_ );
  updateValue< double >( d, names::a, a_ );
  updateValue< double >( d, names::b, b_ );
  updateValue< double >( d, names::c, c_ );
  updateValue< double >( d, names::d, d_ );
  updateValue< bool >(
    d, names::consistent_integration, consistent_integration_ );
  const double h = Time::get_resolution().get_ms();
  if ( not consistent_integration_ && h != 1.0 )
  {
    LOG(
      M_INFO, "Parameters_::set", "Use 1.0 ms as resolution for consistency." );
  }
}

void
nest::izhikevich::State_::get( DictionaryDatum& d, const Parameters_& ) const
{
  def< double >( d, names::U_m, u_ ); // Membrane potential recovery variable
  def< double >( d, names::V_m, v_ ); // Membrane potential
}

void
nest::izhikevich::State_::set( const DictionaryDatum& d, const Parameters_& )
{
  updateValue< double >( d, names::U_m, u_ );
  updateValue< double >( d, names::V_m, v_ );
  updateValue< double >( d, names::I, I_ );

}

nest::izhikevich::Buffers_::Buffers_( izhikevich& n )
  : logger_( n )
{
  post_spikes_.clear();
  post_spikes_.push_back( -10000 );
}

nest::izhikevich::Buffers_::Buffers_( const Buffers_&, izhikevich& n )
  : logger_( n )
{
  post_spikes_.clear();
  post_spikes_.push_back( -10000 );
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::izhikevich::izhikevich()
  : Archiving_Node()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

nest::izhikevich::izhikevich( const izhikevich& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::izhikevich::init_state_( const Node& proto )
{
  const izhikevich& pr = downcast< izhikevich >( proto );
  S_ = pr.S_;
}

void
nest::izhikevich::init_buffers_()
{
  B_.spikes_.clear();          // includes resize
  B_.stdp_izh_spikes_.resize();
  B_.stdp_izh_spikes_.clear();
  B_.currents_.clear();        // includes resize
  B_.logger_.reset();          // includes resize
  Archiving_Node::clear_history();
}

void
nest::izhikevich::calibrate()
{
  B_.logger_.init();
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 */

void
nest::izhikevich::update( Time const& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );
  
  // at start of slice, tell input queue to prepare for delivery
  if ( from == 0 )
    B_.stdp_izh_spikes_.prepare_delivery();

  const double h = Time::get_resolution().get_ms();
  double v_old, u_old;

  for ( long lag = from; lag < to; ++lag )
  {
    // summed contribution of non-STDPIzhConnection spikes
    double I_syn = B_.spikes_.get_value( lag );

    // contributions of STDPIzhConnection spikes
    const long T = origin.get_steps() + lag; // time at start of update step
    double ev_offset;                        // not used here
    double ev_weight;                        // pointer to synapse that requires casting
    bool end_of_refract;                     // not used here
    while ( B_.stdp_izh_spikes_.get_next_spike( T, ev_offset, ev_weight, end_of_refract ) )
    {
      const STDPIzhConnection* syn = reinterpret_cast< STDPIzhConnection* >( static_cast< long >( ev_weight ) );
      I_syn += syn->get_weight();
    }

    S_.combined_current_= S_.I_+I_syn;

    //S_.combined_current_ = round(S_.combined_current_ * 100000000.) / 100000000.;

//    int spike=0;
//    // threshold crossing
//    if ( S_.v_ >= P_.V_th_ )
//    {
//      S_.v_ = P_.c_;
//      S_.u_ = S_.u_ + P_.d_;
//       spike=1;
//     }


    // neuron is never refractory
    // use standard forward Euler numerics in this case
    if ( P_.consistent_integration_ )
    {
      v_old = S_.v_;
      u_old = S_.u_;
      S_.v_ += h * ( 0.04 * v_old * v_old + 5.0 * v_old + 140.0 - u_old + S_.I_
                     + P_.I_e_ ) + I_syn;
      S_.u_ += h * P_.a_ * ( P_.b_ * v_old - u_old );
    }
    // use numerics published in Izhikevich (2003) in this case (not
    // recommended)
    else
    {
      //v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
      //S_.v_ += h * 0.5 * ( (0.04 * S_.v_ +5.) * S_.v_  + 140.0 - S_.u_
      //                     + S_.I_ + P_.I_e_ + I_syn );
      //S_.v_ += h * 0.5 * ( (0.04 * S_.v_ +5.) * S_.v_  + 140.0 - S_.u_
      //                     + S_.I_ + P_.I_e_ + I_syn );
      //S_.u_ += h * P_.a_ * ( P_.b_ * S_.v_ - S_.u_ );

      S_.v_ += 0.5 * ( (0.04 * S_.v_ +5) * S_.v_  + 140. - S_.u_ + S_.combined_current_ );
      S_.v_ += 0.5 * ( (0.04 * S_.v_ +5) * S_.v_  + 140. - S_.u_ + S_.combined_current_ );
      S_.u_ += P_.a_ * ( 0.2 * S_.v_ - S_.u_ );
    }

//				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140.-u[i]+I[i]);
//				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140.-u[i]+I[i]);
//
////				u = u + a.*(0.2*v-u);
//				u[i]+=a[i]*(0.2*v[i]-u[i]);
    //S_.v_ = round(S_.v_ * 10000.) / 10000.;
    //S_.u_ = round(S_.u_ * 10000.) / 10000.;

    // lower bound of membrane potential
    //S_.v_ = ( S_.v_ < P_.V_min_ ? P_.V_min_ : S_.v_ );
    if (S_.v_ < P_.V_min_)
        S_.v_ = P_.V_min_;

    
    // set new input current
    S_.I_ = B_.currents_.get_value( lag );


    // voltage logging
    B_.logger_.record_data( origin.get_steps() + lag );

    if ( (S_.v_ >= P_.V_th_ ))//|| (spike==1) )
    {

    // compute spike time
      set_spiketime( Time::step( origin.get_steps() + lag+1 ) );

      // keeps track of spike history for STDPIzhConnections
      // cleared every syn_update_interval
      B_.post_spikes_.push_back( Time(Time::step( origin.get_steps() + lag+1)).get_ms() );

      SpikeEvent se;
      kernel().event_delivery_manager.send( *this, se, lag );

      S_.v_ = P_.c_;
      S_.u_ = S_.u_ + P_.d_;

    }
  }
}

void
nest::izhikevich::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );
  
  int multiplicity = e.get_multiplicity();

  if ( multiplicity == -1 )
  {
    // multiplicity of -1 indicates a SpikeEvent from an STDPIzhConnection
    B_.stdp_izh_spikes_.add_spike(
      e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      e.get_stamp().get_steps() + e.get_delay() - 1,
      e.get_offset(),
      e.get_weight() );
    // As the event object maybe recycled, reset multiplicity to its default of 1.
    e.set_multiplicity(1);
  }
  else
  {
    // regular SpikeEvent
    B_.spikes_.add_value(
      e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
}

void
nest::izhikevich::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();
  B_.currents_.add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
    w * c );
}

void
nest::izhikevich::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}
