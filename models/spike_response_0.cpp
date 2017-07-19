/*
 *  spike_response_0.cpp
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

#include "spike_response_0.h"

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"
#include "propagator_stability.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::spike_response_0 > nest::spike_response_0::recordablesMap_;

namespace nest
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< spike_response_0 >::create()
{
  // use standard names whereever you can for consistency!
  insert_( names::V_m, &spike_response_0::get_V_m_ );
  insert_( names::weighted_spikes_ex, &spike_response_0::get_weighted_spikes_ex_ );
  insert_( names::weighted_spikes_in, &spike_response_0::get_weighted_spikes_in_ );
  insert_( names::I_syn_ex, &spike_response_0::get_I_syn_ex_ );
  insert_( names::I_syn_in, &spike_response_0::get_I_syn_in_ );
}
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::spike_response_0::Parameters_::Parameters_()
  : tau_mem_( 10.0 )             // ms
  , tau_syn_( 2.0 ) // ms
  , xi_(-0.005)
  , theta_(16.)
  , du_(2.)
  , rho0_(60.)
{
}

nest::spike_response_0::State_::State_()
  : i_0_( 0.0 )
  , i_syn_ex_( 0.0 )
  , i_syn_in_( 0.0 )
  , V_m_( 0.0 )
  , t_last_spike_( -10000.0 )
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::spike_response_0::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::tau_m, tau_mem_ );
  def< double >( d, names::tau_syn, tau_syn_ );
  def< double >( d, "xi", xi_);
  def< double >( d, names::theta, theta_);
  def< double >( d, "du", du_);
  def< double >( d, "rho0", rho0_);
}

double
nest::spike_response_0::Parameters_::set( const DictionaryDatum& d )
{

  updateValue< double >( d, names::tau_m, tau_mem_ );
  updateValue< double >( d, names::tau_syn, tau_syn_ );
  updateValue< double >( d, "xi", xi_);
  updateValue< double >( d, names::theta, theta_);
  updateValue< double >( d, "du", du_);
  updateValue< double >( d, "rho0", rho0_);
  return 0.;
}

void
nest::spike_response_0::State_::get( DictionaryDatum& d, const Parameters_& p ) const
{
  def< double >( d, names::V_m, V_m_ ); // Membrane potential
}

void
nest::spike_response_0::State_::set( const DictionaryDatum& d,
  const Parameters_& p,
  double delta_EL )
{
  if ( updateValue< double >( d, names::V_m, V_m_ ) )
  {
  }
  else
  {
  }
}

nest::spike_response_0::Buffers_::Buffers_( spike_response_0& n )
  : logger_( n )
{
}

nest::spike_response_0::Buffers_::Buffers_( const Buffers_&, spike_response_0& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::spike_response_0::spike_response_0()
  : Archiving_Node()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

nest::spike_response_0::spike_response_0( const spike_response_0& n )
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
nest::spike_response_0::init_state_( const Node& proto )
{
  const spike_response_0& pr = downcast< spike_response_0 >( proto );
  S_ = pr.S_;
}

void
nest::spike_response_0::init_buffers_()
{
  B_.spikes_ex_.clear(); // includes resize
  B_.spikes_in_.clear(); // includes resize
  B_.currents_.clear();  // includes resize
  B_.logger_.reset();
  Archiving_Node::clear_history();
}

void
nest::spike_response_0::calibrate()
{
  B_.currents_.resize( 2 );
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  const double h = Time::get_resolution().get_ms();

  // numbering of state vaiables: i_0 = 0, i_syn_ = 1, V_m_ = 2

  // commented out propagators: forward Euler
  // needed to exactly reproduce Tsodyks network

  // these P are independent
  V_.P11ex_ = std::exp( -h / P_.tau_syn_ );
  // P11ex_ = 1.0-h/tau_syn_;

  V_.P11in_ = std::exp( -h / P_.tau_syn_ );
  // P11in_ = 1.0-h/tau_syn_;

  V_.P22_ = std::exp( -h / P_.tau_mem_ );
  // P22_ = 1.0-h/tau_mem_;

  // these are determined according to a numeric stability criterion
  V_.P21ex_ = propagator_32( P_.tau_syn_, P_.tau_mem_, P_.tau_syn_, h );
  V_.P21in_ = propagator_32( P_.tau_syn_, P_.tau_mem_, P_.tau_syn_, h );

  // P21ex_ = h/C_;
  // P21in_ = h/C_;

  V_.P20_ = P_.tau_mem_ / ( 1.0 - V_.P22_ );
  // P20_ = h/C_;

  // TauR specifies the length of the absolute refractory period as
  // a double in ms. The grid based spike_response_0 can only handle refractory
  // periods that are integer multiples of the computation step size (h).
  // To ensure consistency with the overall simulation scheme such conversion
  // should be carried out via objects of class nest::Time. The conversion
  // requires 2 steps:
  //     1. A time object r is constructed defining  representation of
  //        TauR in tics. This representation is then converted to computation
  //        time steps again by a strategy defined by class nest::Time.
  //     2. The refractory time in units of steps is read out get_steps(), a
  //        member function of class nest::Time.
  //
  // Choosing a TauR that is not an integer multiple of the computation time
  // step h will leed to accurate (up to the resolution h) and self-consistent
  // results. However, a neuron model capable of operating with real valued
  // spike time may exhibit a different effective refractory time.

}

void
nest::spike_response_0::update( const Time& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  // evolve from timestep 'from' to timestep 'to' with steps of h each
  for ( long lag = from; lag < to; ++lag )
  {
    S_.V_m_ = S_.V_m_ * V_.P22_ + S_.i_syn_ex_ * V_.P21ex_
      + S_.i_syn_in_ * V_.P21in_ 
      + P_.xi_ * std::exp( (S_.t_last_spike_ - (origin.get_steps() + lag)) / P_.tau_mem_ ) * Time::get_resolution().get_ms();
    //TODO the last part is Forward Euler..
    //
    
    std::cout << S_.t_last_spike_ << " " << origin.get_steps() + lag << " " << std::exp( (S_.t_last_spike_ - (origin.get_steps() + lag)) / P_.tau_mem_ ) << std::endl;
    // exponential decaying PSCs
    S_.i_syn_ex_ *= V_.P11ex_;
    S_.i_syn_in_ *= V_.P11in_;

    // add evolution of presynaptic input current
    S_.i_syn_ex_ += ( 1. - V_.P11ex_ ) * S_.i_1_;

    // the spikes arriving at T+1 have an immediate effect on the state of the
    // neuron

    V_.weighted_spikes_ex_ = B_.spikes_ex_.get_value( lag );
    V_.weighted_spikes_in_ = B_.spikes_in_.get_value( lag );

    S_.i_syn_ex_ += V_.weighted_spikes_ex_;
    S_.i_syn_in_ += V_.weighted_spikes_in_;

    double rate = P_.rho0_ * std::exp( (S_.V_m_ - P_.theta_) / P_.du_ );

    // rate_ is in Hz, dt in ms, so we have to convert from s to ms
    V_.poisson_dev_.set_lambda(
      Time::get_resolution().get_ms() * rate * 1e-3 );

    librandom::RngPtr rng = kernel().rng_manager.get_rng( get_thread() );
    long n_spikes = V_.poisson_dev_.ldev( rng );

    if ( n_spikes > 0 ) // we must not send events with multiplicity 0
    {
       S_.t_last_spike_ = origin.get_steps() + lag;
       set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

       SpikeEvent se;
       se.set_multiplicity( n_spikes );
       kernel().event_delivery_manager.send( *this, se, lag );
    }

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
nest::spike_response_0::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );

  if ( e.get_weight() >= 0.0 )
  {
    B_.spikes_ex_.add_value( e.get_rel_delivery_steps(
                               kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
  else
  {
    B_.spikes_in_.add_value( e.get_rel_delivery_steps(
                               kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
}

void
nest::spike_response_0::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}
