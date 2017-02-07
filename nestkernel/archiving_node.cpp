/*
 *  archiving_node.cpp
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

/**
 * \file archiving_node.cpp
 * Implementation of archiving_node to record and manage spike history
 * \author Moritz Helias, Abigail Morrison
 * \date april 2006
 */

#include "archiving_node.h"

// Includes from sli:
#include "dictutils.h"

namespace nest
{

// member functions for Archiving_Node

nest::Archiving_Node::Archiving_Node()
  : n_incoming_( 0 )
  , Kminus_( 0.0 )
  , triplet_Kminus_( 0.0 )
  , tau_minus_( 20.0 )
  , tau_minus_triplet_( 110.0 )
  , last_spike_( -1.0 )
  , firing_rate_short_( 0.0 )
  , firing_rate_long_( 0.0 )
  , tau_rate_short_( 100.0 )
  , tau_rate_long_( 10000.0 )
  , Ca_t_( 0.0 )
  , Ca_minus_( 0.0 )
  , tau_Ca_( 10000.0 )
  , beta_Ca_( 0.001 )
  , synaptic_elements_map_()
{
}

nest::Archiving_Node::Archiving_Node( const Archiving_Node& n )
  : Node( n )
  , n_incoming_( n.n_incoming_ )
  , Kminus_( n.Kminus_ )
  , triplet_Kminus_( n.triplet_Kminus_ )
  , tau_minus_( n.tau_minus_ )
  , tau_minus_triplet_( n.tau_minus_triplet_ )
  , last_spike_( n.last_spike_ )
  , firing_rate_short_( n.firing_rate_short_)
  , firing_rate_long_( n.firing_rate_long_)
  , tau_rate_short_( n.tau_rate_short_)
  , tau_rate_long_( n.tau_rate_long_)
  , Ca_t_( n.Ca_t_ )
  , Ca_minus_( n.Ca_minus_ )
  , tau_Ca_( n.tau_Ca_ )
  , beta_Ca_( n.beta_Ca_ )
  , synaptic_elements_map_( n.synaptic_elements_map_ )
{
}

void
Archiving_Node::register_stdp_connection( double t_first_read )
{
  // Mark all entries in the deque, which we will not read in future as read by
  // this input input, so that we savely increment the incoming number of
  // connections afterwards without leaving spikes in the history.
  // For details see bug #218. MH 08-04-22

  for ( std::deque< histentry >::iterator runner = history_.begin();
        runner != history_.end() && runner->t_ <= t_first_read;
        ++runner )
  {
    ( runner->access_counter_ )++;
  }

  n_incoming_++;
}

double
nest::Archiving_Node::get_K_value( double t )
{
  if ( history_.empty() )
  {
    return Kminus_;
  }
  int i = history_.size() - 1;
  while ( i >= 0 )
  {
    if ( t > history_[ i ].t_ )
    {
      return ( history_[ i ].Kminus_
        * std::exp( ( history_[ i ].t_ - t ) / tau_minus_ ) );
    }
    i--;
  }
  return 0;
}

double
nest::Archiving_Node::get_firing_rate_short( double t )
{
  if ( history_.empty() )
  {
    return firing_rate_short_;
  }
  int i = history_.size() - 1;
  while ( i >= 0 )
  {
    if ( t > history_[ i ].t_ )
    {
      return ( history_[ i ].firing_rate_short_
        * std::exp( ( history_[ i ].t_ - t ) / tau_rate_short_) );
    }
    i--;
  }
  return 0;
}


double
nest::Archiving_Node::get_firing_rate_long( double t )
{
  if ( history_.empty() )
  {
    return firing_rate_long_;
  }
  int i = history_.size() - 1;
  while ( i >= 0 )
  {
    if ( t > history_[ i ].t_ )
    {
      return ( history_[ i ].firing_rate_long_
        * std::exp( ( history_[ i ].t_ - t ) / tau_rate_long_) );
    }
    i--;
  }
  return 0;
}

void
nest::Archiving_Node::get_K_values( double t,
  double& K_value,
  double& triplet_K_value )
{
  // case when the neuron has not yet spiked
  if ( history_.empty() )
  {
    triplet_K_value = triplet_Kminus_;
    K_value = Kminus_;
    return;
  }
  // case
  int i = history_.size() - 1;
  while ( i >= 0 )
  {
    if ( t > history_[ i ].t_ )
    {
      triplet_K_value = ( history_[ i ].triplet_Kminus_
        * std::exp( ( history_[ i ].t_ - t ) / tau_minus_triplet_ ) );
      K_value = ( history_[ i ].Kminus_
        * std::exp( ( history_[ i ].t_ - t ) / tau_minus_ ) );
      return;
    }
    i--;
  }

  // we only get here if t< time of all spikes in history)

  // return 0.0 for both K values
  triplet_K_value = 0.0;
  K_value = 0.0;
}


double
nest::Archiving_Node::get_tau_minus()
{
   return tau_minus_;
}

void
nest::Archiving_Node::get_history( double t1,
  double t2,
  std::deque< histentry >::iterator* start,
  std::deque< histentry >::iterator* finish )
{
  *finish = history_.end();
  if ( history_.empty() )
  {
    *start = *finish;
    return;
  }
  else
  {
    std::deque< histentry >::reverse_iterator runner = history_.rbegin();
    while ( ( runner != history_.rend() ) && ( runner->t_ > t1 ) )
      ++runner;
    *start = runner.base();

    std::deque< histentry >::iterator forward_runner = runner.base(); 
    while ( ( forward_runner != history_.end() ) && ( forward_runner->t_ <= t2 ) )
    {
      ( forward_runner->access_counter_ )++;
      ++forward_runner;
    }
    *finish = forward_runner;
  }
}

void
nest::Archiving_Node::set_spiketime( Time const& t_sp, double offset )
{
  const double t_sp_ms = t_sp.get_ms() - offset;
  update_synaptic_elements( t_sp_ms );
  Ca_minus_ += beta_Ca_;

  if ( n_incoming_ )
  {
    // prune all spikes from history which are no longer needed
    // except the penultimate one. we might still need it.
    while ( history_.size() > 1 )
    {
      if ( history_.front().access_counter_ >= n_incoming_ )
      {
        history_.pop_front();
      }
      else
        break;
    }
    // update spiking history
    Kminus_ =
      Kminus_ * std::exp( ( last_spike_ - t_sp_ms ) / tau_minus_ ) + 1.0;
    triplet_Kminus_ = triplet_Kminus_
        * std::exp( ( last_spike_ - t_sp_ms ) / tau_minus_triplet_ )
      + 1.0;
    

    firing_rate_short_ =
      firing_rate_short_ * std::exp( ( last_spike_ - t_sp_ms ) / tau_rate_short_ ) + 1.0 / tau_rate_short_;

    firing_rate_long_ =
      firing_rate_long_ * std::exp( ( last_spike_ - t_sp_ms ) / tau_rate_long_) + 1.0 / tau_rate_long_;

    last_spike_ = t_sp_ms;
    history_.push_back( histentry( last_spike_, Kminus_, triplet_Kminus_, 0, firing_rate_short_, firing_rate_long_ ) );
  }
  else
  {
    last_spike_ = t_sp_ms;
  }
}

void
nest::Archiving_Node::get_status( DictionaryDatum& d ) const
{
  DictionaryDatum synaptic_elements_d;
  DictionaryDatum synaptic_element_d;

  def< double >( d, names::t_spike, get_spiketime_ms() );
  def< double >( d, names::tau_minus, tau_minus_ );
  def< double >( d, "tau_rate_short", tau_rate_short_ );
  def< double >( d, "tau_rate_long", tau_rate_long_ );
  def< double >( d, "firing_rate_short", firing_rate_short_ );
  def< double >( d, "firing_rate_long", firing_rate_long_ );
  def< double >( d, names::Ca, Ca_minus_ );
  def< double >( d, names::tau_Ca, tau_Ca_ );
  def< double >( d, names::beta_Ca, beta_Ca_ );
  def< double >( d, names::tau_minus_triplet, tau_minus_triplet_ );
#ifdef DEBUG_ARCHIVER
  def< int >( d, names::archiver_length, history_.size() );
#endif

  synaptic_elements_d = DictionaryDatum( new Dictionary );
  def< DictionaryDatum >( d, names::synaptic_elements, synaptic_elements_d );
  for ( std::map< Name, SynapticElement >::const_iterator it =
          synaptic_elements_map_.begin();
        it != synaptic_elements_map_.end();
        ++it )
  {
    synaptic_element_d = DictionaryDatum( new Dictionary );
    def< DictionaryDatum >(
      synaptic_elements_d, it->first, synaptic_element_d );
    it->second.get( synaptic_element_d );
  }
}

void
nest::Archiving_Node::set_status( const DictionaryDatum& d )
{
  // We need to preserve values in case invalid values are set
  double new_tau_minus = tau_minus_;
  double new_tau_minus_triplet = tau_minus_triplet_;
  double new_tau_rate_short = tau_rate_short_;
  double new_tau_rate_long = tau_rate_long_;
  double new_firing_rate_short = firing_rate_short_;
  double new_firing_rate_long = firing_rate_long_;
  double new_tau_Ca = tau_Ca_;
  double new_beta_Ca = beta_Ca_;
  updateValue< double >( d, names::tau_minus, new_tau_minus );
  updateValue< double >( d, names::tau_minus_triplet, new_tau_minus_triplet );
  updateValue< double >( d, names::tau_Ca, new_tau_Ca );
  updateValue< double >( d, names::beta_Ca, new_beta_Ca );
  updateValue< double >( d, "tau_rate_short", new_tau_rate_short );
  updateValue< double >( d, "tau_rate_long", new_tau_rate_long );
  updateValue< double >( d, "firing_rate_short", new_firing_rate_short);
  updateValue< double >( d, "firing_rate_long", new_firing_rate_long);

  if ( new_tau_minus <= 0.0 || new_tau_minus_triplet <= 0.0 )
  {
    throw BadProperty( "All time constants must be strictly positive." );
  }

  tau_minus_ = new_tau_minus;
  tau_minus_triplet_ = new_tau_minus_triplet;
  tau_rate_short_ = new_tau_rate_short;
  tau_rate_long_ = new_tau_rate_long;
  firing_rate_short_ = new_firing_rate_short;
  firing_rate_long_ = new_firing_rate_long;

  if ( new_tau_Ca <= 0.0 )
  {
    throw BadProperty( "All time constants must be strictly positive." );
  }
  tau_Ca_ = new_tau_Ca;

  if ( new_beta_Ca <= 0.0 )
  {
    throw BadProperty(
      "For Ca to function as an integrator of the electrical activity, beta_ca "
      "needs to be greater than 0." );
  }
  beta_Ca_ = new_beta_Ca;

  // check, if to clear spike history and K_minus
  bool clear = false;
  updateValue< bool >( d, names::clear, clear );
  if ( clear )
  {
    clear_history();
  }

  if ( not d->known( names::synaptic_elements ) )
  {
    return;
  }
  // we replace the existing synaptic_elements_map_ by the new one
  DictionaryDatum synaptic_elements_d;
  std::pair< std::map< Name, SynapticElement >::iterator, bool > insert_result;

  synaptic_elements_map_ = std::map< Name, SynapticElement >();
  synaptic_elements_d =
    getValue< DictionaryDatum >( d, names::synaptic_elements );

  for ( Dictionary::const_iterator i = synaptic_elements_d->begin();
        i != synaptic_elements_d->end();
        ++i )
  {
    insert_result = synaptic_elements_map_.insert(
      std::pair< Name, SynapticElement >( i->first, SynapticElement() ) );
    ( insert_result.first->second )
      .set( getValue< DictionaryDatum >( synaptic_elements_d, i->first ) );
  }
}

void
nest::Archiving_Node::clear_history()
{
  last_spike_ = -1.0;
  Kminus_ = 0.0;
  triplet_Kminus_ = 0.0;
  history_.clear();
  Ca_minus_ = 0.0;
  Ca_t_ = 0.0;
}


/* ----------------------------------------------------------------
* Get the number of synaptic_elements
* ---------------------------------------------------------------- */
double
nest::Archiving_Node::get_synaptic_elements( Name n ) const
{
  std::map< Name, SynapticElement >::const_iterator se_it;
  se_it = synaptic_elements_map_.find( n );
  double z_value;

  if ( se_it != synaptic_elements_map_.end() )
  {
    z_value = ( se_it->second ).get_z();
    if ( ( se_it->second ).continuous() )
    {
      return z_value;
    }
    else
    {
      return std::floor( z_value );
    }
  }
  else
  {
    return 0.0;
  }
}

int
nest::Archiving_Node::get_synaptic_elements_vacant( Name n ) const
{
  std::map< Name, SynapticElement >::const_iterator se_it;
  se_it = synaptic_elements_map_.find( n );

  if ( se_it != synaptic_elements_map_.end() )
  {
    return se_it->second.get_z_vacant();
  }
  else
  {
    return 0;
  }
}

int
nest::Archiving_Node::get_synaptic_elements_connected( Name n ) const
{
  std::map< Name, SynapticElement >::const_iterator se_it;
  se_it = synaptic_elements_map_.find( n );

  if ( se_it != synaptic_elements_map_.end() )
  {
    return se_it->second.get_z_connected();
  }
  else
  {
    return 0;
  }
}

std::map< Name, double >
nest::Archiving_Node::get_synaptic_elements() const
{
  std::map< Name, double > n_map;

  for ( std::map< Name, SynapticElement >::const_iterator it =
          synaptic_elements_map_.begin();
        it != synaptic_elements_map_.end();
        ++it )
  {
    n_map.insert( std::pair< Name, double >(
      it->first, get_synaptic_elements( it->first ) ) );
  }
  return n_map;
}

void
nest::Archiving_Node::update_synaptic_elements( double t )
{
  assert( t >= Ca_t_ );

  for ( std::map< Name, SynapticElement >::iterator it =
          synaptic_elements_map_.begin();
        it != synaptic_elements_map_.end();
        ++it )
  {
    it->second.update( t, Ca_t_, Ca_minus_, tau_Ca_ );
  }
  // Update calcium concentration
  Ca_minus_ = Ca_minus_ * std::exp( ( Ca_t_ - t ) / tau_Ca_ );
  Ca_t_ = t;
}

void
nest::Archiving_Node::decay_synaptic_elements_vacant()
{
  for ( std::map< Name, SynapticElement >::iterator it =
          synaptic_elements_map_.begin();
        it != synaptic_elements_map_.end();
        ++it )
  {
    it->second.decay_z_vacant();
  }
}

void
nest::Archiving_Node::connect_synaptic_element( Name name, int n )
{
  std::map< Name, SynapticElement >::iterator se_it;
  se_it = synaptic_elements_map_.find( name );

  if ( se_it != synaptic_elements_map_.end() )
  {
    se_it->second.connect( n );
  }
}

} // of namespace nest
