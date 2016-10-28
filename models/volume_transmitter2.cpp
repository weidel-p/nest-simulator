/*
 *  volume_transmitter2.cpp
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

#include "volume_transmitter2.h"

// C++ includes:
#include <numeric>

// Includes from nestkernel:
#include "connector_base.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "spikecounter.h"

// Includes from sli:
#include "arraydatum.h"
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Default constructor defining default parameters
 * ---------------------------------------------------------------- */

nest::volume_transmitter2::Parameters_::Parameters_():
    tau_n_(20.)
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::volume_transmitter2::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, "tau_n", tau_n_);
}

void ::nest::volume_transmitter2::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, "tau_n", tau_n_);
}

/* ----------------------------------------------------------------
 * Default and copy constructor for volume transmitter2
 * ---------------------------------------------------------------- */

nest::volume_transmitter2::volume_transmitter2()
  : Archiving_Node()
  , P_()
{
}

nest::volume_transmitter2::volume_transmitter2( const volume_transmitter2& n )
  : Archiving_Node( n )
  , P_( n.P_ )
{
}

void
nest::volume_transmitter2::init_state_( const Node& )
{
}

void
nest::volume_transmitter2::init_buffers_()
{
  Archiving_Node::clear_history();
  B_.n_ = 0;
  B_.last_update_ = 0;
}

void
nest::volume_transmitter2::calibrate()
{

}

void
nest::volume_transmitter2::update( const Time& t, const long from, const long to )
{ 
  //std::cout << "update from " << from << " to " << to << " time " << t.get_ms() << std::endl;
    
  //for ( long lag = from; lag < to; ++lag )
  //for ( std::vector< SpikeEvent >::iterator e = B_.spikes_.begin();
  //      e != B_.spikes_.end();
  //      ++e )
  //{
  //  set_spiketime( t );
  //}
  //B_.spikes_.clear();
  //
  //B_.n_ = get_K_value(t.get_ms());
  //std::cout << " dopa " << B_.n_ << std::endl;
  //
    double minus_dt = B_.last_update_ - t.get_ms(); 
    B_.n_ *= std::exp( minus_dt / P_.tau_n_ );
    B_.last_update_ = t.get_ms();

    B_.history_.insert(std::pair<double, double>(t.get_ms(), B_.n_) );
}

void
nest::volume_transmitter2::handle( SpikeEvent& e )
{
    //for (int i = 0; i < e.get_multiplicity(); ++i){
    //    //B_.spikes_.push_back( *e.clone() );
    //    //set_spiketime( e.get_stamp() );
    //}
        
    double minus_dt = B_.last_update_ - e.get_stamp().get_ms(); 
    B_.n_ = B_.n_ * std::exp( minus_dt / P_.tau_n_ ) + e.get_multiplicity() / P_.tau_n_;
    B_.last_update_ = e.get_stamp().get_ms();
    B_.history_.insert(std::pair<double, double>(e.get_stamp().get_ms(), B_.n_) );
}

double
nest::volume_transmitter2::get_n(double t)
{
    //TODO is this accurate?
    return B_.history_.find(t)->second;
}

void
nest::volume_transmitter2::get_history( double t1,
  double t2,
  std::map<double, double>::iterator* start,
  std::map<double, double>::iterator* finish )
{
  *finish = B_.history_.end();
  if ( B_.history_.empty() )
  {
    *start = *finish;
    return;
  }
  else
  {
    std::map< double, double >::iterator runner = B_.history_.begin();
    while ( ( runner != B_.history_.end() ) && ( runner->first <= t1 ) )
      ++runner;
    *start = runner;
    while ( ( runner != B_.history_.end() ) && ( runner->first <= t2 ) )
    {
      ++runner;
    }
    *finish = runner;
  }
}

double
nest::volume_transmitter2::get_tau_n()
{
    return P_.tau_n_;
}
