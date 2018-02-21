/*
 *  global_volume_transmitter.cpp
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

#include "global_volume_transmitter.h"

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

nest::global_volume_transmitter::Parameters_::Parameters_()
  : deliver_interval_( 1 ) // in steps of mindelay
  , tau_(20.)
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::global_volume_transmitter::Parameters_::get( DictionaryDatum& d ) const
{
  def< long >( d, names::deliver_interval, deliver_interval_ );
  def< double >( d, "tau", tau_);
}

void ::nest::global_volume_transmitter::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< long >( d, names::deliver_interval, deliver_interval_ );
  updateValue< double >( d, "tau", tau_);
}

/* ----------------------------------------------------------------
 * Default and copy constructor for global_volume transmitter
 * ---------------------------------------------------------------- */

nest::global_volume_transmitter::global_volume_transmitter()
  : Archiving_Node()
  , P_()
{
}

nest::global_volume_transmitter::global_volume_transmitter( const global_volume_transmitter& n )
  : Archiving_Node( n )
  , P_( n.P_ )
{
}

void
nest::global_volume_transmitter::init_state_( const Node& )
{
}

void
nest::global_volume_transmitter::init_buffers_()
{
  B_.trace_ = 0;
  Archiving_Node::clear_history();
}

void
nest::global_volume_transmitter::calibrate()
{
}

void
nest::global_volume_transmitter::update( const Time&, const long from, const long to )
{
 
  B_.trace_ *= std::exp(-(to-from) / P_.tau_);

  // all spikes stored in spikecounter_ are delivered to the target synapses
  if ( ( kernel().simulation_manager.get_slice_origin().get_steps() + to )
      % ( P_.deliver_interval_ * kernel().connection_manager.get_min_delay() )
    == 0 )
  {
    double t_trig =
      Time(
        Time::step( kernel().simulation_manager.get_slice_origin().get_steps()
          + to ) ).get_ms();

      kernel().connection_manager.trigger_update_weight(
        get_gid(), B_.trace_, t_trig );

  }
}

void
nest::global_volume_transmitter::handle( SpikeEvent& e )
{
    B_.trace_ += e.get_multiplicity()/P_.tau_;
}
