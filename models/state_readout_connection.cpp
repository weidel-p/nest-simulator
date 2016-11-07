/*
 *  state_readout_connection.cpp
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

#include "state_readout_connection.h"

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connector_model.h"
#include "event.h"
#include "kernel_manager.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{
//
// Implementation of class StateReadoutCommonProperties.
//

StateReadoutCommonProperties::StateReadoutCommonProperties()
  : CommonSynapseProperties()
  , vt_( 0 )
  , Aplus_( 0.1 )
  , Aminus_( 0.1 )
  , tau_plus_( 20.0 )
  , tau_n_( 200.0 )
  , b_( 0.0 )
  , mean_firing_rate_( 5.0 )
  , n_lower_threshold_( 1.0 )
  , n_upper_threshold_( 2.0 )
  , Wmin_( 0.0 )
  , Wmax_( 200.0 )
{
}

void
StateReadoutCommonProperties::get_status( DictionaryDatum& d ) const
{
  CommonSynapseProperties::get_status( d );

  if ( vt_ != 0 )
    def< long >( d, "vt", vt_->get_gid() );
  else
    def< long >( d, "vt", -1 );

  def< double >( d, "Aplus", Aplus_);
  def< double >( d, "Aminus", Aminus_);
  def< double >( d, "tau_plus", tau_plus_ );
  def< double >( d, "tau_n", tau_n_ );
  def< double >( d, "b", b_ );
  def< double >( d, "mean_firing_rate", mean_firing_rate_ );
  def< double >( d, "n_lower_threshold", n_lower_threshold_ );
  def< double >( d, "n_upper_threshold", n_upper_threshold_ );
  def< double >( d, "Wmin", Wmin_ );
  def< double >( d, "Wmax", Wmax_ );
}

void
StateReadoutCommonProperties::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  CommonSynapseProperties::set_status( d, cm );

  long vtgid;
  if ( updateValue< long >( d, "vt", vtgid ) )
  {
    vt_ = dynamic_cast< volume_transmitter* >(
      kernel().node_manager.get_node( vtgid ) );

    if ( vt_ == 0 )
      throw BadProperty( "Dopamine source must be volume transmitter" );
  }

  updateValue< double >( d, "Aplus", Aplus_);
  updateValue< double >( d, "Aminus", Aminus_);
  updateValue< double >( d, "tau_plus", tau_plus_ );
  updateValue< double >( d, "tau_n", tau_n_ );
  updateValue< double >( d, "b", b_ );
  updateValue< double >( d, "mean_firing_rate", mean_firing_rate_ );
  updateValue< double >( d, "n_lower_threshold", n_lower_threshold_ );
  updateValue< double >( d, "n_upper_threshold", n_upper_threshold_ );
  updateValue< double >( d, "Wmin", Wmin_ );
  updateValue< double >( d, "Wmax", Wmax_ );
}

Node*
StateReadoutCommonProperties::get_node()
{
  if ( vt_ == 0 )
    throw BadProperty(
      "No volume transmitter has been assigned to the synapse." );
  else
    return vt_;
}




} // of namespace nest
