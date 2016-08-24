/*
 *  state_separation.cpp
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

#include "state_separation_connection.h"

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
// Implementation of class StateSeparationCommonProperties.
//

StateSeparationCommonProperties::StateSeparationCommonProperties()
  : CommonSynapseProperties()
  , vt_( 0 )
  , A_( 1.0 )
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
StateSeparationCommonProperties::get_status( DictionaryDatum& d ) const
{
  CommonSynapseProperties::get_status( d );

  if ( vt_ != 0 )
    def< long_t >( d, "vt", vt_->get_gid() );
  else
    def< long_t >( d, "vt", -1 );

  def< double_t >( d, "A", A_);
  def< double_t >( d, "tau_plus", tau_plus_ );
  def< double_t >( d, "tau_n", tau_n_ );
  def< double_t >( d, "b", b_ );
  def< double_t >( d, "mean_firing_rate", mean_firing_rate_ );
  def< double_t >( d, "n_lower_threshold", n_lower_threshold_ );
  def< double_t >( d, "n_upper_threshold", n_upper_threshold_ );
  def< double_t >( d, "Wmin", Wmin_ );
  def< double_t >( d, "Wmax", Wmax_ );
}

void
StateSeparationCommonProperties::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  CommonSynapseProperties::set_status( d, cm );

  long_t vtgid;
  if ( updateValue< long_t >( d, "vt", vtgid ) )
  {
    vt_ = dynamic_cast< volume_transmitter* >(
      kernel().node_manager.get_node( vtgid ) );

    if ( vt_ == 0 )
      throw BadProperty( "Dopamine source must be volume transmitter" );
  }

  updateValue< double_t >( d, "A", A_);
  updateValue< double_t >( d, "tau_plus", tau_plus_ );
  updateValue< double_t >( d, "tau_n", tau_n_ );
  updateValue< double_t >( d, "b", b_ );
  updateValue< double_t >( d, "mean_firing_rate", mean_firing_rate_ );
  updateValue< double_t >( d, "n_lower_threshold", n_lower_threshold_ );
  updateValue< double_t >( d, "n_upper_threshold", n_upper_threshold_ );
  updateValue< double_t >( d, "Wmin", Wmin_ );
  updateValue< double_t >( d, "Wmax", Wmax_ );
}

Node*
StateSeparationCommonProperties::get_node()
{
  if ( vt_ == 0 )
    throw BadProperty(
      "No volume transmitter has been assigned to the synapse." );
  else
    return vt_;
}




} // of namespace nest
