/*
 *  dopa_connection.cpp
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

#include "dopa_connection.h"

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
// Implementation of class DopaCommonProperties.
//

DopaCommonProperties::DopaCommonProperties()
  : CommonSynapseProperties()
  , vt_( 0 )
  , A_( 0.1 )
  , tau_( 100.0 )
  , tau_n_( 20.0 )
  , b_plus_( 0.0 )
  , b_minus_( 0.0 )
  , n_threshold_( 2.0 )
  , Wmin_( 0.0 )
  , Wmax_( 200.0 )
  , LTD_scaling_( 1.0 )
  , tau_decay_( 0.0 )
  , weight0_( 1.0 )
  , heaviside_( true )
  , target_rate_ (8.)
  , scaling_ratio_(0.)
{
}

void
DopaCommonProperties::get_status( DictionaryDatum& d ) const
{
  CommonSynapseProperties::get_status( d );

  if ( vt_ != 0 )
    def< long >( d, names::vt, vt_->get_gid() );
  else
    def< long >( d, names::vt, -1 );

  def< double >( d, names::A, A_);
  def< double >( d, names::tau, tau_ );
  def< double >( d, names::tau_n, tau_n_ );
  def< double >( d, names::b_plus, b_plus_ );
  def< double >( d, names::b_minus, b_minus_ );
  def< double >( d, names::n_threshold, n_threshold_ );
  def< double >( d, names::Wmin, Wmin_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< double >( d, names::LTD_scaling, LTD_scaling_ );
  def< double >( d, names::tau_decay, tau_decay_);
  def< double >( d, names::weight0, weight0_);
  def< bool >( d, names::heaviside, heaviside_);
  def< double >( d, names::target_rate, target_rate_);
  def< double >( d, names::scaling_ratio, scaling_ratio_);
}

void
DopaCommonProperties::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  CommonSynapseProperties::set_status( d, cm );

  long vtgid;
  if ( updateValue< long >( d, names::vt, vtgid ) )
  {
    vt_ = dynamic_cast< global_volume_transmitter* >(
      kernel().node_manager.get_node( vtgid ) );

    if ( vt_ == 0 )
      throw BadProperty( "Dopamine source must be volume transmitter" );
  }


  updateValue< double >( d, names::A, A_);
  updateValue< double >( d, names::tau, tau_);
  updateValue< double >( d, names::tau_n, tau_n_ );
  updateValue< double >( d, names::n_threshold, n_threshold_ );
  updateValue< double >( d, names::b_plus, b_plus_ );
  updateValue< double >( d, names::b_minus, b_minus_ );
  updateValue< double >( d, names::Wmin, Wmin_ );
  updateValue< double >( d, names::Wmax, Wmax_ );
  updateValue< double >( d, names::LTD_scaling, LTD_scaling_);
  updateValue< double >( d, names::tau_decay, tau_decay_);
  updateValue< double >( d, names::weight0, weight0_);
  updateValue< bool >( d, names::heaviside, heaviside_);
  updateValue< double >( d, names::target_rate, target_rate_);
  updateValue< double >( d, names::scaling_ratio, scaling_ratio_);
}

Node*
DopaCommonProperties::get_node()
{
  if ( vt_ == 0 )
    throw BadProperty(
      "No volume transmitter has been assigned to the synapse." );
  else
  {
    return vt_;
  }

}




} // of namespace nest
