/*
 *  stdp_izh_connection.cpp
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

#include "stdp_izh_connection.h"

// C++ includes:
#include <cmath>
#include <limits>

namespace nest
{


STDPIzhCommonProperties::STDPIzhCommonProperties()
  : CommonSynapseProperties()
  , LTP_(0.1)
  , LTD_(-0.12)
  , tau_LTP_(20.0)
  , tau_LTD_(20.0)
  , Wmax_(10.0)
  , tau_syn_update_interval_(10000.)
  , constant_additive_value_(0.01)
{
}

void
STDPIzhCommonProperties::get_status( DictionaryDatum& d ) const
{
  CommonSynapseProperties::get_status( d );
  def< double >( d, "LTP", LTP_);
  def< double >( d, "LTD", LTD_);
  def< double >( d, "tau_LTP", tau_LTP_);
  def< double >( d, "tau_LTD", tau_LTD_);
  def< double >( d, "Wmax", Wmax_);
  def< double >( d, "tau_syn_update_interval", tau_syn_update_interval_ );
  def< double >( d, "constant_additive_value", constant_additive_value_ );
}

void
STDPIzhCommonProperties::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  CommonSynapseProperties::set_status( d, cm );
  updateValue< double >( d, "LTP", LTP_);
  updateValue< double >( d, "LTD", LTD_);
  updateValue< double >( d, "tau_LTP", tau_LTP_);
  updateValue< double >( d, "tau_LTD", tau_LTD_);
  updateValue< double >( d, "Wmax", Wmax_);
  updateValue< double >( d, "tau_syn_update_interval", tau_syn_update_interval_ );
  updateValue< double >( d, "constant_additive_value", constant_additive_value_ );
}


STDPIzhConnection::STDPIzhConnection()
  : ConnectionBase()
  , weight_(1.0)
  , wdev_(0.0)
  , t_last_update_(0.0)
  , t_last_post_spike_(0.0)
  , consistent_integration_(true)
{
  pre_spikes_.clear();
  //pre_spikes_.push_back(-std::numeric_limits<double>::max());
  pre_spikes_.push_back(-10000);

}

STDPIzhConnection::STDPIzhConnection( const STDPIzhConnection& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , wdev_( rhs.wdev_ )
  , t_last_update_( rhs.t_last_update_ )
  , t_last_post_spike_( rhs.t_last_post_spike_ )
  , consistent_integration_( rhs.consistent_integration_ )
{
  pre_spikes_.clear();
  pre_spikes_.push_back(-10000);

}

void
STDPIzhConnection::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, "wdev", wdev_);

  def< bool >( d, "consistent_integration", consistent_integration_);
}

void
STDPIzhConnection::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, "wdev", wdev_);

  updateValue< bool >( d, "consistent_integration", consistent_integration_);
}

void
STDPIzhConnection::time_driven_update( const thread tid, const double t_trig, const CommonPropertiesType& cp )
{

  Node* target = get_target( tid );
  const std::vector< double >& post_spikes_tmp = target->get_post_spikes();

  std::vector< double > post_spikes ( post_spikes_tmp.size() );
  std::copy(post_spikes_tmp.begin(), post_spikes_tmp.end(), post_spikes.begin());

  index i = 1; // index to iterate over post_spikes
  index j = 1; // index to iterate over pre_spikes_

  if (post_spikes[0] == t_last_update_){

      post_spikes.insert(post_spikes.begin(), t_last_post_spike_);
  }

  for ( j = 1; j < pre_spikes_.size() && pre_spikes_[j] <= t_trig; ++j )
  {

    while ( i < post_spikes.size() && post_spikes[i] < pre_spikes_[j] )
    {
      int dt = post_spikes[i] - pre_spikes_[j-1];

      if (post_spikes[i] == pre_spikes_[j])
      {
        dt = 0.;
      }

      wdev_ += cp.LTP_ * std::exp( -dt / cp.tau_LTP_ );
      ++i;
    }

    // facilitation (also for t_pre_spike == t_post_spike)
    // depression (also for t_pre_spike == t_post_spike)
    int dt = pre_spikes_[j] - post_spikes[i-1];
    wdev_ += cp.LTD_ * std::exp( -(dt-1) / cp.tau_LTD_ );
  }
  
  // process remaining postsynaptic spikes in this update interval if there are any
  while ( i < post_spikes.size() && post_spikes[i] < t_trig )
  {
    // facilitation
    int dt = post_spikes[i] - pre_spikes_[j-1];

    wdev_ += cp.LTP_ * std::exp( -dt / cp.tau_LTP_ );

    ++i;
  }
  
  wdev_ *= std::exp( - kernel().simulation_manager.get_syn_update_interval() / cp.tau_syn_update_interval_ );
  weight_ += cp.constant_additive_value_ * (kernel().simulation_manager.get_syn_update_interval() / 1000.) + wdev_;
  if (weight_ > cp.Wmax_)
     weight_ = cp.Wmax_;
  if (weight_ < 0)
     weight_ = 0; 

  // erase all processed presynaptic spikes except the last one
  // due to axonal there might be other pre_spikes left that are relevant only in the next update
  pre_spikes_.erase(pre_spikes_.begin(), pre_spikes_.begin() + (j-1) );

  t_last_update_ = t_trig;
  t_last_post_spike_ = post_spikes[post_spikes.size()-1];

}

} // of namespace nest

