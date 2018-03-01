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

#include "stdp_izh_bitwise_correct_connection.h"

// C++ includes:
#include <cmath>
#include <limits>

namespace nest
{


STDPIzhBitwiseCorrectCommonProperties::STDPIzhBitwiseCorrectCommonProperties()
  : CommonSynapseProperties()
{
  pow_0_95_K_plus_ = new std::vector< double >(); 
  pow_0_95_K_minus_ = new std::vector< double >(); 

  double tmp_0_1 = 0.1;
  double tmp_0_12 = 0.12;
  for (int i = 0; i < 20001; ++i)
  {
    pow_0_95_K_plus_->push_back(tmp_0_1);
    pow_0_95_K_minus_->push_back(tmp_0_12);

    tmp_0_1 *= 0.95;
    tmp_0_12 *= 0.95;
  }
}


STDPIzhBitwiseCorrectConnection::STDPIzhBitwiseCorrectConnection()
  : ConnectionBase()
  , weight_(1.0)
  , wdev_(0.0)
  , K_plus_(0.0)
  , K_minus_(0.0)
  , tau_plus_(20.0)
  , tau_minus_(20.0)
  , lambda_(0.1)
  , alpha_(1.2)
  , Wmax_(10.0)
  , t_last_update_(0.0)
  , t_last_post_spike_(0.0)
  , consistent_integration_(true)
  , plot_(false)
{
  pre_spikes_.clear();
  //pre_spikes_.push_back(-std::numeric_limits<double>::max());
  pre_spikes_.push_back(-10000);

}

STDPIzhBitwiseCorrectConnection::STDPIzhBitwiseCorrectConnection( const STDPIzhBitwiseCorrectConnection& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , wdev_( rhs.wdev_ )
  , K_plus_( rhs.K_plus_ )
  , K_minus_( rhs. K_minus_ )
  , tau_plus_( rhs.tau_plus_ )
  , tau_minus_( rhs. tau_minus_ )
  , lambda_( rhs.lambda_ )
  , alpha_( rhs.alpha_ )
  , Wmax_( rhs.Wmax_ )
  , t_last_update_( rhs.t_last_update_ )
  , t_last_post_spike_( rhs.t_last_post_spike_ )
  , consistent_integration_( rhs.consistent_integration_ )
  , plot_( rhs.plot_)
{
  pre_spikes_.clear();
  pre_spikes_.push_back(-10000);

}

void
STDPIzhBitwiseCorrectConnection::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, "wdev", wdev_);
  def< double >( d, "K_plus", K_plus_);
  def< double >( d, "K_minus", K_minus_);
  def< double >( d, "tau_plus", tau_plus_);
  def< double >( d, "tau_minus", tau_minus_);
  def< double >( d, "lambda", lambda_);
  def< double >( d, "alpha", alpha_);
  def< double >( d, "Wmax", Wmax_);
  def< bool >( d, "consistent_integration", consistent_integration_);
  def< bool >( d, "plot", plot_);
}

void
STDPIzhBitwiseCorrectConnection::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, "wdev", wdev_);
  updateValue< double >( d, "tau_plus", tau_plus_);
  updateValue< double >( d, "tau_minus", tau_minus_);
  updateValue< double >( d, "lambda", lambda_);
  updateValue< double >( d, "alpha", alpha_);
  updateValue< double >( d, "Wmax", Wmax_);
  updateValue< bool >( d, "consistent_integration", consistent_integration_);
  updateValue< bool >( d, "plot", plot_);
}

void
STDPIzhBitwiseCorrectConnection::time_driven_update( const thread tid, const double t_trig, const CommonPropertiesType& cp )
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

      wdev_ += cp.pow_0_95_K_plus_->at(dt);
      
      if (plot_){
        std::cout << "LTP " << pre_spikes_[j-1] << " " << post_spikes[i] << " " << wdev_  << " " << cp.pow_0_95_K_plus_->at(dt) << std::endl;
      }

      ++i;
    }
    int dt = pre_spikes_[j] - post_spikes[i-1];

    wdev_ -= cp.pow_0_95_K_minus_->at(dt-1);
    
    if (plot_){
        std::cout << "LTD " << pre_spikes_[j] << " " << post_spikes[i-1] << " " << wdev_  << " " <<  cp.pow_0_95_K_minus_->at(dt-1) << std::endl;
    }

  }
  
  
  // process remaining postsynaptic spikes in this update interval if there are any
  while ( i < post_spikes.size() && post_spikes[i] < t_trig )
  {
    // facilitation
    int dt = post_spikes[i] - pre_spikes_[j-1];

    wdev_ += cp.pow_0_95_K_plus_->at(dt);

    if (plot_){
        std::cout << "LTP " << pre_spikes_[j-1] << " " << post_spikes[i] << " " << wdev_  << " " << cp.pow_0_95_K_plus_->at(dt) << std::endl;
    }
    ++i;
  }

  wdev_ *= 0.9;
  weight_ += 0.01 + wdev_;
  if (weight_ > Wmax_)
     weight_ = Wmax_;
  if (weight_ < 0)
     weight_ = 0; 


  // erase all processed presynaptic spikes except the last one
  // due to axonal there might be other pre_spikes left that are relevant only in the next update
  pre_spikes_.erase(pre_spikes_.begin(), pre_spikes_.begin() + (j-1) );


  t_last_update_ = t_trig;
  t_last_post_spike_ = post_spikes[post_spikes.size()-1];

}

} // of namespace nest

