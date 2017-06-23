/*
 *  stdp_izh_naive_connection.cpp
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

#include "stdp_izh_naive_connection.h"

// C++ includes:
#include <cmath>
#include <limits>

namespace nest
{

STDPIzhNaiveConnection::STDPIzhNaiveConnection()
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
  , consistent_integration_(true)
{
  pre_spikes_.clear();
  pre_spikes_.push_back(-std::numeric_limits<double>::max());
}

STDPIzhNaiveConnection::STDPIzhNaiveConnection( const STDPIzhNaiveConnection& rhs )
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
  , consistent_integration_( rhs.consistent_integration_ )
{
  pre_spikes_.clear();
  pre_spikes_.push_back(-std::numeric_limits<double>::max());
}

void
STDPIzhNaiveConnection::get_status( DictionaryDatum& d ) const
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
}

void
STDPIzhNaiveConnection::set_status( const DictionaryDatum& d, ConnectorModel& cm )
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
}

void
STDPIzhNaiveConnection::time_driven_update( const thread tid, const double t_trig, const CommonPropertiesType& )
{
  //std::cout << "before update weight = " << weight_ << std::endl;

  Node* target = get_target( tid );
  const std::vector< double >& post_spikes = target->get_post_spikes();

  /* std::cout << "pre_spikes = [ "; */
  /* for ( std::vector< double >::const_iterator it = pre_spikes_.begin(); it < pre_spikes_.end(); ++it ) */
  /* { */
  /*   std::cout << (*it) << ", "; */
  /* } */
  /* std::cout << "]" << std::endl; */

  /* std::cout << "post_spikes = [ "; */
  /* for ( std::vector< double >::const_iterator it = post_spikes.begin(); it < post_spikes.end(); ++it ) */
  /* { */
  /*   std::cout << (*it) << ", "; */
  /* } */
  /* std::cout << "]" << std::endl; */

  index i = 1; // index to iterate over post_spikes
  index j = 1; // index to iterate over pre_spikes_

  // pre_spikes[0] and post_spikes[0] are the times of the last
  // presynaptic and postsynaptic spike in the last update interval, respectively

  for ( j = 1; j < pre_spikes_.size() && pre_spikes_[j] <= t_trig; ++j )
  {
    while ( i < post_spikes.size() && post_spikes[i] <= pre_spikes_[j] )
    {
      // facilitation (also for t_pre_spike == t_post_spike)
      wdev_ += lambda_ * K_plus_ * std::exp( ( pre_spikes_[j-1] - post_spikes[i] ) / tau_plus_ );
      //std::cout << "facilitation t_last_pre = " << pre_spikes_[j-1] << ", t_post = " << post_spikes[i] << ", wdev = " << wdev_ << std::endl;
      K_minus_ = K_minus_ * std::exp( ( post_spikes[i-1] - post_spikes[i] ) / tau_minus_ ) + 1.0;
      ++i;
    }
    
    // depression (also for t_pre_spike == t_post_spike)
    wdev_ -= alpha_ * lambda_ * K_minus_ * std::exp( ( post_spikes[i-1] - pre_spikes_[j] )  / tau_minus_);
    //std::cout << "depression t_last_post = " << post_spikes[i-1] << ", t_pre = " << pre_spikes_[j] << ", wdev = " << wdev_ << std::endl;
    K_plus_ = K_plus_ * std::exp( ( pre_spikes_[j-1] - pre_spikes_[j] ) / tau_plus_ ) + 1.0;
  }
  
  //std::cout << "no more pre_spikes in update interval" << std::endl;
  
  // process remaining postsynaptic spikes in this update interval if there are any
  while ( i < post_spikes.size() && post_spikes[i] <= t_trig )
  {
    // facilitation
    wdev_ += lambda_ * K_plus_ * std::exp( ( pre_spikes_[j-1] - post_spikes[i] ) / tau_plus_ );
    //std::cout << "facilitation t_last_pre = " << pre_spikes_[j-1] << ", t_post = " << post_spikes[i] << ", wdev = " << wdev_ << std::endl;
    K_minus_ = K_minus_ * std::exp( ( post_spikes[i-1] - post_spikes[i] ) / tau_minus_ ) + 1.0;
    ++i;
  }
  
  double w_new = weight_ + wdev_;
  if ( consistent_integration_ )
  {
    wdev_ = 0.0;
  }
  else
  {
    w_new += 0.01;
    wdev_ *= 0.9;
  }

  //std::cout << "before boundary check weight = " << w_new << std::endl;
  if ( w_new > 0.0 )
  {
    if ( w_new < Wmax_ )
    {
      weight_ = w_new;
    }
    else
    {
      weight_ = Wmax_;
    }
  }
  else
  {
    weight_ = 0.0;
  }
  //std::cout << "after update weight = " << std::setprecision(15) << weight_ << std::endl;

  // erase all processed presynaptic spikes except the last one
  // due to axonal there might be other pre_spikes left that are relevant only in the next update
  pre_spikes_.erase(pre_spikes_.begin(), pre_spikes_.begin() + (j-1) );

  /* std::cout << "after update pre_spikes = [ "; */
  /* for ( std::vector< double >::const_iterator it = pre_spikes_.begin(); it < pre_spikes_.end(); ++it ) */
  /* { */
  /*   std::cout << (*it) << ", "; */
  /* } */
  /* std::cout << "]" << std::endl; */

  t_last_update_ = t_trig;
}

} // of namespace nest

