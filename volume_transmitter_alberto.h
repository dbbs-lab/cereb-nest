/*
 *  volume_transmitter_alberto.h
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

#ifndef VOLUME_TRANSMITTER_ALBERTO_H
#define VOLUME_TRANSMITTER_ALBERTO_H

// Includes from nestkernel:
#include "archiving_node.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "spikecounter.h"

// Includes from sli:
#include "namedatum.h"


namespace mynest
{

/** @BeginDocumentation
Name: volume_transmitter_alberto - Node used in combination with neuromodulated synaptic
plasticity. It collects all spikes emitted by the population of neurons
connected to the volume transmitter and transmits the signal to a user-specific
subset of synapses.

Description:
The volume transmitter is used in combination with neuromodulated
synaptic plasticty, plasticity that depends not only on the activity
of the pre- and the postsynaptic neuron but also on a non-local
neuromodulatory third signal. It collects the spikes from all neurons
connected to the volume transmitter and delivers the spikes to a
user-specific subset of synapses.  It is assumed that the
neuromodulatory signal is a function of the spike times of all spikes
emitted by the population of neurons connected to the volume
transmitter.  The neuromodulatory dynamics is calculated in the
synapses itself. The volume transmitter interacts in a hybrid
structure with the neuromodulated synapses. In addition to the
delivery of the neuromodulatory spikes triggered by every pre-synaptic
spike, the neuromodulatory spike history is delivered in discrete time
intervals of a manifold of the minimal synaptic delay. In order to
insure the link between the neuromodulatory synapses and the volume
transmitter, the volume transmitter is passed as a parameter when a
neuromodulatory synapse is defined. The implementation is based on the
framework presented in [1].

Examples:
/volume_transmitter_alberto Create /vol Set
/iaf_psc_alpha Create /pre_neuron Set
/iaf_psc_alpha Create /post_neuron Set
/iaf_psc_alpha Create /neuromod_neuron Set
/stdp_dopamine_synapse  << /vt vol >>  SetDefaults
neuromod_neuron vol Connect
pre_neuron post_neuron /stdp_dopamine_synapse Connect

Parameters:
deliver_interval - time interval given in d_min time steps, in which
                   the volume signal is delivered from the volume
                   transmitter to the assigned synapses

References:
[1] Potjans W, Morrison A and Diesmann M (2010). Enabling functional
    neural circuit simulations with distributed computing of
    neuromodulated plasticity.
    Front. Comput. Neurosci. 4:141. doi:10.3389/fncom.2010.00141

Author: Wiebke Potjans, Abigail Morrison

Remarks: major changes to update function after code revision in Apr 2013 (SK)

Receives: SpikeEvent

SeeAlso: stdp_dopamine_synapse

*/
class ConnectorBase;

class volume_transmitter_alberto : public nest::Archiving_Node
{

public:
  volume_transmitter_alberto();
  volume_transmitter_alberto( const volume_transmitter_alberto& );

  bool
  has_proxies() const
  {
    return false;
  }
  bool
  local_receiver() const
  {
    return false;
  }

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using nest::Node::handle;
  using nest::Node::handles_test_event;

  void handle( nest::SpikeEvent& );

  nest::port handles_test_event( nest::SpikeEvent&, nest::rport );

  void get_status( DictionaryDatum& d ) const;
  void set_status( const DictionaryDatum& d );

  /**
   * Since volume transmitters are duplicated on each thread, and are
   * hence treated just as devices during node creation, we need to
   * define the corresponding setter and getter for local_device_id.
   **/
  void set_local_device_id( const nest::index ldid );
  nest::index get_local_device_id() const;

  const std::vector< nest::spikecounter >& deliver_spikes();

private:
  void init_state_( nest::Node const& );
  void init_buffers_();
  void calibrate();

  void update( const nest::Time&, const long, const long );

  // --------------------------------------------

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {
    Parameters_();
    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum& );
    long deliver_interval_; //!< update interval in d_min time steps
    long vt_num_;
  };

  //-----------------------------------------------

  struct Buffers_
  {
    nest::RingBuffer neuromodulatory_spikes_; //!< buffer to store incoming spikes
    //! vector to store and deliver spikes
    std::vector< nest::spikecounter > spikecounter_;
  };

  Parameters_ P_;
  Buffers_ B_;

  nest::index local_device_id_;
};

inline nest::port
volume_transmitter_alberto::handles_test_event( nest::SpikeEvent&, nest::rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw nest::UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline void
volume_transmitter_alberto::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  nest::Archiving_Node::get_status( d );

  ( *d )[ nest::names::element_type ] = LiteralDatum( nest::names::other );
}

inline void
volume_transmitter_alberto::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d );         // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  nest::Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
}

inline const std::vector< nest::spikecounter >&
volume_transmitter_alberto::deliver_spikes()
{
  return B_.spikecounter_;
}

inline void
volume_transmitter_alberto::set_local_device_id( const nest::index ldid )
{
  local_device_id_ = ldid;
}

inline nest::index
volume_transmitter_alberto::get_local_device_id() const
{
  return local_device_id_;
}

} // namespace

#endif /* #ifndef volume_transmitter_alberto_H */
