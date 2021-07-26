/* generated by template org.nest.nestml.neuron.NeuronClass*/
/*

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

// C++ includes:
#include <limits>
#include <math.h>

// Includes from libnestutil:
#include "numerics.h"
#include "stdio.h"
#include <iostream>

#include <time.h>

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

#include "eglif_cond_alpha_multisyn.h"

std::vector <float> spikes_times_;

/* ----------------------------------------------------------------
* Recordables map
* ---------------------------------------------------------------- */
nest::RecordablesMap<eglif_cond_alpha_multisyn>
    eglif_cond_alpha_multisyn::recordablesMap_;

namespace nest {
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <> void RecordablesMap<eglif_cond_alpha_multisyn>::create() {
  // use standard names whereever you can for consistency!
  /* generated by template org.nest.nestml.function.RecordCallback*/


  insert_("I_dep", &eglif_cond_alpha_multisyn::get_I_dep);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("I_adap", &eglif_cond_alpha_multisyn::get_I_adap);


  insert_("V_th", &eglif_cond_alpha_multisyn::get_V_th);


  //insert_("y_0", &eglif_cond_alpha_multisyn::get_y_0);

  /* generated by template org.nest.nestml.function.RecordCallback*/

  insert_("V_m", &eglif_cond_alpha_multisyn::get_V_m);


  insert_("I_gen", &eglif_cond_alpha_multisyn::get_I_gen);

  insert_("I_input", &eglif_cond_alpha_multisyn::get_I_input);


  insert_("sum_buffer", &eglif_cond_alpha_multisyn::get_sum_buffer);


insert_("init_sign", &eglif_cond_alpha_multisyn::get_init_sign);


  insert_("G1", &eglif_cond_alpha_multisyn::get_G1);


  insert_("DG1", &eglif_cond_alpha_multisyn::get_DG1);


  insert_("G2", &eglif_cond_alpha_multisyn::get_G2);


  insert_("DG2", &eglif_cond_alpha_multisyn::get_DG2);


  insert_("G3", &eglif_cond_alpha_multisyn::get_G3);


  insert_("DG3", &eglif_cond_alpha_multisyn::get_DG3);


  insert_("G4", &eglif_cond_alpha_multisyn::get_G4);


  insert_("DG4", &eglif_cond_alpha_multisyn::get_DG4);

  /* generated by template org.nest.nestml.function.RecordCallback*/



  /* generated by template org.nest.nestml.function.RecordCallback*/

  //insert_("time", &eglif_cond_alpha_multisyn::get_time);

  /* generated by template org.nest.nestml.function.RecordCallback*/



  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the RefractoryCounts with the domain type long

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the r with the domain type long

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the bufferT with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tot_cur with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the receptors with the domain type long

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the C_m with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_m with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the tau_syn with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the t_ref with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the E_L with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the V_reset with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the a with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the b with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the V_th with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the k1 with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the k2 with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_e with the domain type double

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the bT_size with the domain type long

  /* generated by template org.nest.nestml.function.RecordCallback*/

  // ignores the I_tot with the domain type double
}
}

/* ----------------------------------------------------------------
* Default constructors defining default parameters and state
* ---------------------------------------------------------------- */

eglif_cond_alpha_multisyn::Parameters_::Parameters_() {}

eglif_cond_alpha_multisyn::State_::State_() {}

/* ----------------------------------------------------------------
* Parameter and state extractions and manipulation functions
* ---------------------------------------------------------------- */

void eglif_cond_alpha_multisyn::Parameters_::set(
    const DictionaryDatum &__d) {}

void eglif_cond_alpha_multisyn::State_::set(const DictionaryDatum &__d,
                                                   const Parameters_ &p) {}

eglif_cond_alpha_multisyn::Buffers_::Buffers_(
    eglif_cond_alpha_multisyn &n)
    : logger_(n), s_(0), c_(0), e_(0) {}

eglif_cond_alpha_multisyn::Buffers_::Buffers_(
    const Buffers_ &, eglif_cond_alpha_multisyn &n)
    : logger_(n), s_(0), c_(0), e_(0) {}

/* ----------------------------------------------------------------
* Default and copy constructor for node
* ---------------------------------------------------------------- */
// TODO inner components
eglif_cond_alpha_multisyn::eglif_cond_alpha_multisyn()
    : Archiving_Node(), P_(), S_(), B_(*this) {
  recordablesMap_.create();


  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.C_m = 189.79;

 /* generated by template org.nest.nestml.function.MemberInitialization*/


  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.tau_m = 2975.410306906496;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

 //.resize(V_.receptors, 2);

/* generated by template org.nest.nestml.function.MemberInitialization*/
 //P_.E_rev.resize(V_.receptors, 0.0);

  P_.tau_syn1 = 2.0;
  P_.tau_syn2 = 2.0;
  P_.tau_syn3 = 2.0;
  P_.tau_syn4 = 2.0;

  P_.E_rev1 = 0.0;
  P_.E_rev2 = 0.0;
  P_.E_rev3 = 0.0;
  P_.E_rev4 = 0.0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.t_ref = 2;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.E_L = -65.88080929487171;



  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.V_reset = get_E_L();


  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.V_th = ((-45.0));

  P_.Vinit = ((-60.0));

  P_.Vmin = ((-110.0));



  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.k1 = pow((200), (((-1))));

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.k2 = pow((20), (((-1))));


  /* generated by template org.nest.nestml.function.MemberInitialization*/

  // Parameters for exp cum update of Iadap
  P_.a = 117.91230509072211;

  P_.b = 0.5993799602025365;

  P_.c = -105.1597825342601 ;

  P_.ts = std::numeric_limits<int>::max();          // time threshold for depolarization block

  P_.mul = 1.0;

  P_.kadap = 0.0;

  // Additional parameters for enriched update functions
  P_.ith = 0.0;

  P_.sc = 27.678206947038134;

  P_.bet = 0.3279210955457432;

  P_.delta1 = 0.1518130616079973;

  P_.Idep_ini_vr = 0.516088363498455;

  P_.psi1 = 0.12242935394673728;

  P_.alp = 1.0219843023331783;

  P_.istim_min_spiking_exp = 200;


  /* generated by template org.nest.nestml.function.MemberInitialization*/

  P_.I_e = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/


  P_.bT_size = 200;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::I_input] = 0.0;

  S_.y_[State_::I_dep] = 0.0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::I_adap] = 0.0;



 // S_.y_[State_::y_0] = 0.0;

 // S_.y_[State_::y_1] = 0.0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  //S_.y_[State_::V_m] = ((-70.0));
 S_.y_[State_::V_m] = P_.V_reset;
 // S_.y_[State_::V_m] = P_.E_L;


  /* generated by template org.nest.nestml.function.MemberInitialization*/

  //S_.y_[State_::time] = 0;
  S_.time = 0;


/* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::G1] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::DG1] = 0;

/* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::G2] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::DG2] = 0;

/* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::G3] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::DG3] = 0;

/* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::G4] = 0;

  /* generated by template org.nest.nestml.function.MemberInitialization*/

  S_.y_[State_::DG4] = 0;

}

eglif_cond_alpha_multisyn::eglif_cond_alpha_multisyn(
    const eglif_cond_alpha_multisyn &n)
    : Archiving_Node(), P_(n.P_), S_(n.S_), B_(n.B_, *this) {}

/* ----------------------------------------------------------------
* Destructors
* ---------------------------------------------------------------- */

eglif_cond_alpha_multisyn::~eglif_cond_alpha_multisyn() {
  // GSL structs may not have been allocated, so we need to protect destruction
  if (B_.s_)
    gsl_odeiv_step_free(B_.s_);
  if (B_.c_)
    gsl_odeiv_control_free(B_.c_);
  if (B_.e_)
    gsl_odeiv_evolve_free(B_.e_);

}

/* ----------------------------------------------------------------
* Node initialization functions
* ---------------------------------------------------------------- */

void eglif_cond_alpha_multisyn::init_state_(
    const Node &proto) { // TODO inner components

  const eglif_cond_alpha_multisyn &pr =
      downcast<eglif_cond_alpha_multisyn>(proto);
  S_ = pr.S_;
}

/* generated by template org.nest.nestml.function.GSLDifferentiationFunction*/
extern "C" inline int
eglif_cond_alpha_multisyn_dynamics(double, const double y[], double f[],
                                          void *pnode) {

  //const double tIe = nest::kernel().simulation_manager.get_time().get_ms();

  typedef eglif_cond_alpha_multisyn::State_ State_;
  // get access to node so we can almost work as in a member function
  assert(pnode);
  const eglif_cond_alpha_multisyn &node =
      *(reinterpret_cast<eglif_cond_alpha_multisyn *>(pnode));

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y_[].


  // Synaptic current: I_syn = - sum_k g_k (V - E_rev_k).
  double I_syn = 0.0;
  I_syn += y[ State_::G1] * ( node.get_E_rev1() - y[State_::V_m] );
  I_syn += y[ State_::G2] * ( node.get_E_rev2() - y[State_::V_m] );
  I_syn += y[ State_::G3] * ( node.get_E_rev3() - y[State_::V_m] );
  I_syn += y[ State_::G4] * ( node.get_E_rev4() - y[State_::V_m] );


  // Total current
  double I_tot = y[State_::I_dep] - y[State_::I_adap] + node.get_I_e() + node.B_.currents_last_value_ + I_syn;

  f[State_::I_input] = node.B_.currents_last_value_ + I_syn;


  // Model currents
  f[State_::I_dep] = ((-node.get_k1())) * y[State_::I_dep];
  f[State_::I_adap] = ((-node.get_k2())) * y[State_::I_adap] + node.get_kadap()*(y[State_::V_m] - node.get_E_L());


  f[State_::V_m] =
      ((1)) / node.get_tau_m() * (y[State_::V_m] - node.get_E_L()) +
      1 / node.get_C_m() * I_tot ;

 // Conductance dynamics
  // Synaptic conductance derivative dG/dt
  f[State_::DG1] = -y[State_::DG1] / node.get_tau_syn1();
  f[State_::G1] = y[ State_::DG1] - y[ State_::G1] / node.get_tau_syn1();

  f[State_::DG2] = -y[State_::DG2] / node.get_tau_syn2();
  f[State_::G2] = y[ State_::DG2] - y[ State_::G2] / node.get_tau_syn2();

  f[State_::DG3] = -y[State_::DG3] / node.get_tau_syn3();
  f[State_::G3] = y[ State_::DG3] - y[ State_::G3] / node.get_tau_syn3();

  f[State_::DG4] = -y[State_::DG4] / node.get_tau_syn4();
  f[State_::G4] = y[ State_::DG4] - y[ State_::G4] / node.get_tau_syn4();

  return GSL_SUCCESS;
}

void eglif_cond_alpha_multisyn::init_buffers_() {
  B_.get_spikes_r1().clear();   // includes resize
  B_.get_spikes_r2().clear();
  B_.get_spikes_r3().clear();
  B_.get_spikes_r4().clear();
  B_.get_currents().clear(); // includes resize
  B_.logger_.reset();     // includes resize
  Archiving_Node::clear_history();

  int state_size = State_::STATE_VEC_SIZE;

  if (B_.s_ == 0)
    B_.s_ = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE);
  else
    gsl_odeiv_step_reset(B_.s_);

  if (B_.c_ == 0) {
    B_.c_ = gsl_odeiv_control_y_new(1e-6, 0.0);
  } else {
    gsl_odeiv_control_init(B_.c_, 1e-6, 0.0, 1.0, 0.0);
  }

  if (B_.e_ == 0) {
    B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
  } else {
    gsl_odeiv_evolve_reset(B_.e_);
  }


  B_.sys_.function = eglif_cond_alpha_multisyn_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = state_size;
  B_.sys_.params = reinterpret_cast<void *>(this);
}

void eglif_cond_alpha_multisyn::calibrate() {
  B_.logger_.init();

  /* generated by template org.nest.nestml.function.Calibrate*/
  srand(time(NULL));
  V_.rng_ = nest::kernel().rng_manager.get_rng( get_thread() );

  V_.RefractoryCounts =
      nest::Time(nest::Time::ms((double)P_.t_ref)).get_steps();

  /* generated by template org.nest.nestml.function.Calibrate*/

  V_.r = 0;

  V_.receptors = 4;			// Number of receptor ports

  V_.V_th = -50.0;

  V_.I_gen = 0.0;


  V_.sum_buffer = 0.0;


  V_.G0.resize( V_.receptors);

  /* generated by template org.nest.nestml.function.Calibrate*/
  V_.G0[0] = 1.0*numerics::e / P_.tau_syn1;
  V_.G0[1] = 1.0*numerics::e / P_.tau_syn2;
  V_.G0[2] = 1.0*numerics::e / P_.tau_syn3;
  V_.G0[3] = 1.0*numerics::e / P_.tau_syn4;

  /* generated by template org.nest.nestml.function.Calibrate*/

  V_.bufferT.resize(P_.bT_size);
  for (long i = 0; i < get_bT_size(); i++) {
    V_.bufferT[i] = 0;				// Non è qui il problema
  }

  /* generated by template org.nest.nestml.function.Calibrate*/

  V_.tot_cur = 0;

  // Receptors initialization
  //B_.spikes.resize(V_.receptors);

  B_.receptor_types_.resize(4);

  for (long i = 0; i < V_.receptors; i++) {
    B_.receptor_types_[i] = i + 1;			// Assign to ports (receptor_type) numbers from 1 to receptors (= receptors number)
  }

}

/* ----------------------------------------------------------------
* Update and spike handling functions: adapted from the file generate by NESTML to update the buffer of external currents
* ---------------------------------------------------------------- */

/*

 */
void eglif_cond_alpha_multisyn::update(nest::Time const &origin,
                                              const long from, const long to) {

  double step_ = nest::Time::get_resolution().get_ms();
  double IntegrationStep_ = nest::Time::get_resolution().get_ms();
  double t = 0;
  double curr_conv_fact = P_.kadap*(-P_.E_L)/P_.k2;

  V_.old_Iinput = V_.new_Iinput;
  V_.new_Iinput = S_.y_[State_::I_input];
  //std::cout << S_.time << std::endl;

  //if (S_.y_[State_::time] < 2){
  if (S_.time < 2){
	  V_.V_th = get_V_th();
	  S_.y_[State_::V_m] = get_Vinit();

    S_.y_[State_::I_adap] = 0;
    S_.y_[State_::I_dep] = 0;
    }

  if ( S_.y_[State_::V_m] < P_.Vmin){
	  S_.y_[State_::V_m] = get_Vmin();

    }
  if ( S_.y_[State_::I_input] < P_.ith){

        S_.y_[State_::I_dep] = 0;

        S_.y_[State_::I_adap] = curr_conv_fact*((S_.y_[State_::I_input] / P_.sc) / (P_.bet - P_.delta1));

        S_.y_[State_::V_m] = (((S_.y_[State_::I_input] / P_.sc) / (P_.bet - P_.delta1) - 1))*(-P_.E_L);

    } else {
      if (V_.new_Iinput < V_.old_Iinput){
            // To Check if S_.y_[State_::I_adap] is ok; it was Iadap_ini
            S_.y_[State_::I_dep] = curr_conv_fact*(S_.y_[State_::I_adap]-(S_.y_[State_::I_input] / P_.sc) / P_.bet);

            }
          }
/*
  if ( V_.I_gen < 10.0){
	  S_.y_[State_::V_m] = get_Vinit();

    }*/




  for (long lag = from; lag < to; ++lag) {
   // std::cout << S_.y_[State_::y_0] << " " << step_ <<std::endl;
    // TODO this case must be handled uniformly, also NESTReferenceConverter
    // must be adopted
    B_.spikes_r1_last_value_ = get_spikes_r1().get_value(lag);
    // TODO this case must be handled uniformly, also NESTReferenceConverter
    // must be adopted
    B_.spikes_r2_last_value_ = get_spikes_r2().get_value(lag);
    // TODO this case must be handled uniformly, also NESTReferenceConverter
    // must be adopted
    B_.spikes_r3_last_value_ = get_spikes_r3().get_value(lag);
    // TODO this case must be handled uniformly, also NESTReferenceConverter
    // must be adopted
    B_.spikes_r4_last_value_ = get_spikes_r4().get_value(lag);

    B_.currents_last_value_ = get_currents().get_value(lag);

    /* generated by template org.nest.spl.Block*/
    /* generated by template org.nest.spl.Statement*/

    // # Aggiorno il buffer delle correnti esterne e la loro somma

    /* generated by template org.nest.spl.SmallStatement*/
    /* generated by template org.nest.spl.Assignment*/

   if (lag == from){
    //S_.y_[State_::time] += 1;
    S_.time += 1;

    V_.I_gen = B_.currents_last_value_;
    V_.tot_cur = B_.currents_last_value_;

    for (int i = 0; i < P_.bT_size; i++) {
      V_.bufferT[i] = V_.bufferT[i+1];
    }
    int index = P_.bT_size-1;
    V_.bufferT[index] = V_.tot_cur;

    if (V_.bufferT[index]>V_.bufferT[index-1])
    {
      V_.init_sign = nest::Time::get_resolution().get_ms();
    }
    //std::cout<<V_.tot_cur<< " " << P_.bT_size <<std::endl;
    /* generated by template org.nest.spl.Statement*/
    //
    /* generated by template org.nest.spl.SmallStatement*/
    /* generated by template org.nest.spl.Assignment*/
    V_.sum_buffer = 0;
    for (int i = 0; i < P_.bT_size; i++) {
     // V_.sum_buffer += V_.bufferT[i];
     V_.sum_buffer += V_.bufferT[i];		// questa riga funziona
     //S_.y_[State_::sum_buffer] += V_.tot_cur;
    }
    //std::cout<<S_.y_[State_::sum_buffer]<<" "<<V_.bufferT[0]<<" "<<V_.bufferT[9]<<std::endl;
    /* generated by template org.nest.spl.Statement*/
    //
    /* generated by template org.nest.spl.CompoundStatement*/
    /* generated by template org.nest.spl.IfStatement*/

    }


    if (V_.r == 0) {

      /* generated by template org.nest.spl.Block*/
      /* generated by template org.nest.spl.Statement*/
      //
      /* generated by template org.nest.spl.SmallStatement*/
      /* generated by template org.nest.spl.FunctionCall*/
      /* generated by template org.nest.spl.GSLIntegrator*/
      t = 0;

      while (t < step_) {
        const int status =
            gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_,
                                   &B_.sys_,          // system of ODE
                                   &t,                // from t
                                   step_,             // to t <= step
                                   &IntegrationStep_, // integration step size
				   S_.y_);
        if (status != GSL_SUCCESS) {
          throw nest::GSLSolverFailure(get_name(), status);
        }
      }

    } else {
      /* generated by template org.nest.spl.Block*/
      /* generated by template org.nest.spl.Statement*/
      //
      /* generated by template org.nest.spl.SmallStatement*/
      /* generated by template org.nest.spl.Assignment*/
      V_.r = V_.r - 1;

    } /* if end */

    /* generated by template org.nest.spl.Statement*/
    //
    /* generated by template org.nest.spl.CompoundStatement*/
    /* generated by template org.nest.spl.IfStatement*/


    //if (S_.y_[State_::V_m] >= S_.y_[State_::V_th]) {		// Condition for spike generation without stochasticity
      /* generated by template org.nest.spl.Block*/
      /* generated by template org.nest.spl.Statement*/
      //
      /* generated by template org.nest.spl.SmallStatement*/
      /* generated by template org.nest.spl.Assignment*/

	if (S_.y_[State_::V_m] > P_.V_th)
        {									/////////// Spike generation!!!! ///////////
	     //std::cout<<"spike!"<<std::endl;
	     spikes_times_.push_back(S_.time);
	      V_.r = V_.RefractoryCounts;

	      // Update of the buffer of spikes times
	      // compute the mean frequency (1/ISI) of the last 3 spikes
	      // Update of V_th based on the frequency of last 3 spikes: V_th = V_th_start+constant*(f-f_tonic); e.g. constant = a, b, c (gi� definiti nel modello)




	      /* generated by template org.nest.spl.Statement*/
	      //
	      /* generated by template org.nest.spl.SmallStatement*/
	      /* generated by template org.nest.spl.Assignment*/
	      S_.y_[State_::V_m] = P_.V_reset;

		    V_.V_th = P_.V_th;


        // Update Iadap

        /* generated by template org.nest.spl.Statement*/
        //
        /* generated by template org.nest.spl.SmallStatement*/
        /* generated by template org.nest.spl.Assignment*/

        double time_scale = 1 / (-P_.sc / (P_.C_m * P_.E_L));

        if (S_.y_[State_::I_input]<P_.istim_min_spiking_exp){
          double c_aux = 0.8*P_.Idep_ini_vr + (S_.y_[State_::I_input]/(P_.sc)) / P_.bet+(P_.delta1/P_.bet)*(1+P_.V_reset)-P_.a*exp(P_.b*S_.y_[State_::I_input]/1000);

          S_.y_[State_::I_adap] = curr_conv_fact*(c_aux + (P_.a * exp(P_.b*S_.y_[State_::I_input]) * ((nest::Time::get_resolution().get_ms()-V_.init_sign) * time_scale)) / (P_.alp + (nest::Time::get_resolution().get_ms()-V_.init_sign) * time_scale));
        } else {
          S_.y_[State_::I_adap] = curr_conv_fact*(P_.c + (P_.a * exp(P_.b*S_.y_[State_::I_input]) * ((nest::Time::get_resolution().get_ms()-V_.init_sign) * time_scale)) / (P_.alp + (nest::Time::get_resolution().get_ms()-V_.init_sign) * time_scale));
        }

        if ((nest::Time::get_resolution().get_ms()*P_.tau_m )>P_.ts){
          S_.y_[State_::I_adap]/=curr_conv_fact;
          S_.y_[State_::I_adap]*=P_.mul;
          S_.y_[State_::I_adap]+=1;
          S_.y_[State_::I_adap]*=curr_conv_fact;
        }

		/* generated by template org.nest.spl.Statement*/
		//
		/* generated by template org.nest.spl.SmallStatement*/
		/* generated by template org.nest.spl.Assignment*/

      // Update Idep

      if (S_.y_[State_::I_input]<1){
        S_.y_[State_::I_adap]=0;
        S_.y_[State_::I_dep] =0;
      } else {
        S_.y_[State_::I_dep]=curr_conv_fact*P_.Idep_ini_vr;
      }

		/* generated by template org.nest.spl.Statement*/
		//



      /* generated by template org.nest.spl.Statement*/
	      //
	      /* generated by template org.nest.spl.CompoundStatement*/
	      /* generated by template org.nest.spl.IfStatement*/


      /* generated by template org.nest.spl.Statement*/
      //
      /* generated by template org.nest.spl.SmallStatement*/
      /* generated by template org.nest.spl.FunctionCall*/
      set_spiketime(nest::Time::step(origin.get_steps() + lag + 1));
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send(*this, se, lag);
      ;
    } /* if (spike generation condition) end */


    // Spikes on 1st receptor
    S_.y_[ State_::DG1 ] += B_.spikes_r1_last_value_ * V_.G0[ 0 ];
   // std::cout << "spike time: "<<B_.spikes_r1_last_value_ << " DG update: "<<S_.y_[ State_::DG1 ]<< " G0_1" << V_.G0[ 0 ] << std::endl;

    // Spikes on 2nd receptor
    S_.y_[ State_::DG2 ] += B_.spikes_r2_last_value_ * V_.G0[ 1 ];

    // Spikes on 3rd receptor
    S_.y_[ State_::DG3 ] += B_.spikes_r3_last_value_ * V_.G0[ 2 ];

    // Spikes on 4th receptor
    S_.y_[ State_::DG4 ] += B_.spikes_r4_last_value_ * V_.G0[ 3 ];

/*
    for ( size_t i = 0; i < V_.get_receptors(); ++i )
    {
      S_.y_[ State_::DG1 + ( 2 * i ) ] +=
        B_.get_spikes()[ i ].get_value(lag) * V_.G0[ i ]; // add incoming spike
    }  */


    // voltage logging
    B_.logger_.record_data(origin.get_steps() + lag);
  }	// end for loop
}	// end function update

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void eglif_cond_alpha_multisyn::handle(nest::DataLoggingRequest &e) {
  B_.logger_.handle(e);
}

void eglif_cond_alpha_multisyn::handle(nest::SpikeEvent &e) {

  assert(e.get_delay_steps() > 0);
/*
  assert(
  ( e.get_rport() > 0 ) && ( ( size_t ) e.get_rport() <= P_.get_receptors() ) );
*/

   //std::cout << "porta spike: "<< e.get_rport() << " ";
   if (e.get_rport()==1) {
      get_spikes_r1().add_value(
          e.get_rel_delivery_steps(
              nest::kernel().simulation_manager.get_slice_origin()),
          e.get_weight() * e.get_multiplicity());
    }else if (e.get_rport()==2){
      get_spikes_r2().add_value(
          e.get_rel_delivery_steps(
              nest::kernel().simulation_manager.get_slice_origin()),
          e.get_weight() * e.get_multiplicity());
    }else if (e.get_rport()==3){
      get_spikes_r3().add_value(
          e.get_rel_delivery_steps(
              nest::kernel().simulation_manager.get_slice_origin()),
          e.get_weight() * e.get_multiplicity());
    }else{
      get_spikes_r4().add_value(
          e.get_rel_delivery_steps(
              nest::kernel().simulation_manager.get_slice_origin()),
          e.get_weight() * e.get_multiplicity());
    }

}

void eglif_cond_alpha_multisyn::handle(nest::CurrentEvent &e) {
  assert(e.get_delay_steps() > 0);

  const double current = e.get_current();
  const double weight = e.get_weight();

  //std::cout<<current<<std::endl;

  // add weighted current; HEP 2002-10-04
  get_currents().add_value(
      e.get_rel_delivery_steps(
          nest::kernel().simulation_manager.get_slice_origin()),
      weight * current);
}
