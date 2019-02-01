/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



#define PYINSTANCE_EXPORT
#include <iostream>
#include "../molc.h"
#include "openmm_interface.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::OpenMM_Thread_Handler>;

using namespace isolde;

OpenMM_Thread_Handler::OpenMM_Thread_Handler(OpenMM::Context* context)
    : _context(context)
{
    _starting_state = _context->getState(OpenMM::State::Positions + OpenMM::State::Velocities);
    _natoms = _starting_state.getPositions().size();
}

void OpenMM_Thread_Handler::_step_threaded(size_t steps, bool smooth)
{
    try
    {
        _thread_except = nullptr;
        auto start = std::chrono::steady_clock::now();
        _thread_running = true;
        _thread_finished = false;
        size_t steps_done = 0;
        _starting_state = _final_state;
        _smoothing = smooth;
        if (!smooth)
            _smoothed_coords.clear();

        for (; steps_done < steps; )
        {
            size_t these_steps, remaining_steps = steps-steps_done;
            if (remaining_steps > STEPS_PER_VELOCITY_CHECK) {
                these_steps = STEPS_PER_VELOCITY_CHECK;
            } else {
                these_steps = remaining_steps;
            }
            integrator().step(these_steps);
            steps_done += these_steps;
            auto state = _context->getState(OpenMM::State::Positions + OpenMM::State::Velocities);
            auto fast = overly_fast_atoms(state.getVelocities());
            if (fast.size() >0)
            {
                std::cerr << fast.size() << " atoms are moving too fast!" << std::endl;
                _unstable = true;
                break;
            }
            if (smooth)
                _apply_smoothing(state);
        }
        _final_state = _context->getState(OpenMM::State::Positions + OpenMM::State::Velocities);
        auto end = std::chrono::steady_clock::now();
        auto loop_time = end-start;
        if (loop_time < _min_time_per_loop)
            std::this_thread::sleep_for(_min_time_per_loop-loop_time);
        _thread_finished = true;
    } catch (...)
    {
        _thread_except = std::current_exception();
        _thread_finished = true;
    }
}

void OpenMM_Thread_Handler::_apply_smoothing(const OpenMM::State& state)
{
    const auto& coords = state.getPositions();
    if (_smoothed_coords.size() == 0) {
        _smoothed_coords = coords;
        return;
    }
    for (size_t i=0; i<_natoms; ++i)
    {
        auto& smoothed = _smoothed_coords[i];
        const auto& current = coords[i];
        smoothed = current * _smoothing_alpha + smoothed * (1-_smoothing_alpha);
    }
}


void OpenMM_Thread_Handler::_minimize_threaded(const double &tolerance, int max_iterations)
{
    // std::cout << "Starting minimization with tolerance of " << tolerance << " and max iterations per round of " << max_iterations << std::endl;
    try
    {
        _thread_except = nullptr;
        auto start = std::chrono::steady_clock::now();
        _clash = false;
        _smoothed_coords.clear();
        _thread_running = true;
        _thread_finished = false;
        _starting_state = _context->getState(OpenMM::State::Positions + OpenMM::State::Energy);
        _unstable = true;
        double tol = tolerance * _natoms;
        // std::cout << "Initial energy: " << _starting_state.getPotentialEnergy() << " kJ/mol" << std::endl;
        for (size_t i=0; i<1; ++i) {
            // std::cout << "Round " << i << ": " <<std::flush;
            OpenMM::LocalEnergyMinimizer::minimize(*_context, tol, max_iterations);
            _final_state = _context->getState(OpenMM::State::Positions + OpenMM::State::Forces + OpenMM::State::Energy);
            std::cout << _final_state.getPotentialEnergy() << " kJ/mol" << std::endl << std::flush;
            if (_starting_state.getPotentialEnergy() - _final_state.getPotentialEnergy() < tol) {
                _unstable = false;
                break;
            }
        }
        if (max_force(_final_state.getForces()) > MAX_FORCE)
            _clash = true;
        auto end = std::chrono::steady_clock::now();
        auto loop_time = end-start;
        if (loop_time < _min_time_per_loop)
            std::this_thread::sleep_for(_min_time_per_loop-loop_time);
        _thread_finished = true;
        // std::cout << "Finished minimization round" << std::endl;
    } catch (...)
    {
        _thread_except = std::current_exception();
        _thread_finished = true;
    }
}

void OpenMM_Thread_Handler::_reinitialize_context_threaded()
{
    try
    {
        _thread_except = nullptr;
        _thread_running = true;
        _thread_finished = false;
        OpenMM::State current_state = _context->getState(OpenMM::State::Positions + OpenMM::State::Velocities);
        _context->reinitialize();
        _context->setPositions(current_state.getPositions());
        _context->setVelocities(current_state.getVelocities());
        _thread_finished = true;
    } catch (...)
    {
        _thread_except = std::current_exception();
        _thread_finished = true;
    }
}

std::vector<size_t> OpenMM_Thread_Handler::overly_fast_atoms(const std::vector<OpenMM::Vec3>& velocities)
{
    std::vector<size_t> fast_indices;
    size_t i=0;
    for (const auto &v: velocities) {
        if (sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) > MAX_VELOCITY)
            fast_indices.push_back(i);
        i++;
    }
    return fast_indices;
}

std::vector<OpenMM::Vec3> OpenMM_Thread_Handler::get_coords_in_angstroms(const OpenMM::State& state)
{
    finalize_thread();
    const auto &coords_nm = state.getPositions();
    std::vector<OpenMM::Vec3> coords_ang(_natoms);
    auto from = coords_nm.begin();
    auto to = coords_ang.begin();
    for (; from != coords_nm.end(); from++, to++)
    {
        *to = *from * 10.0;
    }
    return coords_ang;
}

std::vector<OpenMM::Vec3> OpenMM_Thread_Handler::get_smoothed_coords_in_angstroms()
{
    finalize_thread();
    if (!_smoothing) {
        throw std::logic_error("Last round of equilibration was not run with smoothing enabled!");
    }
    std::vector<OpenMM::Vec3> coords_ang(_natoms);
    auto from = _smoothed_coords.begin();
    auto to = coords_ang.begin();
    for (; from != _smoothed_coords.end(); from++, to++)
    {
        *to = *from * 10.0;
    }
    return coords_ang;
}

void OpenMM_Thread_Handler::set_coords_in_angstroms(const std::vector<OpenMM::Vec3>& coords_ang)
{
    finalize_thread();
    std::vector<OpenMM::Vec3> coords_nm(_natoms);
    auto from = coords_ang.begin();
    auto to = coords_nm.begin();
    for (; from != coords_nm.end(); from++, to++)
    {
        for (size_t i=0; i<3; ++i) {
            (*to)[i] = (*from)[i]/10.0;
        }
    }
    _context->setPositions(coords_nm);
    _smoothed_coords.clear();
}

void OpenMM_Thread_Handler::set_coords_in_angstroms(double *coords, size_t n)
{
    finalize_thread();
    if (n != natoms())
        throw std::logic_error("Number of input atoms does not match number in simulation!");
    std::vector <OpenMM::Vec3> coords_nm(n);
    for (auto &v: coords_nm) {
        for (size_t i=0; i<3; ++i)
            v[i] = (*coords++)/10.0;
    }
    _context->setPositions(coords_nm);
    _smoothed_coords.clear();
}

double OpenMM_Thread_Handler::max_force(const std::vector<OpenMM::Vec3>& forces) const
{
    double max_force = 0;

    for (const auto &fv: forces)
    {
        double f_mag = 0;
        for (size_t i=0; i<3; ++i)
            f_mag += fv[i]*fv[i];
        f_mag = sqrt(f_mag);
        max_force = f_mag>max_force ? f_mag : max_force;
    }
    return max_force;
}

void OpenMM_Thread_Handler::set_smoothing_alpha(const double &alpha)
{
    if (alpha < SMOOTHING_ALPHA_MIN)
        _smoothing_alpha = SMOOTHING_ALPHA_MIN;
    else if (alpha > SMOOTHING_ALPHA_MAX)
        _smoothing_alpha = SMOOTHING_ALPHA_MAX;
    else
        _smoothing_alpha = alpha;
}

// PYTHON INTERFACE BELOW

SET_PYTHON_INSTANCE(openmm_thread_handler, OpenMM_Thread_Handler)
GET_PYTHON_INSTANCES(openmm_thread_handler, OpenMM_Thread_Handler)

extern "C" EXPORT void*
openmm_thread_handler_new(void *context)
{
    OpenMM::Context *c = static_cast<OpenMM::Context *>(context);
    try {
        OpenMM_Thread_Handler *h = new OpenMM_Thread_Handler(c);
        return h;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
openmm_thread_handler_delete(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        delete h;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
openmm_thread_handler_num_atoms(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        return h->natoms();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT double
openmm_thread_handler_smoothing_alpha(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        return h->smoothing_alpha();
    } catch (...) {
        molc_error();
        return -1;
    }
}

extern "C" EXPORT void
set_openmm_thread_handler_smoothing_alpha(void *handler, double alpha)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->set_smoothing_alpha(alpha);
    } catch (...) {
        molc_error();
    }
}



extern "C" EXPORT void
openmm_thread_handler_step(void *handler, size_t steps, npy_bool smooth)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->step_threaded(steps, smooth);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
openmm_thread_handler_minimize(void *handler, double tolerance, int max_iterations)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->minimize_threaded(tolerance, max_iterations);
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT npy_bool
openmm_thread_handler_thread_finished(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        return h->thread_finished();
    } catch (...) {
        molc_error();
        return false;
    }
}

extern "C" EXPORT void
openmm_thread_handler_finalize_thread(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->finalize_thread();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT npy_bool
openmm_thread_handler_unstable(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        return h->unstable();
    } catch(...) {
        molc_error();
        return false;
    }
}

extern "C" EXPORT npy_bool
openmm_thread_handler_clashing(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        return h->clash_detected();
    } catch(...) {
        molc_error();
        return false;
    }
}

extern "C" EXPORT void
openmm_thread_handler_unstable_atoms(void *handler, size_t n, npy_bool *unstable)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        auto fast_indices = h->overly_fast_atoms(h->final_state().getVelocities());
        for (auto i: fast_indices)
            unstable[i] = true;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
openmm_thread_handler_last_coords(void *handler, size_t n, double *coords)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        auto sim_coords = h->get_coords_in_angstroms(h->initial_state());
        if (n != h->natoms())
            throw std::logic_error("Mismatch between number of atoms and output array size!");
        for (const auto &coord: sim_coords) {
            for (size_t i=0; i<3; ++i)
                *coords++ = coord[i];
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
openmm_thread_handler_current_coords(void *handler, size_t n, double *coords)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        auto sim_coords = h->get_coords_in_angstroms(h->final_state());
        if (n != h->natoms())
            throw std::logic_error("Mismatch between number of atoms and output array size!");
        for (const auto &coord: sim_coords) {
            for (size_t i=0; i<3; ++i)
                *coords++ = coord[i];
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
openmm_thread_handler_smoothed_coords(void *handler, size_t n, double *coords)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        if (n != h->natoms())
            throw std::logic_error("Mismatch between number of atoms and output array size!");
        auto sim_coords = h->get_smoothed_coords_in_angstroms();
        for (const auto &coord: sim_coords) {
            for (size_t i=0; i<3; ++i)
                *coords++ = coord[i];
        }
    } catch (...) {
        molc_error();
    }

}

extern "C" EXPORT void
set_openmm_thread_handler_current_coords(void *handler, size_t n, double *coords)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->set_coords_in_angstroms(coords, n);
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT double
openmm_thread_handler_min_thread_period(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        return h->get_minimum_thread_time_in_ms();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
set_openmm_thread_handler_min_thread_period(void *handler, double time_ms)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->set_minimum_thread_time_in_ms(time_ms);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
openmm_thread_handler_reinitialize_context_and_keep_state(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->reinitialize_context_and_keep_state();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
openmm_thread_handler_reinitialize_context_and_keep_state_threaded(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->reinitialize_context_threaded();
    } catch (...) {
        molc_error();
    }
}


// extern "C" EXPORT void
// openmm_thread_handler_initial_positions(void *handler, size_t n, double *coords)
