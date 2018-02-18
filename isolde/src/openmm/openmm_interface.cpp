#define PYINSTANCE_EXPORT
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

void OpenMM_Thread_Handler::_step_threaded(size_t steps)
{
    _thread_running = true;
    _thread_finished = false;
    size_t steps_done = 0;
    _starting_state = _final_state;
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
        auto state = _context->getState(OpenMM::State::Velocities);
        if (overly_fast_atoms(state.getVelocities()).size() >0)
        {
            _unstable = true;
            break;
        }
    }
    _final_state = _context->getState(OpenMM::State::Positions + OpenMM::State::Velocities);
    _thread_finished = true;
}

void OpenMM_Thread_Handler::_minimize_threaded()
{
    _thread_running = true;
    _thread_finished = false;
    _starting_state = _context->getState(OpenMM::State::Positions + OpenMM::State::Energy);
    OpenMM::LocalEnergyMinimizer::minimize(*_context, MIN_TOLERANCE, MAX_MIN_ITERATIONS);
    _final_state = _context->getState(OpenMM::State::Positions + OpenMM::State::Energy);
    if (_starting_state.getPotentialEnergy() - _final_state.getPotentialEnergy() > MIN_TOLERANCE)
        _unstable = true;
    else
        _unstable = false;
    _thread_finished = true;
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
    for (; from != coords_ang.end(); from++, to++)
    {
        for (size_t i=0; i<3; ++i)
            (*to)[i] = (*from)[i]*10.0;
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

extern "C" EXPORT void
openmm_thread_handler_step(void *handler, size_t steps)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->step_threaded(steps);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
openmm_thread_handler_minimize(void *handler)
{
    OpenMM_Thread_Handler *h = static_cast<OpenMM_Thread_Handler *>(handler);
    try {
        h->minimize_threaded();
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



// extern "C" EXPORT void
// openmm_thread_handler_initial_positions(void *handler, size_t n, double *coords)
