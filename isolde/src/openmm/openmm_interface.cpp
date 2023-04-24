/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "openmm_interface.h"
#include "minimize.h"


using namespace isolde;

OpenmmThreadHandler::OpenmmThreadHandler(OpenMM::Context* context)
    : _context(context)
{
    _starting_state = _context->getState(OpenMM::State::Positions | OpenMM::State::Velocities);
    _natoms = _starting_state.getPositions().size();
    _smoothed_coords.resize(_natoms);
    _atom_fixed.resize(_natoms);
    const auto& system = _context->getSystem();
    for (size_t i=0; i<system.getNumParticles(); ++i)
    {
        _atom_fixed[i] = (system.getParticleMass(i)==0.0);
    }
}

void OpenmmThreadHandler::_step_threaded(size_t steps)
{
    try
    {
        // _thread_except = nullptr;
        auto start = std::chrono::steady_clock::now();
        // _thread_running = true;
        // _thread_finished = false;
        size_t steps_done = 0;
        _starting_state = _final_state;

        while (steps_done < steps)
        {
            size_t these_steps, remaining_steps = steps-steps_done;
            auto& integr = integrator();
            integr.setCurrentIntegrator(MAIN);
            if (remaining_steps > STEPS_PER_VELOCITY_CHECK) {
                these_steps = STEPS_PER_VELOCITY_CHECK;
                integr.step(STEPS_PER_VELOCITY_CHECK);
                integr.setCurrentIntegrator(VELOCITY_CHECK);
                integr.step(1);
                auto& vcheck = static_cast<OpenMM::CustomIntegrator&>(integr.getIntegrator(VELOCITY_CHECK));
                if (vcheck.getGlobalVariable(0) > 0.0) {
                    std::cerr << int(vcheck.getGlobalVariableByName("fast_count")) << " atoms are moving too fast!" << std::endl;
                    _unstable = true;
                    break;
                }
                integr.setCurrentIntegrator(SMOOTH);
                integr.step(1);
            } else {
                these_steps = remaining_steps;
                integr.step(these_steps);
            }

            steps_done += these_steps;
        }
        _final_state = _context->getState(OpenMM::State::Positions | OpenMM::State::Velocities);
        _apply_smoothing(static_cast<OpenMM::CustomIntegrator&>(integrator().getIntegrator(SMOOTH)), _final_state);
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

void OpenmmThreadHandler::_apply_smoothing(OpenMM::CustomIntegrator& igr, const OpenMM::State& state)
{
    igr.getPerDofVariableByName("smoothed", _smoothed_coords);
    // Integrators ignore fixed atoms, so we need to copy those coords from the State
    const auto& unsmoothed = state.getPositions();
    for (size_t i=0; i<_natoms; ++i)
    {
        if (_atom_fixed[i])
        {
            _smoothed_coords[i] = unsmoothed[i];
        }
    }

}


void OpenmmThreadHandler::_minimize_threaded(const double &tolerance, int max_iterations)
{
    // std::cout << "Starting minimization with tolerance of " << tolerance << " and max iterations per round of " << max_iterations << std::endl;
    try
    {
        // _thread_except = nullptr;
        auto start = std::chrono::steady_clock::now();
        _clash = false;
        reset_smoothing();
        // _thread_running = true;
        // _thread_finished = false;
        _starting_state = _context->getState(OpenMM::State::Positions | OpenMM::State::Energy);
        _min_converged = false;
        double tol = tolerance * _natoms;
        // std::cout << "Initial energy: " << _starting_state.getPotentialEnergy() << " kJ/mol" << std::endl;
        auto result = isolde::LocalEnergyMinimizer::minimize(*_context, tol, max_iterations);
        if (result == isolde::LocalEnergyMinimizer::SUCCESS)
        {
            // Minimisation has converged to within the desired tolerance,
            // and all constraints are satisfied.
            _final_state = _context->getState(OpenMM::State::Positions | OpenMM::State::Forces | OpenMM::State::Energy);
            _min_converged = true;
            _unstable = false;
        } else if (result == isolde::LocalEnergyMinimizer::DID_NOT_CONVERGE) {
            // Minimisation ongoing. Just leave _min_converged = false, but
            // let ISOLDE have the new coordinates.
            _final_state = _context->getState(OpenMM::State::Positions | OpenMM::State::Forces | OpenMM::State::Energy);
        } else // if (result < 0)
        {
            // Minimisation failed. Revert the model to its initial state
            // and let ISOLDE point out problem areas to the user.
            _clash = true;
            _final_state = _context->getState(OpenMM::State::Positions | OpenMM::State::Forces | OpenMM::State::Energy);
            // _final_state = _starting_state;
        }
        //if (_min_converged && max_force(_final_state.getForces()) > MAX_FORCE)
        if (_min_converged && max_force(_context->getSystem(), _final_state) > MAX_FORCE)
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

void OpenmmThreadHandler::_reinitialize_context_threaded()
{
    try
    {
        _thread_except = nullptr;
        // _thread_running = true;
        // _thread_finished = false;
        OpenMM::State current_state = _context->getState(OpenMM::State::Positions | OpenMM::State::Velocities);
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

std::vector<size_t> OpenmmThreadHandler::overly_fast_atoms(const std::vector<OpenMM::Vec3>& velocities) const
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

std::vector<OpenMM::Vec3> OpenmmThreadHandler::get_coords_in_angstroms(const OpenMM::State& state)
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

std::vector<OpenMM::Vec3> OpenmmThreadHandler::get_smoothed_coords_in_angstroms()
{
    finalize_thread();
    std::vector<OpenMM::Vec3> coords_ang(_natoms);
    auto from = _smoothed_coords.begin();
    auto to = coords_ang.begin();
    for (; from != _smoothed_coords.end(); from++, to++)
    {
        *to = *from * 10.0;
    }
    return coords_ang;
}

void OpenmmThreadHandler::set_coords_in_angstroms(const std::vector<OpenMM::Vec3>& coords_ang)
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
    reset_smoothing();
}

void OpenmmThreadHandler::set_coords_in_angstroms(double *coords, size_t n)
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
    reset_smoothing();
}

double OpenmmThreadHandler::max_force(const std::vector<OpenMM::Vec3>& forces) const
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

// get maximum force, ignoring massless (fixed) particles
double OpenmmThreadHandler::max_force(const OpenMM::System& system, const OpenMM::State& state) const
{
    double max_force = 0;
    auto forces = state.getForces();
    for (size_t i=0; i<system.getNumParticles(); ++i)
    {
        if (_atom_fixed[i]) continue;
        const auto& f = forces[i];
        double f_mag = 0;
        for (size_t j=0; j<3; ++j)
            f_mag += f[j]*f[j];
        f_mag = sqrt(f_mag);
        if (f_mag > max_force)
            max_force = f_mag;
    }
    return max_force;
}


void OpenmmThreadHandler::set_smoothing_alpha(const double &alpha)
{
    if (alpha < SMOOTHING_ALPHA_MIN)
        _smoothing_alpha = SMOOTHING_ALPHA_MIN;
    else if (alpha > SMOOTHING_ALPHA_MAX)
        _smoothing_alpha = SMOOTHING_ALPHA_MAX;
    else
        _smoothing_alpha = alpha;
}

void OpenmmThreadHandler::reset_smoothing()
{
    auto& integr = static_cast<OpenMM::CustomIntegrator&>(integrator().getIntegrator(SMOOTH));
    integr.setGlobalVariableByName("reset_smooth", 1.0);
}

// PYTHON INTERFACE BELOW


namespace py=pybind11;
using OTH = OpenmmThreadHandler;

PYBIND11_MODULE(_openmm_async, m) {
    m.doc() = "Manager for running and updating an OpenMM simulation in an asynchronous "
        "threaded manner.";
    py::class_<OTH>(m, "OpenmmThreadHandler")
        .def(py::init([](uintptr_t context)
        {
            auto c = reinterpret_cast<OpenMM::Context*>(context);
            OTH* th = new OTH(c);
            return std::unique_ptr<OTH>(th);
        }))
        .def_property_readonly("num_atoms", &OTH::natoms)
        .def_property("smoothing_alpha", &OTH::smoothing_alpha, &OTH::set_smoothing_alpha)
        .def("step", &OTH::step_threaded)
        .def("minimize", &OTH::minimize_threaded)
        .def("thread_finished", &OTH::thread_finished)
        .def("finalize_thread", &OTH::finalize_thread)
        .def("unstable", &OTH::unstable)
        .def_property_readonly("minimization_converged", &OTH::converged)
        .def_property_readonly("clashing", &OTH::clash_detected)
        .def_property_readonly("unstable_atoms", [](const OTH& self){
            auto fast_indices = self.overly_fast_atoms(self.final_state().getVelocities());
            auto ret = py::array_t<bool>(self.natoms());
            ret[py::make_tuple(py::ellipsis())]=false;
            bool* ptr = (bool*)ret.request().ptr;
            for (const auto& i: fast_indices) ptr[i]=true;
            return ret;
        })
        .def_property_readonly("last_coords", [](OTH& self) {
            auto ret = py::array_t<double>({self.natoms(), (size_t)3});
            double* ptr = (double*)ret.request().ptr;
            auto sim_coords = self.get_coords_in_angstroms(self.initial_state());
            for (const auto& coord: sim_coords) {
                for (size_t i=0; i<3; ++i)
                    *ptr++ = coord[i];
            }
            return ret;
        })
        .def_property("current_coords", 
        [](OTH& self) {
            auto ret = py::array_t<double>({self.natoms(), (size_t)3});
            double* ptr = (double*)ret.request().ptr;
            auto sim_coords = self.get_coords_in_angstroms(self.final_state());
            for (const auto& coord: sim_coords) {
                for (size_t i=0; i<3; ++i)
                    *ptr++ = coord[i];
            }
            return ret;
        },
        [](OTH& self, const py::array_t<double> coords) {
            double* ptr = (double*)coords.request().ptr;
            self.set_coords_in_angstroms(ptr, coords.shape(0));
        }
        
        )
        .def_property_readonly("smoothed_coords", [](OTH& self) {
            auto ret = py::array_t<double>({self.natoms(), (size_t)3});
            double* ptr = (double*)ret.request().ptr;
            auto sim_coords = self.get_smoothed_coords_in_angstroms();
            for (const auto& coord: sim_coords) {
                for (size_t i=0; i<3; ++i)
                    *ptr++ = coord[i];
            }
            return ret;
        })
        .def_property("min_thread_period", &OTH::get_minimum_thread_time_in_ms, &OTH::set_minimum_thread_time_in_ms)
        .def("reinitialize_context_and_keep_state", &OTH::reinitialize_context_and_keep_state)
        .def("reinitialize_context_and_keep_state_threaded", &OTH::reinitialize_context_threaded)
        ;

};
