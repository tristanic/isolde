#ifndef ISOLDE_OPENMM
#define ISOLDE_OPENMM

#include <thread>
#include <cstddef>
#include <OpenMM.h>
#include <pyinstance/PythonInstance.declare.h>

namespace isolde
{

class OpenMM_Thread_Handler: public pyinstance::PythonInstance<OpenMM_Thread_Handler>
{
public:
    OpenMM_Thread_Handler() {}
    ~OpenMM_Thread_Handler() { if (_thread_running) _thread.join(); }
    /*! Rather annoyingly, we have to set the temperature explicitly here
     *  since the Integrator base class doesn't provide a virtual
     *  getTemperature() method
     *
     */
    OpenMM_Thread_Handler(OpenMM::Context* context);

    OpenMM::Integrator& integrator() { return _context->getIntegrator();}

    /*! Runs the desired number of steps in a separate CPU thread, checking every
     *  ten steps to make sure that velocities remain under control. If excessive
     *  atomic velocities are detected, the thread will stop at that point and
     *  set unstable() to true. The initial and final states (containing coordinates
     *  and velocities) are accessible via initial_state() and final_state()
     *  respectively. Instability can be checked for using unstable(). The
     *  indices of any excessively fast-moving atoms can be retrieved using
     *  overly_fast_atoms().
     */
    void step_threaded(size_t steps)
    {
        _thread_safety_check();
        if (_unstable)
            throw std::logic_error("The last round had atoms moving dangerously fast. Fix the issues and minimise first.");
        _thread = std::thread(&OpenMM_Thread_Handler::_step_threaded, this, steps);
    }

    void minimize_threaded()
    {
        _thread_safety_check();
        _thread = std::thread(&OpenMM_Thread_Handler::_minimize_threaded, this);
    }

    std::vector<size_t> overly_fast_atoms(const std::vector<OpenMM::Vec3>& velocities);

    std::vector<OpenMM::Vec3> get_coords_in_angstroms(const OpenMM::State& state);
    void set_coords_in_angstroms(const std::vector<OpenMM::Vec3>& coords_ang);
    void finalize_thread()
    {
        if (_thread_running)
            _thread.join();
        _thread_running = false;
    }

    const OpenMM::State& initial_state() const { _thread_finished_check(); return _starting_state; }
    const OpenMM::State& final_state() const { _thread_finished_check(); return _final_state; }

    bool thread_finished() const { return _thread_finished; }
    bool thread_running() const { return _thread_running; }
    bool unstable() const { return _unstable; }
    size_t natoms() const { return _natoms; }

private:
    OpenMM::Context* _context;
    OpenMM::State _starting_state;
    OpenMM::State _final_state;
    std::thread _thread;
    size_t _natoms;
    bool _thread_running = false;
    bool _thread_finished = true;
    bool _unstable = false;

    const double MAX_VELOCITY = 15; //nm ps-1 (15,000 m/s)
    const double MIN_TOLERANCE = 1.0; //kJ mol-1
    const size_t MAX_MIN_ITERATIONS = 500;
    const size_t STEPS_PER_VELOCITY_CHECK = 10;

    void _thread_safety_check() const {
        if (_thread_running) {
            throw std::logic_error("You need to finalize the old thread first!");
        }
    }
    void _thread_finished_check() const {
        if (!_thread_finished) {
            throw std::logic_error("This function is not available while a thread is running!");
        }
    }
    void _step_threaded(size_t steps);
    void _minimize_threaded();
};

} //namespace isolde

#endif //ISOLDE_OPENMM
