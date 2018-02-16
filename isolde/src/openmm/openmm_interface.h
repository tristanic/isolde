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
    OpenMM_Thread_Handler(OpenMM::Integrator* integrator)
        : _integrator(integrator)
        {}

    void step_threaded(size_t steps)
    {
        if (_thread_running)
            throw std::runtime_error("Wait for the current thread to complete first!");
        _thread = std::thread(&OpenMM_Thread_Handler::_step_threaded, this, steps);
    }

    void finalize_thread()
    {
        _thread.join();
        _thread_running = false;
    }

    bool thread_finished() const { return _thread_finished; }
    bool thread_running() const { return _thread_running; }

private:
    OpenMM::Integrator* _integrator;
    std::thread _thread;
    bool _thread_running = false;
    bool _thread_finished = true;

    void _step_threaded(size_t steps)
    {
        _thread_running = true;
        _thread_finished = false;
        _integrator->step(steps);
        _thread_finished = true;
    }
};

} //namespace isolde

#endif //ISOLDE_OPENMM
