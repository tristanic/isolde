#define PYINSTANCE_EXPORT
#include "../molc.h"
#include "openmm_interface.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::OpenMM_Thread_Handler>;

using namespace isolde;

SET_PYTHON_INSTANCE(openmm_thread_handler, OpenMM_Thread_Handler)
GET_PYTHON_INSTANCES(openmm_thread_handler, OpenMM_Thread_Handler)

extern "C" EXPORT void*
openmm_thread_handler_new(void *integrator)
{
    OpenMM::Integrator *i = static_cast<OpenMM::Integrator *>(integrator);
    try {
        OpenMM_Thread_Handler *h = new OpenMM_Thread_Handler(i);
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
