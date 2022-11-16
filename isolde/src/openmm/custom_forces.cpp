/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */




#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <OpenMM.h>

namespace py=pybind11;
PYBIND11_MODULE(_openmm_force_ext, m) {
    m.doc() = "Extra functions providing C++ array-style getters/setters for select OpenMM custom forces.";

    m.def("cmaptorsionforce_add_torsions", [](uintptr_t force, py::array_t<int> map_indices, py::array_t<int> phi_indices, py::array_t<int> psi_indices) {
        auto f = reinterpret_cast<OpenMM::CMAPTorsionForce *>(force);
        auto n = map_indices.shape(0);
        py::array_t<int> force_indices(n);
        auto fi = (int*)force_indices.request().ptr;
        for (size_t i=0; i<n; ++i) {
            *(fi++) = f->addTorsion(map_indices.at(i), phi_indices.at(i,0), phi_indices.at(i,1), phi_indices.at(i,2), phi_indices.at(i,3), psi_indices.at(i,0), psi_indices.at(i,1), psi_indices.at(i,2), psi_indices.at(i,3));
        }
        return force_indices;
    })
    .def("customcompoundbondforce_add_bonds", [](uintptr_t force, py::array_t<int>p_indices, py::array_t<double> params) {
        auto f = reinterpret_cast<OpenMM::CustomCompoundBondForce *>(force);
        int n_params = f->getNumPerBondParameters();
        std::vector<double> param_vec(n_params);
        int n_particles = f->getNumParticlesPerBond();
        std::vector<int> particles(n_particles);
        auto n = p_indices.shape(0);
        py::array_t<int> force_indices(n);
        auto fi = (int*)force_indices.request().ptr; 
        for (size_t i=0; i<n; ++i)
        {
            for (int j=0; j<n_params; ++j)
                param_vec[j]=params.at(i,j);
            if (n_particles > 1) {
                for (int j=0; j<n_particles; ++j)
                    particles[j] = p_indices.at(i,j);
            } else {
                particles[0] = p_indices.at(i);
            }
            *(fi++) = f->addBond(particles, param_vec);
        }
        return force_indices;
    })
    .def("customcompoundbondforce_update_bond_parameters", [](uintptr_t force, py::array_t<int> indices, py::array_t<double> params) {
        auto f = reinterpret_cast<OpenMM::CustomCompoundBondForce *>(force);
        auto n_params = f->getNumPerBondParameters();
        std::vector<double> param_vec(n_params);
        std::vector<int> particles(f->getNumParticlesPerBond());
        for (size_t i=0; i < indices.shape(0); ++i)
        {
            auto index = indices.at(i);
            f->getBondParameters(index, particles, param_vec);
            for (int j=0; j<n_params; ++j)
                param_vec[j] = params.at(i,j);
            f->setBondParameters(index, particles, param_vec);
        }
    })
    .def("custombondforce_add_bonds", [](uintptr_t force, py::array_t<int> p_indices, py::array_t<double> params) {
        auto f = reinterpret_cast<OpenMM::CustomBondForce *>(force);
        auto n_params = f->getNumPerBondParameters();
        std::vector<double> param_vec(n_params);
        int p1, p2;
        auto n = p_indices.shape(0);
        py::array_t<int> force_indices(n);
        auto fi = (int*)force_indices.request().ptr;
        for (size_t i=0; i<n; ++i)
        {
            for (int j=0; j<n_params; ++j)
                param_vec[j] = params.at(i,j);
            p1 = p_indices.at(i,0);
            p2 = p_indices.at(i,1);
            *(fi++) = f->addBond(p1, p2, param_vec);
        }
        return force_indices;
    })
    .def("custombondforce_update_bond_parameters", [](uintptr_t force, py::array_t<int> indices, py::array_t<double> params) {
        auto f = reinterpret_cast<OpenMM::CustomBondForce *>(force);
        int n_params = f->getNumPerBondParameters();
        std::vector<double> param_vec(n_params);
        int p1, p2;
        auto n = indices.shape(0);
        for (size_t i=0; i<n; ++i)
        {
            int index = indices.at(i);
            f-> getBondParameters(index, p1, p2, param_vec);
            for (int j=0; j<n_params; ++j)
                param_vec[j] = params.at(i,j);
            f->setBondParameters(index, p1, p2, param_vec);
        }
    })
    .def("customexternalforce_add_particles", [](uintptr_t force, py::array_t<int> particle_indices, py::array_t<double> params) {
        auto f = reinterpret_cast<OpenMM::CustomExternalForce *>(force);
        auto n_params = f->getNumPerParticleParameters();
        auto n = particle_indices.shape(0);
        std::vector<double> param_vec(n_params);
        py::array_t<int> force_indices(n);
        auto fi = (int*)force_indices.request().ptr;
        for (size_t i=0; i<n; ++i)
        {
            for (int j=0; j<n_params; ++j)
                param_vec[j] = params.at(i,j);
            *(fi++) = f->addParticle(particle_indices.at(i), param_vec);
        }
        return force_indices;
    })
    .def("customexternalforce_update_particle_parameters", [](uintptr_t force, py::array_t<int> indices, py::array_t<double> params) {
        auto f = reinterpret_cast<OpenMM::CustomExternalForce *>(force);
        auto n_params = f->getNumPerParticleParameters();
        auto n = indices.shape(0);
        int particle;
        std::vector<double> param_vec(n_params);
        for (size_t i=0; i<n; ++i) {
            int index = indices.at(i);
            f->getParticleParameters(index, particle, param_vec);
            for (int j=0; j<n_params; ++j)
                param_vec[j] = params.at(i,j);
            f->setParticleParameters(index, particle, param_vec);
        }
    })
    .def("customtorsionforce_add_torsions", [](uintptr_t force, py::array_t<int> particle_indices, py::array_t<double> params) {
        auto f = reinterpret_cast<OpenMM::CustomTorsionForce *>(force);
        auto n_params = f->getNumPerTorsionParameters();
        auto n = particle_indices.shape(0);
        std::vector<double> param_vec(n_params);
        py::array_t<int> force_indices(n);
        auto fi = (int*)force_indices.request().ptr;
        auto pi = particle_indices;
        for (size_t i=0; i<n; ++i)
        {
            for (int j=0; j<n_params; ++j)
                param_vec[j] = params.at(i,j);
            *(fi++) = f->addTorsion(pi.at(i,0), pi.at(i,1), pi.at(i,2), pi.at(i,3), param_vec);
        }
        return force_indices;
    })
    .def("customtorsionforce_update_torsion_parameters", [](uintptr_t force, py::array_t<int> indices, py::array_t<double> params) {
        auto f = reinterpret_cast<OpenMM::CustomTorsionForce *>(force);
        auto n_params = f->getNumPerTorsionParameters();
        auto n = indices.shape(0);
        std::vector<double> param_vec(n_params);
        int p1, p2, p3, p4;
        for (size_t i=0; i<n; ++i)
        {
            auto index = indices.at(i);
            f->getTorsionParameters(index, p1, p2, p3, p4, param_vec);
            for (int j=0; j<n_params; ++j) {
                param_vec[j] = params.at(i,j);
            }
            f->setTorsionParameters(index, p1, p2, p3, p4, param_vec);
        }
    })
    ;
};

