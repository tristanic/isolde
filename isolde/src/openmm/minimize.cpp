#include "minimize.h"
#include "lbfgs.h"
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace OpenMM;

struct MinimizerData
{
    Context& context;
    double k;
    bool checkLargeForces;
    bool largeForceEncountered=false;
    MinimizerData(Context& context, double k) : context(context), k(k) {
        std::string platformName = context.getPlatform().getName();
        checkLargeForces = (platformName == "CUDA" || platformName == "OpenCL");
    }
    ~MinimizerData() {}
};

static double computeForcesAndEnergy(Context& context,
        const std::vector<Vec3>& positions, lbfgsfloatval_t *g)
{
    context.setPositions(positions);
    context.computeVirtualSites();
    State state = context.getState(State::Forces | State::Energy);
    const std::vector<Vec3>& forces = state.getForces();
    const System& system = context.getSystem();
    for (int i = 0; i < forces.size(); i++) {
        if (system.getParticleMass(i) == 0) {
            g[3*i] = 0.0;
            g[3*i+1] = 0.0;
            g[3*i+2] = 0.0;
        }
        else {
            g[3*i] = -forces[i][0];
            g[3*i+1] = -forces[i][1];
            g[3*i+2] = -forces[i][2];
        }
    }
    return state.getPotentialEnergy();
}

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
    MinimizerData* data = reinterpret_cast<MinimizerData*>(instance);
    Context& context = data->context;
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();

    // Compute the force and energy for this configuration.

    std::vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++)
        positions[i] = Vec3(x[3*i], x[3*i+1], x[3*i+2]);
    double energy = computeForcesAndEnergy(context, positions, g);
    if (data->checkLargeForces) {
        // The CUDA and OpenCL platforms accumulate forces in fixed point, so they
        // can't handle very large forces.  If we encounter an infinite or NaN
        // value, immediately drop out with an error - something is *very* wrong
        // with the model. The maximum force the GPU can calculate is a little
        // over 2e9 kJ/mol/nm. If we encounter a force greater than 2e9 it's a
        // sign of trouble, but not necessarily disastrous - as long as such
        // cases are limited to only a few positions along a given line search
        // the minimiser will almost always recover gracefully. The default
        // OpenMM minimiser would fall back to calculating forces on the CPU
        // in this case, but this is unworkable for ISOLDE - in the presence of
        // MDFF, energy calculations on the CPU take on the order of minutes
        // rather than a few seconds, not conducive to an interactive
        // application. Instead, we'll just note that the out-of-range force
        // occurred, and then let the user know only if the minimiser fails to
        // converge.

        for (int i = 0; i < 3*numParticles; i++) {
            lbfgsfloatval_t grad = g[i];
            // If infinite or NaN, throw our hands up immediately. Something's
            // badly wrong.
            if (grad - grad != 0) {
                throw std::out_of_range("Infinite or NaN value encountered!");
            }
            // Otherwise, just note that an out-of-range force occurred, but see
            // if we can successfully converge anyway (this is usually true if
            // it's only one or two entries on a given line search).
            else if (!(std::fabs(grad) < 2e9)) {
                data->largeForceEncountered = true;
            }
        }

    }

    // Add harmonic forces for any constraints.

    int numConstraints = system.getNumConstraints();
    const double& k = data->k;
    for (int i = 0; i < numConstraints; i++) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        Vec3 delta = positions[particle2]-positions[particle1];
        double r2 = delta.dot(delta);
        double r = sqrt(r2);
        delta *= 1/r;
        double dr = r-distance;
        double kdr = k*dr;
        energy += 0.5*kdr*dr;
        g[3*particle1] -= kdr*delta[0];
        g[3*particle1+1] -= kdr*delta[1];
        g[3*particle1+2] -= kdr*delta[2];
        g[3*particle2] += kdr*delta[0];
        g[3*particle2+1] += kdr*delta[1];
        g[3*particle2+2] += kdr*delta[2];
    }
    return energy;
}

int isolde::LocalEnergyMinimizer::minimize(Context& context, double tolerance,
        int maxIterations) {
    int ret = 0;
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    double constraintTol = context.getIntegrator().getConstraintTolerance();
    double workingConstraintTol = std::max(1e-4, constraintTol);
    double k = 100/workingConstraintTol;
    lbfgsfloatval_t *x = lbfgs_malloc(numParticles*3);
    if (x == NULL)
        throw OpenMMException("LocalEnergyMinimizer: Failed to allocate memory");
    try {

        // Initialize the minimizer.

        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);
        if (!context.getPlatform().supportsDoublePrecision())
            param.xtol = 1e-7;
        param.max_iterations = maxIterations;
        param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;

        // Make sure the initial configuration satisfies all constraints.

        context.applyConstraints(workingConstraintTol);

        // Record the initial positions and determine a normalization constant for scaling the tolerance.

        std::vector<Vec3> initialPos = context.getState(State::Positions).getPositions();
        double norm = 0.0;
        for (int i = 0; i < numParticles; i++) {
            const Vec3& p = initialPos[i];
            x[3*i] = p[0];
            x[3*i+1] = p[1];
            x[3*i+2] = p[2];
            norm += p.dot(p);
        }
        norm /= numParticles;
        norm = (norm < 1 ? 1 : sqrt(norm));
        param.epsilon = tolerance/norm;

        // Appears to improve the stability of the minimization, without too
        // much impact on performance.
        param.m = 15;
        param.ftol = 1e-6;

        /* In tests repeatedly minimizing from coordinates that have been
         * randomly perturbed by up to +/-0.4A along each axis, this leads to
         * calls to lbfgs() exiting early (returning LBFGSERR_MAXIMUMSTEP) about
         * 1 time in 10. That seems fairly harmless, whereas higher values
         * occasionally allow the minimizer to stray into wildly wrong
         * configurations that it can struggle to recover from.
        */
        param.max_step = 1.0;


        // Repeatedly minimize, steadily increasing the strength of the springs until all constraints are satisfied.

        double prevMaxError = 1e10;
        MinimizerData data(context, k);
        int tries=0, maxTries=3;
        while (true) {
            // Perform the minimization.

            lbfgsfloatval_t fx;
            int result = -1024;
            try {
                result = lbfgs(numParticles*3, x, &fx, evaluate, NULL, &data, &param);
            } catch (std::out_of_range) {
                ret = INFINITE_OR_NAN_FORCE;
                break;
            }
            if (LBFGSERR_MAXIMUMSTEP == result)
            {
                tries++;
                if (tries < maxTries)
                    // Try a couple more times to see if resetting the
                    // minimizer helps find a better path from the current
                    // coordinates.
                    continue;
                else
                    // Continue on, and (where applicable) try again with
                    // tighter spring constants on constrained bonds.
                    tries=0;
            }

            // Check whether all constraints are satisfied.

            std::vector<Vec3> positions = context.getState(State::Positions).getPositions();
            int numConstraints = system.getNumConstraints();
            double maxError = 0.0;
            for (int i = 0; i < numConstraints; i++) {
                int particle1, particle2;
                double distance;
                system.getConstraintParameters(i, particle1, particle2, distance);
                Vec3 delta = positions[particle2]-positions[particle1];
                double r = sqrt(delta.dot(delta));
                double error = fabs(r-distance);
                if (error > maxError)
                    maxError = error;
            }
            if (maxError <= workingConstraintTol)
            {
                if (LBFGSERR_MAXIMUMITERATION == result)
                    ret = DID_NOT_CONVERGE;
                break; // All constraints are satisfied.
            }
            context.setPositions(initialPos);
            if (maxError >= prevMaxError)
            {
                if (data.largeForceEncountered)
                {
                    ret = CONSTRAINT_VIOLATION_LARGE_FORCE;
                } else {
                    ret = CONSTRAINT_VIOLATION_NO_LARGE_FORCE;
                }
                break; // Further tightening the springs doesn't seem to be helping, so just give up.
            }
            prevMaxError = maxError;
            data.k *= 10;
            if (maxError > 100*workingConstraintTol) {
                // We've gotten far enough from a valid state that we might have trouble getting
                // back, so reset to the original positions.

                for (int i = 0; i < numParticles; i++) {
                    x[3*i] = initialPos[i][0];
                    x[3*i+1] = initialPos[i][1];
                    x[3*i+2] = initialPos[i][2];
                }
            }
        }
    }
    catch (...) {
        lbfgs_free(x);
        throw;
        return UNSPECIFIED_ERROR;
    }
    lbfgs_free(x);
    return ret;

    // If necessary, do a final constraint projection to make sure they are satisfied
    // to the full precision requested by the user.

    if (constraintTol < workingConstraintTol)
        context.applyConstraints(workingConstraintTol);
}
