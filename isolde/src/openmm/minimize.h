/**
 * @Author: Tristan Croll <tic20>
 * @Date:   11-Jun-2019
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 11-Jun-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#pragma once

/*  This is mostly a copy of OpenMM's own L-BFGS minimiser, with a few tweaks to
 *  suit ISOLDE's purposes. Most importantly, when using a GPU platform it
 *  avoids falling back to the CPU at all costs, since this blows out timeframes
 *  from seconds to many minutes in the presence of MDFF. Instead, if it
 *  encounters out-of-gamut forces **and** fails to successfully converge it
 *  will simply bail out. ISOLDE's existing code will find the offending atoms
 *  and highlight them for the user.
 */

#include <OpenMM.h>

namespace isolde
{

class LocalEnergyMinimizer
{

public:
    /**
     * Search for a new set of particle positions that represent a local potential energy minimum.
     * On exit, the Context will have been updated with the new positions.
     *
     * @param context        a Context specifying the System to minimize and the initial particle positions
     * @param tolerance      this specifies how precisely the energy minimum must be located.  Minimization
     *                       will be halted once the root-mean-square value of all force components reaches
     *                       this tolerance.  The default value is 10.
     * @param maxIterations  the maximum number of iterations to perform.  If this is 0, minimation is continued
     *                       until the results converge without regard to how many iterations it takes.  The
     *                       default value is 0.
     */
    static int minimize(OpenMM::Context& context, double tolerance = 10, int maxIterations = 0);
    enum {
        UNSPECIFIED_ERROR = -1024,
        INFINITE_OR_NAN_FORCE = -3,
        CONSTRAINT_VIOLATION_NO_LARGE_FORCE = -2,
        CONSTRAINT_VIOLATION_LARGE_FORCE = -1,
        SUCCESS = 0,
        DID_NOT_CONVERGE=1
    };

};

} // namespace isolde
