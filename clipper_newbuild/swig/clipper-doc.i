
// File: index.xml

// File: classclipper_1_1datatypes_1_1ABCD.xml


%feature("docstring") clipper::datatypes::ABCD "

Reflection data type: Hendrickson-Lattman coeff.  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::ABCD::ABCD "
";

%feature("docstring") clipper::datatypes::ABCD::ABCD "
";

%feature("docstring") clipper::datatypes::ABCD::d "
";

%feature("docstring") clipper::datatypes::ABCD::d "
";

%feature("docstring") clipper::datatypes::ABCD::c "
";

%feature("docstring") clipper::datatypes::ABCD::c "
";

%feature("docstring") clipper::datatypes::ABCD::a "
";

%feature("docstring") clipper::datatypes::ABCD::a "
";

%feature("docstring") clipper::datatypes::ABCD::friedel "
";

%feature("docstring") clipper::datatypes::ABCD::data_export "
";

%feature("docstring") clipper::datatypes::ABCD::b "
";

%feature("docstring") clipper::datatypes::ABCD::b "
";

%feature("docstring") clipper::datatypes::ABCD::shift_phase "
";

%feature("docstring") clipper::datatypes::ABCD::type "
";

%feature("docstring") clipper::datatypes::ABCD::data_size "
";

%feature("docstring") clipper::datatypes::ABCD::missing "
";

%feature("docstring") clipper::datatypes::ABCD::set_null "
";

%feature("docstring") clipper::datatypes::ABCD::data_import "
";

%feature("docstring") clipper::datatypes::ABCD::data_names "
";

// File: classclipper_1_1Array2d.xml


%feature("docstring") clipper::Array2d "

Simple 2-d array class.  

C++ includes: clipper_types.h
";

%feature("docstring") clipper::Array2d::resize "

resize  
";

%feature("docstring") clipper::Array2d::resize "

resize  
";

%feature("docstring") clipper::Array2d::cols "

number of cols  
";

%feature("docstring") clipper::Array2d::Array2d "

null constructor  
";

%feature("docstring") clipper::Array2d::Array2d "

constructor  
";

%feature("docstring") clipper::Array2d::Array2d "

constructor  
";

%feature("docstring") clipper::Array2d::rows "

number of rows  
";

%feature("docstring") clipper::Array2d::size "

size  
";

// File: classclipper_1_1Atom.xml


%feature("docstring") clipper::Atom "

Atom class.  

This class defines a minimal atom object providing only those properties
required for an electron density calculation. A template constructor allows it
to be constructed from any other object with appropriate properties.  

C++ includes: coords.h
";

%feature("docstring") clipper::Atom::set_u_aniso_orth "

set u_aniso  
";

%feature("docstring") clipper::Atom::element "

get atom element name: e.g. \"C\", \"N\", \"Zn2+\"  
";

%feature("docstring") clipper::Atom::set_element "

set element  
";

%feature("docstring") clipper::Atom::coord_orth "

get atom orthogonal (Angstrom) coordinate  
";

%feature("docstring") clipper::Atom::set_coord_orth "

set coord_orth  
";

%feature("docstring") clipper::Atom::u_iso "

get atom orthogonal isotropic U value  
";

%feature("docstring") clipper::Atom::transform "

apply a rotation-translation operator (RTop) to the atom  

The coordinates and U_aniso_orth are transformed. The sigmas are not, since
without the full variance-covariance matrix this transformation is impossible.  

Parameters
----------
* `rt` :  
    The operator to apply.  
";

%feature("docstring") clipper::Atom::set_occupancy "

set occupancy  
";

%feature("docstring") clipper::Atom::u_aniso_orth "

get atom orthogonal anisotropic U value  
";

%feature("docstring") clipper::Atom::set_u_iso "

set u_iso  
";

%feature("docstring") clipper::Atom::null "

return null atom  
";

%feature("docstring") clipper::Atom::Atom "

null constructor  
";

%feature("docstring") clipper::Atom::Atom "

Constructor: from atom-like object.  
";

%feature("docstring") clipper::Atom::occupancy "

get atom occupancy  
";

%feature("docstring") clipper::Atom::is_null "

test for null atom: atom is null is coord is null  
";

// File: classclipper_1_1Atom__list.xml


%feature("docstring") clipper::Atom_list "

Atom list class.  

This class defines a minimal atom list object providing only those properties
required for an electron density calculation. It is a trivial derivation from
std::vector<Atom>. In addition a template constructor allows it to be
constructed from any other object with appropriate properties.  

C++ includes: coords.h
";

%feature("docstring") clipper::Atom_list::Atom_list "

null constructor  
";

%feature("docstring") clipper::Atom_list::Atom_list "

constructor to generate a list of empty Atom objects to be filled later  
";

%feature("docstring") clipper::Atom_list::Atom_list "

constructor: from std::vector<Atom>  
";

%feature("docstring") clipper::Atom_list::Atom_list "

Constructor: from vector-like list of atom-like objects.  
";

// File: classclipper_1_1AtomSF.xml


%feature("docstring") clipper::AtomSF "

Atomic scattering factor object.  

Deprecated
This class has been replaced by AtomShapeFn, which is smaller, faster, and more
capable. This class is now a wrapper for that class.  

C++ includes: atomsf.h
";

%feature("docstring") clipper::AtomSF::rho_iso "
";

%feature("docstring") clipper::AtomSF::rho_aniso "
";

%feature("docstring") clipper::AtomSF::AtomSF "
";

%feature("docstring") clipper::AtomSF::AtomSF "
";

%feature("docstring") clipper::AtomSF::f_aniso "
";

%feature("docstring") clipper::AtomSF::init "
";

%feature("docstring") clipper::AtomSF::init "
";

%feature("docstring") clipper::AtomSF::f_iso "
";

// File: classclipper_1_1AtomShapeFn.xml


%feature("docstring") clipper::AtomShapeFn "

Atomic shape function object.  

The atomic scattering factor object is instantiated for each atom in turn,
giving the atom parameters: position, element, occupancy and the isotropic or
anisotropic U-value. (See clipper::Util for conversion from B-factors.). The
methods of the class may then be called to return the scattering in reciprocal
space or density in real space using either isotropic or anistropic models as
required.  

If the atom only has an isotropic U, the faster isotropic methods will be used
where available.  

This implementation uses the coefficients from Waasmaier & Kirfel (1995), Acta
Cryst. A51, 416-431. The source data can be found at: ftp://wrzx02.rz.uni-
wuerzburg.de/pub/local/Crystallography/sfac.dat  

C++ includes: atomsf.h
";

%feature("docstring") clipper::AtomShapeFn::rho_grad "

return Agarwal density gradients as a function of coordinate  

Return the Agarwal gradients of the density with respect to tha atomic
parameters as a function of position in real space in electrons. The parameter
list is defined by assignment to agarwal_params().  

Parameters
----------
* `xyz` :  
    Position in real space.  rho The density in electrons.  
* `grad` :  
    Vector gradient in electrons (pre-size for best performance).  
";

%feature("docstring") clipper::AtomShapeFn::rho_grad "

Deprecated
return Agarwal density gradients as a function of coordinate  
";

%feature("docstring") clipper::AtomShapeFn::AtomShapeFn "

null constructor  
";

%feature("docstring") clipper::AtomShapeFn::AtomShapeFn "

constructor: from atom object  

If the atom has an anisotropic U (even if the values are isotropic), then it is
initialised as anisotropic, otherwise it is isotropic.  

Parameters
----------
* `atom` :  
    The atom object.  
";

%feature("docstring") clipper::AtomShapeFn::AtomShapeFn "

constructor: from coord, element, isotropic U, occupancy  

The atom is initialised as isotropic.  

Parameters
----------
* `xyz` :  
    The atom coordinate.  
* `element` :  
    The atom element.  
* `u_iso` :  
    The isotropic U-value.  
* `occ` :  
    The occupancy.  
";

%feature("docstring") clipper::AtomShapeFn::AtomShapeFn "

constructor: from coord, element, anisotropic U, occupancy  

The atom is initialised as anisotropic.  

Parameters
----------
* `xyz` :  
    The atom coordinate.  
* `element` :  
    The atom element.  
* `u_aniso` :  
    The anisotropic U-value.  
* `occ` :  
    The occupancy.  
";

%feature("docstring") clipper::AtomShapeFn::f "

return scattering factor as a function of reflection posn  

Return the scattering factor as a function of position in reciprocal space in
electrons.  

Parameters
----------
* `rfl` :  
    Position in reciprocal space.  

Returns
-------
The scattering factor in electrons.  
";

%feature("docstring") clipper::AtomShapeFn::f "

return (isotropic) scattering factor as a function of resolution  

Return the scattering factor as a function of position in reciprocal space in
electrons.  

Parameters
----------
* `invresolsq` :  
    Inverse resolution squared in inverse Angstroms squared.  

Returns
-------
The scattering factor in electrons.  
";

%feature("docstring") clipper::AtomShapeFn::agarwal_params "

define parameters for Agarwal gradient/curvature calcs  
";

%feature("docstring") clipper::AtomShapeFn::init "

initialiser: from atom object  

If the atom has an anisotropic U (even if the values are isotropic), then it is
initialised as anisotropic, otherwise it is isotropic.  

Parameters
----------
* `atom` :  
    The atom object.  
";

%feature("docstring") clipper::AtomShapeFn::init "

initialiser: from coord, element, isotropic U, occupancy  

The atom is initialised as isotropic.  

Parameters
----------
* `xyz` :  
    The atom coordinate.  
* `element` :  
    The atom element.  
* `u_iso` :  
    The isotropic U-value.  
* `occ` :  
    The occupancy.  
";

%feature("docstring") clipper::AtomShapeFn::init "

initialiser: from coord, element, anisotropic U, occupancy  

The atom is initialised as anisotropic.  

Parameters
----------
* `xyz` :  
    The atom coordinate.  
* `element` :  
    The atom element.  
* `u_aniso` :  
    The anisotropic U-value.  
* `occ` :  
    The occupancy.  
";

%feature("docstring") clipper::AtomShapeFn::rho_curv "

return Agarwal density gradient/curvature as a function of coordinate  

Return the Agarwal gradients of the density with respect to tha atomic
parameters as a function of position in real space in electrons. The parameter
list is defined by assignment to agarwal_params().  

Parameters
----------
* `xyz` :  
    Position in real space.  rho The density in electrons.  
* `grad` :  
    Vector gradient in electrons (pre-size for best performance).  
* `curv` :  
    Matrix curvature in electrons (pre-size for best performance).  
";

%feature("docstring") clipper::AtomShapeFn::rho "

return electron density as a function of coordinate  

Return the density as a function of position in real space in electrons.  

Parameters
----------
* `xyz` :  
    Position in real space.  
* `The` :  
    density in electrons.  
";

%feature("docstring") clipper::AtomShapeFn::rho "

return (isotropic) electron density as a function of radius  

Return the density as a function of position in real space in electrons.  

Parameters
----------
* `invresolsq` :  
    Radius squared in Angstroms squared.  
* `The` :  
    density in electrons.  
";

// File: classclipper_1_1BasisFn__aniso__gaussian.xml


%feature("docstring") clipper::BasisFn_aniso_gaussian "

simple anisotropic Gaussian basis function  

This class provides a anisotropic Gaussian basis function.  

C++ includes: resol_basisfn.h
";

%feature("docstring") clipper::BasisFn_aniso_gaussian::BasisFn_aniso_gaussian "

constructor:  
";

%feature("docstring") clipper::BasisFn_aniso_gaussian::fderiv_coord "

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_aniso_gaussian::scale "

return the scale factor corresponding to the Gaussian parameters  
";

%feature("docstring") clipper::BasisFn_aniso_gaussian::u_aniso_orth "

return the anisotropic U corresponding to the Gaussian parameters  
";

%feature("docstring") clipper::BasisFn_aniso_gaussian::fderiv "

the value of the resolution function (override for speed)  

the derivatives of the resolution function w.r.t. the parameters  
";

// File: classclipper_1_1BasisFn__base.xml


%feature("docstring") clipper::BasisFn_base "

abstract base class for resolution function basis functions  

A basis function must be able to return its value and derivatives for any given
HKL.  

Optionally, performance can be improved by returning a flag to indicate if the
value of the basis function for a given reflection is linearly dependent on the
values of the parameter, and a value indicating whether the curvature matrix
takes an N-diagonal form.  

**NOTE:** for performance reasons the derivatives are returned as a reference to
an internal object, so if you store a reference to the result (which is also
good for performance, it will be overwritten on the next call. If this matters,
store a copy rather than a reference.  

C++ includes: resol_fn.h
";

%feature("docstring") clipper::BasisFn_base::num_params "

the number of parameters of this basis function  
";

%feature("docstring") clipper::BasisFn_base::f "

the value of the resolution function  
";

%feature("docstring") clipper::BasisFn_base::BasisFn_base "

null constructor  
";

%feature("docstring") clipper::BasisFn_base::BasisFn_base "

constructor: takes number of parameters  
";

%feature("docstring") clipper::BasisFn_base::type "

the type of the function: optionally used to improve convergence  

Defaults to GENERAL, which will always work. If the basis function is linearly
dependent on the parameters, override this with a function returning LINEAR for
improved performance. See the provided basis functions for examples.  

Returns
-------
The function type enumeration.  
";

%feature("docstring") clipper::BasisFn_base::fderiv "

the value of the resolution function and its first two derivatives  
";

%feature("docstring") clipper::BasisFn_base::num_diagonals "

number of non-zero diagonals in the upper triangle of the curvatures  

Defaults to 0, which will always work. If the basis function has compact support
among the parameters, i.e. the value for any HKL depends only on a few
parameters, then set this to the number of non-zero diagonals in the upper
triangle of the matrix, i.e. 1 for a diagonal matrix, 2 for a tri-diagonal
matrix etc.  

Returns
-------
The number of non-zero upper diagonals, or zero for full-matrix.  
";

// File: classclipper_1_1BasisFn__binner.xml


%feature("docstring") clipper::BasisFn_binner "

simple binning basis function  

This class bins reflections on the basis of resolution, i.e. it generates a
resolution function from spherical shells.  

C++ includes: resol_basisfn.h
";

%feature("docstring") clipper::BasisFn_binner::f_s "

the value of the resolution function (override for speed)  
";

%feature("docstring") clipper::BasisFn_binner::fderiv_s "

the derivative of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_binner::f "

the value of the resolution function (override for speed)  
";

%feature("docstring") clipper::BasisFn_binner::fderiv "

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_binner::type "

the type of the function: optionally used to improve convergence  
";

%feature("docstring") clipper::BasisFn_binner::BasisFn_binner "

constructor: include whole reflection list in histogram  
";

%feature("docstring") clipper::BasisFn_binner::BasisFn_binner "

constructor: include only non-missing reflections in histogram  
";

%feature("docstring") clipper::BasisFn_binner::num_diagonals "

number of non-zero diagonals in the upper triangle of the curvatures  
";

// File: classclipper_1_1BasisFn__expcubic.xml


%feature("docstring") clipper::BasisFn_expcubic "

simple Expcubic basis function  

This class provides a Expcubic basis function.  

C++ includes: resol_basisfn.h
";

%feature("docstring") clipper::BasisFn_expcubic::fderiv_s "

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_expcubic::fderiv "

the value of the resolution function (override for speed)  

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_expcubic::BasisFn_expcubic "

constructor  
";

// File: classclipper_1_1BasisFn__gaussian.xml


%feature("docstring") clipper::BasisFn_gaussian "

simple Gaussian basis function  

This class provides a Gaussian basis function.  

C++ includes: resol_basisfn.h
";

%feature("docstring") clipper::BasisFn_gaussian::fderiv "

the value of the resolution function (override for speed)  

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_gaussian::fderiv_s "

the value of the resolution function  

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_gaussian::u_iso "

return the isotropic U corresponding to the Gaussian parameters  
";

%feature("docstring") clipper::BasisFn_gaussian::scale "

return the scale factor corresponding to the Gaussian parameters  
";

%feature("docstring") clipper::BasisFn_gaussian::BasisFn_gaussian "

constructor:  
";

// File: classclipper_1_1BasisFn__linear.xml


%feature("docstring") clipper::BasisFn_linear "

simple linear basis function  

This class fits a piecewise linear function through reflections on the basis of
resolution.  

C++ includes: resol_basisfn.h
";

%feature("docstring") clipper::BasisFn_linear::fderiv_s "

the derivative of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_linear::f "

the value of the resolution function (override for speed)  
";

%feature("docstring") clipper::BasisFn_linear::fderiv "

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_linear::num_diagonals "

number of non-zero diagonals in the upper triangle of the curvatures  
";

%feature("docstring") clipper::BasisFn_linear::type "

the type of the function: optionally used to improve convergence  
";

%feature("docstring") clipper::BasisFn_linear::f_s "

the value of the resolution function (override for speed)  
";

%feature("docstring") clipper::BasisFn_linear::BasisFn_linear "

constructor: include whole reflection list in histogram  
";

%feature("docstring") clipper::BasisFn_linear::BasisFn_linear "

constructor: include only non-missing reflections in histogram  
";

// File: classclipper_1_1BasisFn__log__aniso__gaussian.xml


%feature("docstring") clipper::BasisFn_log_aniso_gaussian "

simple anisotropic Gaussian basis function  

This class provides a anisotropic Gaussian basis function. i.e. a general
quadratic function of resolution. Use this in conjunction with a Log-target
function to get a fast estimate to a Gaussian fit.  

C++ includes: resol_basisfn.h
";

%feature("docstring") clipper::BasisFn_log_aniso_gaussian::type "

the type of the function: optionally used to improve convergence  
";

%feature("docstring") clipper::BasisFn_log_aniso_gaussian::scale "

return the scale factor corresponding to the Gaussian parameters  
";

%feature("docstring") clipper::BasisFn_log_aniso_gaussian::fderiv "

the value of the resolution function (override for speed)  

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_log_aniso_gaussian::BasisFn_log_aniso_gaussian "

constructor:  
";

%feature("docstring") clipper::BasisFn_log_aniso_gaussian::u_aniso_orth "

return the anisotropic U corresponding to the Gaussian parameters  
";

%feature("docstring") clipper::BasisFn_log_aniso_gaussian::fderiv_coord "

the derivatives of the resolution function w.r.t. the parameters  
";

// File: classclipper_1_1BasisFn__log__gaussian.xml


%feature("docstring") clipper::BasisFn_log_gaussian "

simple log Gaussian basis function  

This class provides a Log Gaussian basis function. i.e. a quadratic function of
resolution. Use this in conjunction with a Log-target function to get a fast
estimate to a Gaussian fit.  

C++ includes: resol_basisfn.h
";

%feature("docstring") clipper::BasisFn_log_gaussian::fderiv "

the value of the resolution function (override for speed)  

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_log_gaussian::fderiv_s "

the value of the resolution function  

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_log_gaussian::u_iso "

return the isotropic U corresponding to the Gaussian parameters  
";

%feature("docstring") clipper::BasisFn_log_gaussian::BasisFn_log_gaussian "

constructor:  
";

%feature("docstring") clipper::BasisFn_log_gaussian::type "

the type of the function: optionally used to improve convergence  
";

%feature("docstring") clipper::BasisFn_log_gaussian::scale "

return the scale factor corresponding to the Gaussian parameters  
";

// File: classclipper_1_1BasisFn__spline.xml


%feature("docstring") clipper::BasisFn_spline "

simple smooth basis function  

This class fits a Bspline through reflections on the basis of resolution.  

C++ includes: resol_basisfn.h
";

%feature("docstring") clipper::BasisFn_spline::fderiv_s "

the derivative of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_spline::f "

the value of the resolution function (override for speed)  
";

%feature("docstring") clipper::BasisFn_spline::fderiv "

the derivatives of the resolution function w.r.t. the parameters  
";

%feature("docstring") clipper::BasisFn_spline::num_diagonals "

number of non-zero diagonals in the upper triangle of the curvatures  
";

%feature("docstring") clipper::BasisFn_spline::type "

the type of the function: optionally used to improve convergence  
";

%feature("docstring") clipper::BasisFn_spline::f_s "

the value of the resolution function (override for speed)  
";

%feature("docstring") clipper::BasisFn_spline::BasisFn_spline "

constructor: include whole reflection list in histogram  
";

%feature("docstring") clipper::BasisFn_spline::BasisFn_spline "

constructor: include only non-missing reflections in histogram  
";

// File: classclipper_1_1CCell.xml


%feature("docstring") clipper::CCell "

CCell container.  

CCell: This has a name and a cell. It overrides the base cell for any objects
below it.  

C++ includes: container_types.h
";

%feature("docstring") clipper::CCell::init "

initialiser: from Cell  

The object is initialised, and children are updated.  

Parameters
----------
* `cell_` :  
    The value to give to the contained object.  
";

%feature("docstring") clipper::CCell::CCell "

constructor: make null object or top object in tree  
";

%feature("docstring") clipper::CCell::CCell "

constructor: make child object  
";

// File: classclipper_1_1Cell.xml


%feature("docstring") clipper::Cell "

Cell object.  

The Cell class is the fully functional description of the unit cell. In addition
to the cell parameters, it stores derived information including the cell volume,
orthogonalising and fractionalising matrices, and the metric tensors.  

C++ includes: cell.h
";

%feature("docstring") clipper::Cell::beta_deg "
";

%feature("docstring") clipper::Cell::gamma_star "

get gamma*  
";

%feature("docstring") clipper::Cell::equals "

test equality with another cell  

Two cells disagree if the difference in their orthogonalisation matrices is
sufficient to map a reflection from one cell onto a different reflection in the
other cell at the given tolerance, which is the resolution of the reflection in
Angstroms.  

Parameters
----------
* `other` :  
    The other cell to compare.  
* `tol` :  
    The tolerance, in Angstroms.  
";

%feature("docstring") clipper::Cell::a_star "

get a*  
";

%feature("docstring") clipper::Cell::alpha "
";

%feature("docstring") clipper::Cell::beta_star "

get beta*  
";

%feature("docstring") clipper::Cell::is_null "

test if object has been initialised  

Returns
-------
true if the object has not been initalised.  
";

%feature("docstring") clipper::Cell::Cell "

null constructor: must initialise later  
";

%feature("docstring") clipper::Cell::Cell "

constructor: takes a Cell descriptor  
";

%feature("docstring") clipper::Cell::metric_real "

return real space metric tensor  
";

%feature("docstring") clipper::Cell::b_star "

get b*  
";

%feature("docstring") clipper::Cell::format "
";

%feature("docstring") clipper::Cell::init "

initialiser  

Initialise the Cell object from a cell description.  

Parameters
----------
* `cell_` :  
    The cell descirption.  
";

%feature("docstring") clipper::Cell::descr "

return cell dimensions  
";

%feature("docstring") clipper::Cell::alpha_deg "
";

%feature("docstring") clipper::Cell::a "
";

%feature("docstring") clipper::Cell::c "
";

%feature("docstring") clipper::Cell::b "
";

%feature("docstring") clipper::Cell::gamma "
";

%feature("docstring") clipper::Cell::debug "
";

%feature("docstring") clipper::Cell::beta "
";

%feature("docstring") clipper::Cell::volume "

return cell volume  
";

%feature("docstring") clipper::Cell::gamma_deg "
";

%feature("docstring") clipper::Cell::alpha_star "

get alpha*  
";

%feature("docstring") clipper::Cell::c_star "

get c*  
";

%feature("docstring") clipper::Cell::matrix_frac "

return fractionalisation matrix  
";

%feature("docstring") clipper::Cell::metric_reci "

return reciprocal space metric tensor  
";

%feature("docstring") clipper::Cell::matrix_orth "

return orthogonalisation matrix  
";

// File: classclipper_1_1Cell__descr.xml


%feature("docstring") clipper::Cell_descr "

cell description (automatically converts to radians)  

The cell description is a compact description of a cell, containing just the
cell parameters. It is usually used to construct a full Cell object, which
provides the expected functionality.  

C++ includes: cell.h
";

%feature("docstring") clipper::Cell_descr::a "

get a  
";

%feature("docstring") clipper::Cell_descr::c "

get c  
";

%feature("docstring") clipper::Cell_descr::b "

get b  
";

%feature("docstring") clipper::Cell_descr::alpha "

get alpha  
";

%feature("docstring") clipper::Cell_descr::alpha_deg "

get alpha in degrees  

Returns
-------
The cell angle in degrees  
";

%feature("docstring") clipper::Cell_descr::gamma_deg "

get gamma in degrees  

Returns
-------
The cell angle in degrees  
";

%feature("docstring") clipper::Cell_descr::Cell_descr "

null constructor  
";

%feature("docstring") clipper::Cell_descr::Cell_descr "

constructor: from cell parameters  

Parameters
----------
* `a` :  
    A axis in Angstroms.  
* `b` :  
    B axis in Angstroms.  
* `c` :  
    C axis in Angstroms.  
* `alpha` :  
    Angle between B and C axes in radians or degrees, default=90  
* `beta` :  
    Angle between A and C axes in radians or degrees, default=90  
* `gamma` :  
    Angle between A and C axes in radians or degrees, default=90  
";

%feature("docstring") clipper::Cell_descr::beta_deg "

get alpha in degrees  

Returns
-------
The cell angle in degrees  
";

%feature("docstring") clipper::Cell_descr::beta "

get beta  
";

%feature("docstring") clipper::Cell_descr::format "

return formatted String representation  

Returns
-------
A string describing the cell  
";

%feature("docstring") clipper::Cell_descr::gamma "

get gamma  
";

// File: classclipper_1_1CGrid__sampling.xml


%feature("docstring") clipper::CGrid_sampling "

CGrid_sampling container.  

CGrid_sampling: This has a name and a grid sampling It overrides the grid
sampling for any objects below it.  

C++ includes: container_types.h
";

%feature("docstring") clipper::CGrid_sampling::CGrid_sampling "

constructor: make null object or top object in tree  

The top object in a tree is initialised from a known grid.  

Parameters
----------
* `name` :  
    The object name.  
* `grid` :  
    The grid sampling.  
";

%feature("docstring") clipper::CGrid_sampling::CGrid_sampling "

constructor: make child object  

The normal form for a child object - spacegroup and cell inherited.  

Parameters
----------
* `parent` :  
    The objects parent.  
* `name` :  
    The object name.  
* `rate` :  
    The Shannon rate (default 1.5).  
";

%feature("docstring") clipper::CGrid_sampling::CGrid_sampling "

constructor: make child object with explicit value  

This is still a child object but is initialised directly.  

Parameters
----------
* `parent` :  
    The objects parent.  
* `name` :  
    The object name.  
* `grid` :  
    The grid sampling.  
";

%feature("docstring") clipper::CGrid_sampling::init "

initialiser: from sampling rate  

The object is initialised if the appropriate parent objects are available, and
children are updated.  

Parameters
----------
* `spacegroup` :  
    The spacegroup.  
* `cell` :  
    The cell.  
* `resolution` :  
    The resolution.  
* `rate_` :  
    The Shannon rate (If <1 previous value is used, default 1.5).  
";

%feature("docstring") clipper::CGrid_sampling::init "

initialiser: from Grid_sampling  

The object is initialised, and children are updated.  

Parameters
----------
* `grid_sampling_` :  
    The value to give to the contained object.  
";

%feature("docstring") clipper::CGrid_sampling::update "

hierarchical update  

Hierarchical update. If this object is uninitialised, an attempt is made to
initialise the object using information from its parents in the hierarchy. The
childen of the object are then updated.  
";

// File: classclipper_1_1CHKL__data.xml


%feature("docstring") clipper::CHKL_data "

Reflection data list container.  

CHKL_data: This is the list object containing the actual data. It must be
indexed by a parent HKL list.  

C++ includes: container_hkl.h
";

%feature("docstring") clipper::CHKL_data::CHKL_data "

null constructor  
";

%feature("docstring") clipper::CHKL_data::CHKL_data "

constructor: inherit datalist and cell  

The object is constructed at the given location in the hierarchy. An attempt is
made to initialise the object using information from its parents in the
hierarchy.  

Parameters
----------
* `parent` :  
    An object in the hierarchy (usually the parent of the new object).  
* `name` :  
    The path from `parent` to the new object (usually just the name of the new
    object).  
";

%feature("docstring") clipper::CHKL_data::init "

initialiser: supply or inherit hkl list, and cell  

An attempt is made to initialise the object using information from the supplied
parameters, or if they are Null, from its parents in the hierarchy.  

Parameters
----------
* `hkl_info` :  
    The reflection list object for this datalist.  
* `cell` :  
    The cell object for this datalist.  
";

%feature("docstring") clipper::CHKL_data::init "

initialiser: from spacegroup, cell, and HKL_sampling  
";

%feature("docstring") clipper::CHKL_data::update "

hierarchical update  

Hierarchical update. If this object is uninitialised, an attempt is made to
initialise the object using information from its parents in the hierarchy. The
childen of the object are then updated.  

The data list is also synchronized with the parent reflection list.  
";

// File: classclipper_1_1CHKL__info.xml


%feature("docstring") clipper::CHKL_info "

HKL list and indexing object container.  

CHKL_info: This is the reflection list object for reflection data objects to
reference in order to identify their data entries.  

C++ includes: container_hkl.h
";

%feature("docstring") clipper::CHKL_info::init "

initialiser: supply or inherit spacegroup, cell and resolution  

An attempt is made to initialise the object using information from the supplied
parameters, or if they are Null, from its parents in the hierarchy.  

Parameters
----------
* `spacegroup` :  
    The spacegroup.  
* `cell` :  
    The cell.  
* `resolution` :  
    The resolution.  
* `generate` :  
    Generate reflection list if true.  
";

%feature("docstring") clipper::CHKL_info::update "

hierarchical update  

Hierarchical update. If this object is uninitialised, an attempt is made to
initialise the object using information from its parents in the hierarchy. The
childen of the object are then updated.  
";

%feature("docstring") clipper::CHKL_info::CHKL_info "

constructor: make null object or top object in tree  

Construct and initialise as the top object in a tree.  

Parameters
----------
* `spacegroup` :  
    The spacegroup.  
* `cell` :  
    The cell.  
* `resolution` :  
    The resolution.  
* `generate` :  
    Generate reflection list if true.  
";

%feature("docstring") clipper::CHKL_info::CHKL_info "

constructor: inherit spacegroup, cell and resolution  

The object is constructed at the given location in the hierarchy. An attempt is
made to initialise the object using information from its parents in the
hierarchy.  

Parameters
----------
* `parent` :  
    An object in the hierarchy (usually the parent of the new object).  
* `name` :  
    The path from `parent` to the new object (usually just the name of the new
    object).  
";

%feature("docstring") clipper::CHKL_info::generate_hkl_list "

synthesize hkl list and update children  

The reflection list is sythesized to match the given spacegroup, cell, and
resolution, and a hierarchical update is triggered to update the sizes of the
reflection lists for all dependent HKL_data objects.  
";

// File: classclipper_1_1ClipperInstance.xml


%feature("docstring") clipper::ClipperInstance "
";

%feature("docstring") clipper::ClipperInstance::destroy "

VERY DANGEROUS, DO NOT USE.  
";

%feature("docstring") clipper::ClipperInstance::~ClipperInstance "
";

%feature("docstring") clipper::ClipperInstance::hkl_data_cache "
";

%feature("docstring") clipper::ClipperInstance::ClipperInstance "
";

%feature("docstring") clipper::ClipperInstance::xmap_cache "
";

%feature("docstring") clipper::ClipperInstance::spacegroup_cache "
";

%feature("docstring") clipper::ClipperInstance::util "
";

// File: classclipper_1_1ClipperInstantiator.xml


%feature("docstring") clipper::ClipperInstantiator "
";

%feature("docstring") clipper::ClipperInstantiator::~ClipperInstantiator "
";

%feature("docstring") clipper::ClipperInstantiator::ClipperInstantiator "
";

%feature("docstring") clipper::ClipperInstantiator::instance "
";

// File: classclipper_1_1CNXmap.xml


%feature("docstring") clipper::CNXmap "

Non-Crystallographic map container.  

CNXmap: This is a non-crystallographic map. Since it does not exist in a
crystallographic frame, it does not inherit anything.  

C++ includes: container_map.h
";

%feature("docstring") clipper::CNXmap::CNXmap "

null constructor  
";

%feature("docstring") clipper::CNXmap::CNXmap "

constructor:  
";

// File: classclipper_1_1CNXmap__operator.xml


%feature("docstring") clipper::CNXmap_operator "

Non-Crystallographic map operator container.  

CNXmap: This is an operator relating a non-crystallographic map into a
crystallgraphic frame. It can inherit the crystallographic cell and grid
sampling.  

C++ includes: container_map.h
";

%feature("docstring") clipper::CNXmap_operator::update "

hierarchical update  

Hierarchical update. If this object is uninitialised, an attempt is made to
initialise the object using information from its parents in the hierarchy. The
childen of the object are then updated.  
";

%feature("docstring") clipper::CNXmap_operator::CNXmap_operator "

null constructor  
";

%feature("docstring") clipper::CNXmap_operator::CNXmap_operator "

constructor: do not initialise  

The object is not initialised.  

Parameters
----------
* `parent` :  
    The objects parent.  
* `name` :  
    The object name.  
";

%feature("docstring") clipper::CNXmap_operator::CNXmap_operator "

constructor: inherit cell and grid  

The object is initialised if the appropriate parent objects are available, and
children are updated.  

Parameters
----------
* `parent` :  
    The objects parent.  
* `name` :  
    The object name.  
* `nxmap` :  
    The non-crystal map object.  
* `nxop` :  
    The orth. operator mapping the NXmap into the crystal frame.  
";

%feature("docstring") clipper::CNXmap_operator::init "

initialier: supply or inherit cell, grid, NXmap, RTop_orth  

An attempt is made to initialise the object using information from the supplied
parameters, or if they are Null, from its parents in the hierarchy.  

Parameters
----------
* `cell` :  
    The unit cell for the crystallographic frame.  
* `grid` :  
    The grid sampling for the crystallographic frame.  
* `nxmap` :  
    The non-crystal map object.  
* `nxop` :  
    The orth. operator mapping the NXmap into the crystal frame.  
";

// File: classclipper_1_1Map__index__sort_1_1Compare__density.xml

// File: structclipper_1_1Compare__grid.xml


%feature("docstring") clipper::Compare_grid "
";

// File: classclipper_1_1datatypes_1_1Compute__abcd__from__phifom.xml


%feature("docstring") clipper::datatypes::Compute_abcd_from_phifom "

Compute from Phi_fom to ABCD ( C = D = 0 )  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1datatypes_1_1Compute__add__abcd.xml


%feature("docstring") clipper::datatypes::Compute_add_abcd "

Add two ABCD datalists.  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1datatypes_1_1Compute__add__fphi.xml


%feature("docstring") clipper::datatypes::Compute_add_fphi "

Add two F_phi datalists.  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1datatypes_1_1Compute__diff__fsigf__from__fsigfano.xml


%feature("docstring") clipper::datatypes::Compute_diff_fsigf_from_fsigfano "

Compute from F_sigF_anom to F_sigF (difference structure factor)  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1datatypes_1_1Compute__EsigE__from__FsigF.xml


%feature("docstring") clipper::datatypes::Compute_EsigE_from_FsigF "

Compute from F_sigF to E_sigE.  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1datatypes_1_1Compute__fphi__from__fsigf__phifom.xml


%feature("docstring") clipper::datatypes::Compute_fphi_from_fsigf_phifom "

Compute from F_sigF+Phi_fom to F_phi.  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1datatypes_1_1Compute__FsigF.xml


%feature("docstring") clipper::datatypes::Compute_FsigF "

Compute from F_sigF to F_sigF.  

Use this to get F_sigF from F_sigF of a different precision or from F_sigF_ano.  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1datatypes_1_1Compute__mean__fsigf__from__fsigfano.xml


%feature("docstring") clipper::datatypes::Compute_mean_fsigf_from_fsigfano "

Compute from F_sigF_anom to F_sigF (mean structure factor)  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1datatypes_1_1Compute__neg__fphi.xml


%feature("docstring") clipper::datatypes::Compute_neg_fphi "

Negate F_phi (i.e. advance phase by pi)  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1datatypes_1_1Compute__phifom__from__abcd.xml


%feature("docstring") clipper::datatypes::Compute_phifom_from_abcd "

Compute from ABCD to Phi_fom by phase integration (loses bimodality)  

C++ includes: hkl_compute.h
";

%feature("docstring") clipper::datatypes::Compute_phifom_from_abcd::Compute_phifom_from_abcd "
";

// File: classclipper_1_1datatypes_1_1Compute__scale__u.xml


%feature("docstring") clipper::datatypes::Compute_scale_u "

Deprecated
Apply scale and U to any scalable datatype  

C++ includes: hkl_compute.h
";

%feature("docstring") clipper::datatypes::Compute_scale_u::Compute_scale_u "

constructor: takes scale, U value  

DEPRECATED: This operator applies the scale factor against intensities and the U
value against magnitudes, which is counterintuitive. Compute_scale_u_iso is more
intutive.  

Construct conversion operator to scale list according to the formula I_new =
s.exp( b.|h|^2/2 ) I_old or F_new^2 = s.exp( b.|h|^2/2 ) F_old^2 where |h| =
invresolsq.  

Parameters
----------
* `s` :  
    The intensity scale factor.  
* `u` :  
    The temperature factor (U-value).  
";

// File: classclipper_1_1datatypes_1_1Compute__scale__u__aniso.xml


%feature("docstring") clipper::datatypes::Compute_scale_u_aniso "

Apply scale and U to any scalable datatype.  

C++ includes: hkl_compute.h
";

%feature("docstring") clipper::datatypes::Compute_scale_u_aniso::Compute_scale_u_aniso "

constructor: takes scale, U value  

Construct conversion operator to scale list according to the formula I_new = s^2
.exp( 4 ^2 h^T U h ) I_old or F_new = s.exp( 2 ^2 h^T U h ) F_old  

Parameters
----------
* `s` :  
    The scale factor.  
* `u` :  
    The temperature factor (U-value).  
";

// File: classclipper_1_1datatypes_1_1Compute__scale__u__iso.xml


%feature("docstring") clipper::datatypes::Compute_scale_u_iso "

Apply scale and U to any scalable datatype.  

C++ includes: hkl_compute.h
";

%feature("docstring") clipper::datatypes::Compute_scale_u_iso::Compute_scale_u_iso "

constructor: takes scale, U value  

Construct conversion operator to scale list according to the formula I_new = s^2
.exp( 4 ^2 u.|h|^2 ) I_old or F_new = s.exp( 2 ^2 u.|h|^2 ) F_old where |h|^2 =
invresolsq.  

Parameters
----------
* `s` :  
    The scale factor.  
* `u` :  
    The temperature factor (U-value).  
";

// File: classclipper_1_1datatypes_1_1Compute__sub__fphi.xml


%feature("docstring") clipper::datatypes::Compute_sub_fphi "

Subtract two F_phi datalists.  

C++ includes: hkl_compute.h
";

// File: classclipper_1_1Container.xml


%feature("docstring") clipper::Container "

Definition for a generic container Object.  

Container is a definition for a generic container object with a name, parents,
and children. Any object that wants to be part of the tree simply subclasses
this class. The class also implements search and move objects. The tree is
navigate using unix-like pathnames. A recursive update method can be overridden
to update content after altering the hierarchy.  

The top container in a tree is created by passing Container() as its parent.  

C++ includes: container.h
";

%feature("docstring") clipper::Container::find_path_ptr "

find an object using a directory-like path (NULL on fail)  
";

%feature("docstring") clipper::Container::parent_ptr "

get the parent of this object (NULL on fail)  
";

%feature("docstring") clipper::Container::is_destroyed_with_parent "

is this object to be destroyed when parent is destroyed?  
";

%feature("docstring") clipper::Container::set_destroyed_with_parent "

set this object to be destroyed when parent is destroyed  
";

%feature("docstring") clipper::Container::Container "

constructor: make null object or top object in a tree  
";

%feature("docstring") clipper::Container::Container "

constructor: from any other member and a relative path  
";

%feature("docstring") clipper::Container::~Container "

destructor: virtual  
";

%feature("docstring") clipper::Container::debug "
";

%feature("docstring") clipper::Container::set_name "

set the name of this tree object  
";

%feature("docstring") clipper::Container::update "

update: hierarchical content update function  
";

%feature("docstring") clipper::Container::name "

get the name of this tree object  
";

%feature("docstring") clipper::Container::move "

'move' method moves this object to somewhere else in the hierarchy  
";

%feature("docstring") clipper::Container::has_parent "

test if this object has a parent  
";

%feature("docstring") clipper::Container::path "

get the path of this tree object  
";

%feature("docstring") clipper::Container::parent_of_type_ptr "

search up the tree for a parent of the specified type (NULL on fail)  
";

%feature("docstring") clipper::Container::ultimate_parent "

get the ultimate parent of this object - the top of the tree  
";

%feature("docstring") clipper::Container::ultimate_parent "

get the ultimate parent of this object - the top of the tree  
";

%feature("docstring") clipper::Container::num_children "

return number of children  
";

%feature("docstring") clipper::Container::parent "

get the parent of this object  
";

%feature("docstring") clipper::Container::parent "

get the parent of this object  
";

%feature("docstring") clipper::Container::child "

get the i'th child of this object  
";

%feature("docstring") clipper::Container::child "

get the i'th child of this object  
";

// File: classclipper_1_1Coord__frac.xml


%feature("docstring") clipper::Coord_frac "

fractional (cell) coordinates  

C++ includes: coords.h
";

%feature("docstring") clipper::Coord_frac::transform "

return transformed coordinate  
";

%feature("docstring") clipper::Coord_frac::lattice_copy_unit "

return lattice copy in unit box (0...1,0...1,0...1)  
";

%feature("docstring") clipper::Coord_frac::lattice_copy_near "

return lattice copy near the specified coordinate  
";

%feature("docstring") clipper::Coord_frac::symmetry_copy_near "

return symmetry copy near the specified coordinate  
";

%feature("docstring") clipper::Coord_frac::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::Coord_frac::u "

get u  
";

%feature("docstring") clipper::Coord_frac::w "

get w  
";

%feature("docstring") clipper::Coord_frac::v "

get v  
";

%feature("docstring") clipper::Coord_frac::lengthsq "

return square of length of vector in Angstroms  

Returns
-------
The squared length in Angstroms squared  
";

%feature("docstring") clipper::Coord_frac::coord_map "

fractional-grid coordinate conversion  

Parameters
----------
* `g` :  
    The grid concerned  

Returns
-------
The transformed coordinate.  
";

%feature("docstring") clipper::Coord_frac::lattice_copy_zero "

return lattice copy nearest origin  
";

%feature("docstring") clipper::Coord_frac::coord_orth "

fractional-orthogonal coordinate conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed coordinate.  
";

%feature("docstring") clipper::Coord_frac::coord_grid "

fractional-grid coordinate conversion  

Parameters
----------
* `g` :  
    The grid concerned  

Returns
-------
The transformed coordinate.  
";

%feature("docstring") clipper::Coord_frac::Coord_frac "

null constructor  
";

%feature("docstring") clipper::Coord_frac::Coord_frac "

constructor: copy/convert  
";

%feature("docstring") clipper::Coord_frac::Coord_frac "

constructor: from u,v,w  
";

// File: classclipper_1_1Coord__grid.xml


%feature("docstring") clipper::Coord_grid "

Grid coordinate.  

C++ includes: coords.h
";

%feature("docstring") clipper::Coord_grid::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::Coord_grid::last "

test if done in storage order (see index())  

Test whether this coordinate has been incremented using next() beyond the end of
the specified zero based grid.  

Parameters
----------
* `g` :  
    The grid concerned.  
";

%feature("docstring") clipper::Coord_grid::last "

test if done in storage order (see index())  

Test whether this coordinate has been incremented using next() beyond the end of
the specified non-zero based grid.  

Parameters
----------
* `g` :  
    The grid concerned.  
";

%feature("docstring") clipper::Coord_grid::index "

grid indexing operator  

Return the index in a 1-d array corresponding to this coordinate for a zero
based grid.  

Parameters
----------
* `g` :  
    The grid concerned.  

Returns
-------
The corresponding index.  
";

%feature("docstring") clipper::Coord_grid::coord_frac "

convert to Coord_frac using given Grid_sampling  

Fractional coordinate is not normalised onto range 0..1  

Parameters
----------
* `g` :  
    The grid concerned  

Returns
-------
The fractional coordinate  
";

%feature("docstring") clipper::Coord_grid::coord_map "

convert to Coord_map  

Returns
-------
The non-integer coordinate.  
";

%feature("docstring") clipper::Coord_grid::u "

get u  
";

%feature("docstring") clipper::Coord_grid::u "

set u  
";

%feature("docstring") clipper::Coord_grid::w "

get w  
";

%feature("docstring") clipper::Coord_grid::w "

set w  
";

%feature("docstring") clipper::Coord_grid::v "

get v  
";

%feature("docstring") clipper::Coord_grid::v "

set v  
";

%feature("docstring") clipper::Coord_grid::Coord_grid "

null constructor  
";

%feature("docstring") clipper::Coord_grid::Coord_grid "

constructor: copy/convert  
";

%feature("docstring") clipper::Coord_grid::Coord_grid "

constructor: from u,v,w  
";

%feature("docstring") clipper::Coord_grid::Coord_grid "

constructor: from a grid and an index in that grid  
";

%feature("docstring") clipper::Coord_grid::next "

increment in storage order (see index())  

guaranteed to increment index(g) by 1  

The grid coordinate is incremented efficiently in a manner which is exaclty
equivalent to increasing index() by 1 in a zero based grid.  

Parameters
----------
* `g` :  
    The grid with which this increment is synchronised.  
";

%feature("docstring") clipper::Coord_grid::next "

increment in storage order (see index())  

guaranteed to increment index(g) by 1  

The grid coordinate is incremented efficiently in a manner which is exaclty
equivalent to increasing index() by 1 in a non-zero based grid.  

Parameters
----------
* `g` :  
    The grid with which this increment is synchronised.  
";

%feature("docstring") clipper::Coord_grid::deindex "

grid deindexing operator  

Return the coordinate corresponding to a given index in a zero based grid.  

Parameters
----------
* `g` :  
    The grid concerned.  

Returns
-------
The corresponding coordinate.  
";

%feature("docstring") clipper::Coord_grid::transform "

return transformed coordinate  
";

%feature("docstring") clipper::Coord_grid::unit "

reduce to unit box: (0..nu-1, 0..nv-1, 0..nw-1)  

Parameters
----------
* `g` :  
    The grid concerned  

Returns
-------
The transformed coordinate.  
";

// File: classclipper_1_1Coord__map.xml


%feature("docstring") clipper::Coord_map "

map coordinate: this is like Coord_grid, but non-integer  

C++ includes: coords.h
";

%feature("docstring") clipper::Coord_map::v "

get v  
";

%feature("docstring") clipper::Coord_map::w "

get w  
";

%feature("docstring") clipper::Coord_map::u "

get u  
";

%feature("docstring") clipper::Coord_map::Coord_map "

null constructor  
";

%feature("docstring") clipper::Coord_map::Coord_map "

constructor: copy/convert  
";

%feature("docstring") clipper::Coord_map::Coord_map "

constructor: from Coord_grid  
";

%feature("docstring") clipper::Coord_map::Coord_map "

constructor: from u,v,w  
";

%feature("docstring") clipper::Coord_map::ceil "

return integer Coord_grid above this coordinate  
";

%feature("docstring") clipper::Coord_map::coord_frac "

grid-fractional coordinate conversion  

Parameters
----------
* `g` :  
    The grid concerned  

Returns
-------
The transformed coordinate.  
";

%feature("docstring") clipper::Coord_map::floor "

return integer Coord_grid below this coordinate  
";

%feature("docstring") clipper::Coord_map::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::Coord_map::coord_grid "

return integer Coord_grid nearest this coordinate  
";

// File: classclipper_1_1Coord__orth.xml


%feature("docstring") clipper::Coord_orth "

orthogonal (Angstrom) coordinates  

C++ includes: coords.h
";

%feature("docstring") clipper::Coord_orth::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::Coord_orth::lengthsq "

return square of length of vector in Angstroms  

Returns
-------
The squared length in Angstroms squared  
";

%feature("docstring") clipper::Coord_orth::Coord_orth "

null constructor  
";

%feature("docstring") clipper::Coord_orth::Coord_orth "

constructor: copy/convert  
";

%feature("docstring") clipper::Coord_orth::Coord_orth "

constructor: from x,y,z  
";

%feature("docstring") clipper::Coord_orth::Coord_orth "

constructor: from 3 coords and bond length, angle, torsion  

The coordinate is calculated which extends the sequence of coordinates x1, x2,
x3 with the specified distance to x3, angle to x2,x3, and torsion to x1,x2,x3.  

Parameters
----------
* `x1` :  
    First coordinate.  
* `x2` :  
    Second coordinate.  
* `x3` :  
    Third coordinate.  
* `length` :  
    x3-new bond length in Angstroms.  
* `angle` :  
    x2-x3-new opening angle in Radians.  
* `torsion` :  
    x1-x2-x3-new torsion angle in Radians.  
";

%feature("docstring") clipper::Coord_orth::length "

Return length of vector between two coord orths.  

Returns
-------
The bond length x1-x2 in Angstroms.  
";

%feature("docstring") clipper::Coord_orth::torsion "

Return torsion between four coord orths.  

Returns
-------
The bond torsion x1-x2-x3-x4 in Radians.  
";

%feature("docstring") clipper::Coord_orth::x "

get x  
";

%feature("docstring") clipper::Coord_orth::y "

get y  
";

%feature("docstring") clipper::Coord_orth::z "

get z  
";

%feature("docstring") clipper::Coord_orth::transform "

return transformed coordinate  
";

%feature("docstring") clipper::Coord_orth::angle "

Return angle between three coord orths.  

Returns
-------
The bond angle x1-x2-x3 in Radians.  
";

%feature("docstring") clipper::Coord_orth::coord_frac "

orthogonal-fractional coordinate conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed coordinate.  
";

// File: classclipper_1_1Coord__reci__frac.xml


%feature("docstring") clipper::Coord_reci_frac "

fractional reciprocal coordinate (i.e. non-integer hkl)  

C++ includes: coords.h
";

%feature("docstring") clipper::Coord_reci_frac::vs "

get v*  
";

%feature("docstring") clipper::Coord_reci_frac::Coord_reci_frac "

null constructor  
";

%feature("docstring") clipper::Coord_reci_frac::Coord_reci_frac "

constructor: copy/convert  
";

%feature("docstring") clipper::Coord_reci_frac::Coord_reci_frac "

constructor: from u,v,w  
";

%feature("docstring") clipper::Coord_reci_frac::Coord_reci_frac "

constructor: from HKL  
";

%feature("docstring") clipper::Coord_reci_frac::ws "

get w*  
";

%feature("docstring") clipper::Coord_reci_frac::transform "

return transformed coordinate  
";

%feature("docstring") clipper::Coord_reci_frac::invresolsq "

return inverse resolution squared for this reflection in given cell  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The inverse resolution squared.  
";

%feature("docstring") clipper::Coord_reci_frac::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::Coord_reci_frac::coord_reci_orth "

fractional-orthogonal reciprocal space coordinate conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed coordinate.  
";

%feature("docstring") clipper::Coord_reci_frac::us "

get u*  
";

%feature("docstring") clipper::Coord_reci_frac::hkl "

round to HKL  
";

// File: classclipper_1_1Coord__reci__orth.xml


%feature("docstring") clipper::Coord_reci_orth "

orthogonal reciprocal coordinate (length of which is invresolsq)  

C++ includes: coords.h
";

%feature("docstring") clipper::Coord_reci_orth::xs "

get x*  
";

%feature("docstring") clipper::Coord_reci_orth::zs "

get z*  
";

%feature("docstring") clipper::Coord_reci_orth::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::Coord_reci_orth::Coord_reci_orth "

null constructor  
";

%feature("docstring") clipper::Coord_reci_orth::Coord_reci_orth "

constructor: copy/convert  
";

%feature("docstring") clipper::Coord_reci_orth::Coord_reci_orth "

constructor: from x*,y*,z*  
";

%feature("docstring") clipper::Coord_reci_orth::transform "

return transformed coordinate  
";

%feature("docstring") clipper::Coord_reci_orth::coord_reci_frac "

orthogonal-fractional reciprocal space coordinate conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed coordinate.  
";

%feature("docstring") clipper::Coord_reci_orth::invresolsq "

return inverse resolution squared for this coord  

Returns
-------
The inverse resolution squared.  
";

%feature("docstring") clipper::Coord_reci_orth::ys "

get y*  
";

// File: classclipper_1_1CResolution.xml


%feature("docstring") clipper::CResolution "

Resolution container.  

CResolution: This has a name and a resolution. It overrides the base resolution
for any objects below it.  

C++ includes: container_types.h
";

%feature("docstring") clipper::CResolution::init "

initialiser: from Resolution  

The object is initialised, and children are updated.  

Parameters
----------
* `resolution_` :  
    The value to give to the contained object.  
";

%feature("docstring") clipper::CResolution::CResolution "

constructor: make null object or top object in tree  
";

%feature("docstring") clipper::CResolution::CResolution "

constructor: make child object  
";

// File: classclipper_1_1CSpacegroup.xml


%feature("docstring") clipper::CSpacegroup "

Spacegroup container.  

CSpacegroup: This has a name and a spacegroup. It overrides the base spacegroup
for any objects below it.  

C++ includes: container_types.h
";

%feature("docstring") clipper::CSpacegroup::init "

initialiser: from Spacegroup  

The object is initialised, and children are updated.  

Parameters
----------
* `spacegroup_` :  
    The value to give to the contained object.  
";

%feature("docstring") clipper::CSpacegroup::CSpacegroup "

constructor: make null object or top object in tree  
";

%feature("docstring") clipper::CSpacegroup::CSpacegroup "

constructor: make child object  
";

// File: classclipper_1_1Curv__frac.xml


%feature("docstring") clipper::Curv_frac "

fractional (cell) curvatures, with respect to fractional u,v,w  

C++ includes: derivs.h
";

%feature("docstring") clipper::Curv_frac::Curv_frac "

null constructor  
";

%feature("docstring") clipper::Curv_frac::Curv_frac "

constructor: copy/convert  
";

%feature("docstring") clipper::Curv_frac::curv_map "

fractional-grid derivative conversion  

Parameters
----------
* `g` :  
    The grid concerned  

Returns
-------
The transformed derivative.  
";

%feature("docstring") clipper::Curv_frac::curv_orth "

fractional-orthogonal derivative conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed derivative.  
";

// File: classclipper_1_1Curv__map.xml


%feature("docstring") clipper::Curv_map "

map coordinate curvatures, with respect to grid u,v,w  

C++ includes: derivs.h
";

%feature("docstring") clipper::Curv_map::curv_frac "

grid-fractional derivative conversion  

Parameters
----------
* `g` :  
    The grid concerned  

Returns
-------
The transformed derivative.  
";

%feature("docstring") clipper::Curv_map::Curv_map "

null constructor  
";

%feature("docstring") clipper::Curv_map::Curv_map "

constructor: copy/convert  
";

// File: classclipper_1_1Curv__orth.xml


%feature("docstring") clipper::Curv_orth "

orthogonal (Angstom) curvatures, with respect to orthogonal x,y,z  

C++ includes: derivs.h
";

%feature("docstring") clipper::Curv_orth::curv_frac "

orthogonal-fractional derivative conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed derivative.  
";

%feature("docstring") clipper::Curv_orth::Curv_orth "

null constructor  
";

%feature("docstring") clipper::Curv_orth::Curv_orth "

constructor: copy/convert  
";

// File: classclipper_1_1CXmap.xml


%feature("docstring") clipper::CXmap "

Crystallographic map container.  

CXmap: This is a crystallographic map.  

C++ includes: container_map.h
";

%feature("docstring") clipper::CXmap::CXmap "

null constructor  
";

%feature("docstring") clipper::CXmap::CXmap "

constructor: inherit spacegroup, cell and grid  

The object is constructed at the given location in the hierarchy. An attempt is
made to initialise the object using information from its parents in the
hierarchy.  

Parameters
----------
* `parent` :  
    An object in the hierarchy (usually the parent of the new object).  
* `name` :  
    The path from `parent` to the new object (usually just the name of the new
    object).  
";

%feature("docstring") clipper::CXmap::update "

hierarchical update  

Hierarchical update. If this object is uninitialised, an attempt is made to
initialise the object using information from its parents in the hierarchy. The
childen of the object are then updated.  
";

%feature("docstring") clipper::CXmap::init "

initialiser: supply or inherit spacegroup, cell and grid  

An attempt is made to initialise the object using information from the supplied
parameters, or if they are Null, from its parents in the hierarchy.  

Parameters
----------
* `spacegroup` :  
    The spacegroup for the map.  
* `name` :  
    The cell for the map.  
* `grid` :  
    The grid for the map.  
";

// File: classclipper_1_1datatypes_1_1D__sigD.xml


%feature("docstring") clipper::datatypes::D_sigD "

Deprecated
Anomalous difference data type: D + sigD  

Provided for i/o compatibbility with legacy code only. Do not use.  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::D_sigD::scale "

this type is scalable - apply magnitude scale factor  
";

%feature("docstring") clipper::datatypes::D_sigD::sigd "
";

%feature("docstring") clipper::datatypes::D_sigD::sigd "
";

%feature("docstring") clipper::datatypes::D_sigD::d "
";

%feature("docstring") clipper::datatypes::D_sigD::d "
";

%feature("docstring") clipper::datatypes::D_sigD::missing "
";

%feature("docstring") clipper::datatypes::D_sigD::data_size "
";

%feature("docstring") clipper::datatypes::D_sigD::shift_phase "
";

%feature("docstring") clipper::datatypes::D_sigD::set_null "
";

%feature("docstring") clipper::datatypes::D_sigD::friedel "
";

%feature("docstring") clipper::datatypes::D_sigD::D_sigD "
";

%feature("docstring") clipper::datatypes::D_sigD::D_sigD "
";

%feature("docstring") clipper::datatypes::D_sigD::data_export "
";

%feature("docstring") clipper::datatypes::D_sigD::data_names "
";

%feature("docstring") clipper::datatypes::D_sigD::type "
";

%feature("docstring") clipper::datatypes::D_sigD::data_import "
";

// File: classclipper_1_1Datatype__base.xml


%feature("docstring") clipper::Datatype_base "

Reflection data type objects.  

A class from which data type objects are usually derived  

To define a new type for use in an HKL_data structure, subclass this class and
override the following methods:  

*   constructor - initialises to NAN  
*   type() - returns type name, which is a list of the contained data  
*   friedel() - applies Fridel transformation  
*   shift_phase() - applies phase shift transformation  
*   missing() - checks if data is present  
*   data_size() - number of data elements in this type  
*   data_names() - names of data elements in this type  
*   data_export() - conversion to array (for I/O)  
*   data_import() - conversion from array (for I/O)  
*   scale() - (OPTIONAL) apply magnitude scale factor to data  

note: polymorphism is NOT used here because virtual tables would be to expensive
    for every individual reflection, both in terms of space and cpu cycles.  

C++ includes: hkl_data.h
";

// File: classclipper_1_1datatypes_1_1E__sigE.xml


%feature("docstring") clipper::datatypes::E_sigE "

Reflection data type: E + sigE.  

This is not strictly a type for storing E values, but rather a type for storing
any sturcture factor magnitude-like quantity which has already had a symmetry
enhancement factor (epsilon) removed from it. E's are most commonly stored in
this form, wheras F's and U's are not. You can compute corrected F's from
uncorrected F's using:  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::E_sigE::type "
";

%feature("docstring") clipper::datatypes::E_sigE::shift_phase "
";

%feature("docstring") clipper::datatypes::E_sigE::set_null "
";

%feature("docstring") clipper::datatypes::E_sigE::E_sigE "
";

%feature("docstring") clipper::datatypes::E_sigE::E_sigE "
";

%feature("docstring") clipper::datatypes::E_sigE::E_pl "
";

%feature("docstring") clipper::datatypes::E_sigE::data_export "
";

%feature("docstring") clipper::datatypes::E_sigE::missing "
";

%feature("docstring") clipper::datatypes::E_sigE::E "
";

%feature("docstring") clipper::datatypes::E_sigE::E "
";

%feature("docstring") clipper::datatypes::E_sigE::scale "

this type is scalable - apply magnitude scale factor  
";

%feature("docstring") clipper::datatypes::E_sigE::sigE_mi "
";

%feature("docstring") clipper::datatypes::E_sigE::data_import "
";

%feature("docstring") clipper::datatypes::E_sigE::E_mi "
";

%feature("docstring") clipper::datatypes::E_sigE::data_names "
";

%feature("docstring") clipper::datatypes::E_sigE::sigE_pl "
";

%feature("docstring") clipper::datatypes::E_sigE::cov "
";

%feature("docstring") clipper::datatypes::E_sigE::sigE "
";

%feature("docstring") clipper::datatypes::E_sigE::sigE "
";

%feature("docstring") clipper::datatypes::E_sigE::data_size "
";

%feature("docstring") clipper::datatypes::E_sigE::friedel "
";

// File: classclipper_1_1datatypes_1_1E__sigE__ano.xml


%feature("docstring") clipper::datatypes::E_sigE_ano "

Reflection data type: E(+) E(+) sigE(+) sigE(-) cov+-.  

see datatypes::E_sigE  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::E_sigE_ano::scale "

this type is scalable - apply magnitude scale factor  
";

%feature("docstring") clipper::datatypes::E_sigE_ano::data_export "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::E "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::missing "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::type "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::data_names "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::shift_phase "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::friedel "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::cov "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::cov "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::set_null "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::E_sigE_ano "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::E_mi "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::E_mi "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::E_pl "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::E_pl "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::sigE_pl "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::sigE_pl "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::data_import "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::data_size "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::sigE "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::sigE_mi "
";

%feature("docstring") clipper::datatypes::E_sigE_ano::sigE_mi "
";

// File: classclipper_1_1Euler.xml


%feature("docstring") clipper::Euler "

Euler angle class.  

C++ includes: rotation.h
";

%feature("docstring") clipper::Euler::Euler "

constructor: null  
";

%feature("docstring") clipper::Euler::Euler "

constructor: from specified angles  
";

%feature("docstring") clipper::Euler::Euler "

constructor: from rotation  
";

%feature("docstring") clipper::Euler::rotation "

return rotation  
";

%feature("docstring") clipper::Euler::gamma "

return gamma  
";

%feature("docstring") clipper::Euler::format "

return formatted String representation  
";

%feature("docstring") clipper::Euler::beta "

return beta  
";

%feature("docstring") clipper::Euler::alpha "

return alpha  
";

// File: classclipper_1_1Euler__ccp4.xml


%feature("docstring") clipper::Euler_ccp4 "

Euler_ccp4 angle class.  

C++ includes: rotation.h
";

%feature("docstring") clipper::Euler_ccp4::Euler_ccp4 "

constructor: null  
";

%feature("docstring") clipper::Euler_ccp4::Euler_ccp4 "

constructor: from specified angles  
";

%feature("docstring") clipper::Euler_ccp4::beta "

return beta  
";

%feature("docstring") clipper::Euler_ccp4::gamma "

return gamma  
";

%feature("docstring") clipper::Euler_ccp4::format "

return formatted String representation  
";

%feature("docstring") clipper::Euler_ccp4::alpha "

return alpha  
";

// File: classclipper_1_1datatypes_1_1F__phi.xml


%feature("docstring") clipper::datatypes::F_phi "

Reflection data type: F + phi model or map coeff (e.g. Fcalc, Fbest)  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::F_phi::data_size "
";

%feature("docstring") clipper::datatypes::F_phi::data_names "
";

%feature("docstring") clipper::datatypes::F_phi::friedel "
";

%feature("docstring") clipper::datatypes::F_phi::resolve "

resolve along phase direction  
";

%feature("docstring") clipper::datatypes::F_phi::set_null "
";

%feature("docstring") clipper::datatypes::F_phi::shift_phase "
";

%feature("docstring") clipper::datatypes::F_phi::type "
";

%feature("docstring") clipper::datatypes::F_phi::a "

read real part  
";

%feature("docstring") clipper::datatypes::F_phi::b "

read imag part  
";

%feature("docstring") clipper::datatypes::F_phi::phi "
";

%feature("docstring") clipper::datatypes::F_phi::phi "
";

%feature("docstring") clipper::datatypes::F_phi::missing "
";

%feature("docstring") clipper::datatypes::F_phi::f "
";

%feature("docstring") clipper::datatypes::F_phi::f "
";

%feature("docstring") clipper::datatypes::F_phi::norm "

tidy up so that real part is positive and phase 0...twopi  
";

%feature("docstring") clipper::datatypes::F_phi::scale "

this type is scalable - apply magnitude scale factor  
";

%feature("docstring") clipper::datatypes::F_phi::data_export "
";

%feature("docstring") clipper::datatypes::F_phi::data_import "
";

%feature("docstring") clipper::datatypes::F_phi::F_phi "
";

%feature("docstring") clipper::datatypes::F_phi::F_phi "
";

%feature("docstring") clipper::datatypes::F_phi::F_phi "

convert from complex  
";

// File: classclipper_1_1datatypes_1_1F__sigF.xml


%feature("docstring") clipper::datatypes::F_sigF "

Reflection data type: F + sigF.  

Note that F_sigF also has methods for returning f_pl(), sigf_pl(), f_mi,
sigf_mi(), so you can use this type in any template type where you would use
F_sigF_ano.  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::F_sigF::type "
";

%feature("docstring") clipper::datatypes::F_sigF::f_pl "
";

%feature("docstring") clipper::datatypes::F_sigF::shift_phase "
";

%feature("docstring") clipper::datatypes::F_sigF::F_sigF "
";

%feature("docstring") clipper::datatypes::F_sigF::F_sigF "
";

%feature("docstring") clipper::datatypes::F_sigF::set_null "
";

%feature("docstring") clipper::datatypes::F_sigF::data_size "
";

%feature("docstring") clipper::datatypes::F_sigF::missing "
";

%feature("docstring") clipper::datatypes::F_sigF::sigf "
";

%feature("docstring") clipper::datatypes::F_sigF::sigf "
";

%feature("docstring") clipper::datatypes::F_sigF::sigf_mi "
";

%feature("docstring") clipper::datatypes::F_sigF::sigf_pl "
";

%feature("docstring") clipper::datatypes::F_sigF::data_import "
";

%feature("docstring") clipper::datatypes::F_sigF::data_names "
";

%feature("docstring") clipper::datatypes::F_sigF::cov "
";

%feature("docstring") clipper::datatypes::F_sigF::friedel "
";

%feature("docstring") clipper::datatypes::F_sigF::f_mi "
";

%feature("docstring") clipper::datatypes::F_sigF::f "
";

%feature("docstring") clipper::datatypes::F_sigF::f "
";

%feature("docstring") clipper::datatypes::F_sigF::data_export "
";

%feature("docstring") clipper::datatypes::F_sigF::scale "

this type is scalable - apply magnitude scale factor  
";

// File: classclipper_1_1datatypes_1_1F__sigF__ano.xml


%feature("docstring") clipper::datatypes::F_sigF_ano "

Reflection data type: F(+) F(+) sigF(+) sigF(-) cov+-.  

Note that F_sigF_ano also has methods for returning f(), sigf(), so you can use
this type in any template type where you would use F_sigF.  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::F_sigF_ano::f "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::F_sigF_ano "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::friedel "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::f_pl "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::f_pl "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::type "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::sigf_mi "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::sigf_mi "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::cov "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::cov "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::missing "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::data_size "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::shift_phase "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::set_null "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::sigf_pl "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::sigf_pl "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::data_import "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::scale "

this type is scalable - apply magnitude scale factor  
";

%feature("docstring") clipper::datatypes::F_sigF_ano::data_export "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::sigf "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::f_mi "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::f_mi "
";

%feature("docstring") clipper::datatypes::F_sigF_ano::data_names "
";

// File: classclipper_1_1BasisFn__base_1_1Fderiv.xml


%feature("docstring") clipper::BasisFn_base::Fderiv "

object holding the basis function and its first two derivatives  

C++ includes: resol_fn.h
";

%feature("docstring") clipper::BasisFn_base::Fderiv::Fderiv "

null constructor  
";

%feature("docstring") clipper::BasisFn_base::Fderiv::Fderiv "
";

// File: classclipper_1_1FFTmap.xml


%feature("docstring") clipper::FFTmap "

FFTmap: P1 map with symmetry used for calculating FFTs.  

The FFTmap is represented in P1 in memory. However, it also has a spacegroup,
and the contained data remains consistent with this spacegroup at all times. It
has three states - unassigned, real-space, and reciprocal space. In real space
it contains real map data. In reciprocal space it holds a hemisphere of complex
structure factors, with the Friedels duplicated on the zero section.  

The user should be able to ignore all the issues of spacegroup symmetry, Friedel
opposites, and storage order.  

C++ includes: fftmap.h
";

%feature("docstring") clipper::FFTmap::debug "
";

%feature("docstring") clipper::FFTmap::get_recip_data "

get reciprocal space data  

The data value for the given HKL, or the conjugate of its Friedel opposite if
required, is returned. The symmetry related copies of the data are ignored.  

Parameters
----------
* `rfl` :  
    The HKL of the data to be returned.  
* `fphi` :  
    The value, as a magnitude and phase of type `ffttype`  
";

%feature("docstring") clipper::FFTmap::get_recip_data "

get reciprocal space data (No error checking)  

No error is produced if the space is wrong.  

Parameters
----------
* `rfl` :  
    The HKL of the data to be returned.  

Returns
-------
The value, as magnitude and phase of type `ffttype`  
";

%feature("docstring") clipper::FFTmap::set_recip_data "

set reciprocal space data  

The data value for the given HKL, or the conjugate of its Friedel opposite if
required, is set. All the symmetry related copies of the data, and any Friedel
copies in the zero section, are also set.  

Parameters
----------
* `rfl` :  
    The HKL of the data to be set.  
* `fphi` :  
    The value, as a magnitude and phase of type `ffttype`  
";

%feature("docstring") clipper::FFTmap::spacegroup "

get the spacegroup  
";

%feature("docstring") clipper::FFTmap::reset "

Reset.  

Reset the space and zero all the data, if necessary.  
";

%feature("docstring") clipper::FFTmap::grid_sampling "

get the cell grid  
";

%feature("docstring") clipper::FFTmap::FFTmap "

Null constructor.  

For later initialisation: see init()  
";

%feature("docstring") clipper::FFTmap::FFTmap "

Constructor: takes spacegroup, cell, grid.  

Construct an FFTmap for a given spacegroup, cell, and grid. The map values are
initialised to zero.  

The FFTmap is initially in neither real nor reciprocal spce, however as soon as
one of the 'set' methods is called, it will be defined as in either real or
reciprocal space until the next fft.  

Parameters
----------
* `spacegroup` :  
    The spacegroup.  
* `cell` :  
    The cell, used for scaling.  
* `grid_sam` :  
    The grid sampling of the unit cell.  
* `precalc` :  
    Perform slow precalculation to get faster FFT. (default: no)  
";

%feature("docstring") clipper::FFTmap::cell "

get the cell  
";

%feature("docstring") clipper::FFTmap::init "

initialiser  

Initialise an FFTmap for a given spacegroup, cell, and grid. The map values are
initialised to zero.  

The FFTmap is initially in neither real nor reciprocal spce, however as soon as
one of the 'set' methods is called, it will be defined as in either real or
reciprocal space until the next fft.  

Parameters
----------
* `spacegroup` :  
    The spacegroup.  
* `cell` :  
    The cell, used for scaling.  
* `grid_sam` :  
    The grid sampling of the unit cell.  
* `precalc` :  
    Perform slow precalculation to get faster FFT. This adds a penalty of about
    4s on Linux for the first FFT of any grid and direction. Subsequent FFTs
    will be faster. Set to true for programs which will use many FFTs. default:
    false.  
";

%feature("docstring") clipper::FFTmap::fft_h_to_x "

Transform to real space.  

The data is transformed from recirocal to real space. A scale factor of 1/v
(where v is the cell volume) is applied. If the FFTmap is already in real space,
no action is taken.  
";

%feature("docstring") clipper::FFTmap::fft_map_to_rfl "

calculate reflection-like object from map-like object  

Fill this FFTmap object from a map object, transform it, and fill the given
reflection object from the FFTmap. This will work for any reflection data object
which implements a HKL_reference_index, and every map data object which
implements a Map_reference_index.  

For the results to be sensible, the spacegroup, cell and grids should match.
(The map will be zeroed if necessary).  

Parameters
----------
* `x` :  
    The source map object.  
* `h` :  
    The target reflection data object.  
";

%feature("docstring") clipper::FFTmap::get_real_data "

get real space data  

The data value for the given grid coordinate is returned. Symmetry related
copies are ignored.  

Parameters
----------
* `c` :  
    The coordinate of the data to be returned.  
* `datum` :  
    The value of the data.  
";

%feature("docstring") clipper::FFTmap::get_real_data "

get real space data (No error checking)  

No error is produced if the space is wrong.  

Parameters
----------
* `c` :  
    The grid coordinate of the data to be returned.  

Returns
-------
The value, as type `ffttype`  
";

%feature("docstring") clipper::FFTmap::fft_x_to_h "

Transform to reciprocal space.  

The data is transformed from real to recirocal space. A scale factor of v/n
(where v is the cell volume and n the number of grid points) is applied. If the
FFTmap is already in reciproal space, no action is taken.  
";

%feature("docstring") clipper::FFTmap::set_real_data "

set real space data  

The data value for the given grid coordinate is set. All the symmetry related
copies of the data are also set.  

Parameters
----------
* `c` :  
    The coordinate of the data to be set.  
* `datum` :  
    The value of the data.  
";

%feature("docstring") clipper::FFTmap::fft_rfl_to_map "

calculate map-like object from reflection-like object  

Fill this FFTmap object from a reflection object, transform it, and fill the
given map object from the FFTmap. This will work for any reflection data object
which implements a HKL_reference_index, and every map data object which
implements a Map_reference_index.  

For the results to be sensible, the spacegroup, cell and grids should match.
(The map will be zeroed if necessary).  

Parameters
----------
* `h` :  
    The source reflection data object.  
* `x` :  
    The target map object.  
";

// File: classclipper_1_1FFTmap__base.xml


%feature("docstring") clipper::FFTmap_base "
";

// File: classclipper_1_1FFTmap__p1.xml


%feature("docstring") clipper::FFTmap_p1 "

FFTmap_p1: low level P1 map used for calculating FFTs.  

This is a pure real P1 map, with an extra section in reciprocal space to allow
generation of the full set of resiprocal space magnitudes. Access is by
Coord_grid in both spaces, and indices must be non-negative and in range. The
first and last sections along the half-length direction only have half the
elements stored, the contents of the other half is ignored.  

C++ includes: fftmap.h
";

%feature("docstring") clipper::FFTmap_p1::real_data "

get real space data  
";

%feature("docstring") clipper::FFTmap_p1::real_data "

set real space data  
";

%feature("docstring") clipper::FFTmap_p1::debug "
";

%feature("docstring") clipper::FFTmap_p1::fft_h_to_x "

Transform to real space.  

The data is transformed from recirocal to real space. If the FFTmap_p1 is
already in real space, no action is taken.  

Parameters
----------
* `Scale` :  
    factor to apply (normally 1/cell_volume).  
";

%feature("docstring") clipper::FFTmap_p1::grid_reci "

Return reciprocal space grid (i.e. half real grid + 1 section).  

The reciprocal grid is half-length, plus one section, in the w direction. The
remainder of the grid may be generated by Hermitian symmetry. When accessing
data with reci_data, the coordinate should always be in this grid. Some points
in this grid are redundent, see FFTmap_p1::uniq_reci().  

Returns
-------
The reciprocal grid.  
";

%feature("docstring") clipper::FFTmap_p1::get_hkl "

get reciprocal space data: slow form with hemisphere check  

This form returns the data for an HKL. The HKL is converted into a grid
reference, and the data, or if necessary the conjugate of the opposite, is
returned.  

Parameters
----------
* `hkl` :  
    The HKL of the data.  
";

%feature("docstring") clipper::FFTmap_p1::set_hkl "

set reciprocal space data: slow form with hemisphere check  

This form returns the data for an HKL. The HKL is converted into a grid
reference, and the data, and if necessary the conjugate of the opposite, is set.  

Parameters
----------
* `hkl` :  
    The HKL of the data.  
";

%feature("docstring") clipper::FFTmap_p1::reset "

Reset.  

Reset the space and zero all the data, if necessary.  
";

%feature("docstring") clipper::FFTmap_p1::FFTmap_p1 "

Null constructor.  

For later initialisation: see init()  
";

%feature("docstring") clipper::FFTmap_p1::FFTmap_p1 "

Copy constructor.  
";

%feature("docstring") clipper::FFTmap_p1::FFTmap_p1 "

Constructor: takes grid.  

Construct an FFTmap_p1 for a given spacegroup, cell, and grid. The map values
are initialised to zero.  

The FFTmap_p1 is initially in neither real nor reciprocal spce, however as soon
as one of the 'set' methods is called, it will be defined as in either real or
reciprocal space until the next fft.  

Parameters
----------
* `grid_sam` :  
    The grid sampling of the unit cell.  
* `type` :  
    Can be FFTmap_p1::Measure, FFTmap_p1::Estimate. Measure performs slow
    precalculation (first time only) to get faster FFT.  
";

%feature("docstring") clipper::FFTmap_p1::uniq_reci "

Test whether a coordinate is in the valid part of the recip. grid.  

The w=0 and w=nw/2 sections contain some duplicated points related by a cetre of
symmetry. On of these is considered to be significant, and the other redundent.
This function returns 'true' for the significant point.  

note: For some calculations it may be quicker to set the whole grid than call
    this function for every coordinate.  

Parameters
----------
* `c` :  
    The coordinate to test. Must be in grid_reci().  

Returns
-------
true if the coordinate is for a significant point.  
";

%feature("docstring") clipper::FFTmap_p1::default_type "

set/get default optimisation type  
";

%feature("docstring") clipper::FFTmap_p1::cplx_data "

get reciprocal space data  
";

%feature("docstring") clipper::FFTmap_p1::cplx_data "

set reciprocal space data  
";

%feature("docstring") clipper::FFTmap_p1::init "

initialiser: takes grid  

Initialise an FFTmap_p1 for a given spacegroup, cell, and grid. The map values
are initialised to zero.  

The FFTmap_p1 is initially in neither real nor reciprocal spce, however as soon
as one of the 'set' methods is called, it will be defined as in either real or
reciprocal space until the next fft.  

Parameters
----------
* `grid_sam` :  
    The grid sampling of the unit cell.  
* `type` :  
    Can be FFTmap_p1::Measure, FFTmap_p1::Estimate. Measure performs slow
    precalculation (first time only) to get faster FFT.  
";

%feature("docstring") clipper::FFTmap_p1::grid_real "

Return real space grid.  

Returns
-------
The grid sampling of the real space grid.  
";

%feature("docstring") clipper::FFTmap_p1::fft_x_to_h "

Transform to reciprocal space.  

The data is transformed from real to recirocal space. If the FFTmap_p1 is
already in reciproal space, no action is taken.  

Parameters
----------
* `Scale` :  
    factor to apply (in addition to 1/N_grid factor) (normally cell_volume).  
";

// File: classclipper_1_1FFTmap__sparse__p1__base.xml


%feature("docstring") clipper::FFTmap_sparse_p1_base "

base type for sparse P1 fft maps  

C++ includes: fftmap_sparse.h
";

%feature("docstring") clipper::FFTmap_sparse_p1_base::grid_real "

get real grid sampling  
";

%feature("docstring") clipper::FFTmap_sparse_p1_base::~FFTmap_sparse_p1_base "

Destructor.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_base::init "

initialiser: takes grid  

Initialise an FFTmap_sparse_p1_base for a grid.  

Parameters
----------
* `grid_sam` :  
    The grid sampling of the unit cell.  
* `type` :  
    Can be FFTmap_sparse_base::Measure, ::Estimate. Measure performs slow
    precalculation (first time only) to get faster FFT.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_base::default_type "

set/get default optimisation type  
";

%feature("docstring") clipper::FFTmap_sparse_p1_base::grid_reci "

get reciprocal grid  
";

// File: classclipper_1_1FFTmap__sparse__p1__hx.xml


%feature("docstring") clipper::FFTmap_sparse_p1_hx "

FFTmap_sparse_p1_hx: low level sparse P1 map used for calculating FFTs.  

This version computes sparse Hermititan...real FFTs.  

By specifying what parts of the map are needed in advance, it is possible to
perform highly optimised FFTs, including some of the benefits of symmetry.  

C++ includes: fftmap_sparse.h
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::FFTmap_sparse_p1_hx "

Null constuctor.  

For later initialisation: see init()  
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::FFTmap_sparse_p1_hx "

Constructor: takes grid.  

Construct an FFTmap_sparse_p1_hx for a given grid.  

Parameters
----------
* `grid_sam` :  
    The grid sampling of the unit cell.  
* `type` :  
    Can be FFTmap_sparse_base::Measure, ::Estimate. Measure performs slow
    precalculation (first time only) to get faster FFT.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::cplx_data "

set reciprocal space data (internal use)  
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::real_data "

get real space data ( uvw must be in grid_real() )  

( uvw must be in grid_sampling(), and have been requested )  

Parameters
----------
* `uvw` :  
    The coordinate to get.  

Returns
-------
The real value at that coordinate.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::require_real_data "

express need for real space data  

The given Coord_grid will be required in the final map. ( uvw must be in
grid_sampling() )  

Parameters
----------
* `uvw` :  
    The coordinate to require.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::grid_reci "
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::grid_real "
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::init "
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::set_hkl "

set reciprocal space data by hkl  

Friedel opposites are handled correctly  

Parameters
----------
* `hkl` :  
    The HKL to set.  
* `f` :  
    The complex value to set.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_hx::fft_h_to_x "

Transform to real space.  

The 'require' functions must have been called first to mark the required data in
the target space. (Source space requirements are inferred automatically).  
";

// File: classclipper_1_1FFTmap__sparse__p1__xh.xml


%feature("docstring") clipper::FFTmap_sparse_p1_xh "

FFTmap_sparse_p1_xh: low level sparse P1 map used for calculating FFTs.  

This version computes sparse real...Hermititan FFTs.  

By specifying what parts of the map are needed in advance, it is possible to
perform highly optimised FFTs, including some of the benefits of symmetry.  

C++ includes: fftmap_sparse.h
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::FFTmap_sparse_p1_xh "

Null constuctor.  

For later initialisation: see init()  
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::FFTmap_sparse_p1_xh "

Constructor: takes grid.  

Construct an FFTmap_sparse_p1_xh for a given grid.  

Parameters
----------
* `grid_sam` :  
    The grid sampling of the unit cell.  
* `type` :  
    Can be FFTmap_sparse_base::Measure, ::Estimate. Measure performs slow
    precalculation (first time only) to get faster FFT.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::require_hkl "

express need for reciprocal space data by hkl  

Friedel opposites are handled correctly  

Parameters
----------
* `hkl` :  
    The HKL required.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::cplx_data "

get reciprocal space data (internal use)  

( hkl must be in grid_reci(), and have been requested )  

Parameters
----------
* `uvw` :  
    The coordinate to get.  

Returns
-------
The complex value at that coordinate.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::real_data "

set real space data ( uvw must be in grid_real() )  

( uvw must be in grid_real() )  
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::require_cplx_data "

express need for reciprocal space data (internal use)  

The given Coord_grid will be required in the final reflections. ( uvw must be in
grid_reci() )  

Parameters
----------
* `uvw` :  
    The coordinate to require.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::grid_reci "
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::get_hkl "

get reciprocal space data by hkl  

Friedel opposites are handled correctly  

Parameters
----------
* `hkl` :  
    The required.  
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::grid_real "
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::init "
";

%feature("docstring") clipper::FFTmap_sparse_p1_xh::fft_x_to_h "

Transform to real space.  

The 'require' functions must have been called first to mark the required data in
the target space. (Source space requirements are inferred automatically).  
";

// File: classclipper_1_1datatypes_1_1Flag.xml


%feature("docstring") clipper::datatypes::Flag "

Reflection data type: Free-R flag.  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::Flag::set_null "
";

%feature("docstring") clipper::datatypes::Flag::friedel "
";

%feature("docstring") clipper::datatypes::Flag::data_export "
";

%feature("docstring") clipper::datatypes::Flag::flag "
";

%feature("docstring") clipper::datatypes::Flag::flag "
";

%feature("docstring") clipper::datatypes::Flag::data_import "
";

%feature("docstring") clipper::datatypes::Flag::data_names "
";

%feature("docstring") clipper::datatypes::Flag::type "
";

%feature("docstring") clipper::datatypes::Flag::data_size "
";

%feature("docstring") clipper::datatypes::Flag::missing "
";

%feature("docstring") clipper::datatypes::Flag::Flag "
";

%feature("docstring") clipper::datatypes::Flag::Flag "
";

%feature("docstring") clipper::datatypes::Flag::shift_phase "
";

// File: classclipper_1_1datatypes_1_1Flag__bool.xml


%feature("docstring") clipper::datatypes::Flag_bool "

Reflection data type: boolean (false = missing)  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::Flag_bool::Flag_bool "
";

%feature("docstring") clipper::datatypes::Flag_bool::data_names "
";

%feature("docstring") clipper::datatypes::Flag_bool::missing "
";

%feature("docstring") clipper::datatypes::Flag_bool::set_null "
";

%feature("docstring") clipper::datatypes::Flag_bool::flag "
";

%feature("docstring") clipper::datatypes::Flag_bool::flag "
";

%feature("docstring") clipper::datatypes::Flag_bool::data_export "
";

%feature("docstring") clipper::datatypes::Flag_bool::data_import "
";

%feature("docstring") clipper::datatypes::Flag_bool::type "
";

%feature("docstring") clipper::datatypes::Flag_bool::data_size "
";

%feature("docstring") clipper::datatypes::Flag_bool::friedel "
";

%feature("docstring") clipper::datatypes::Flag_bool::shift_phase "
";

// File: classclipper_1_1Generic__ordinal.xml


%feature("docstring") clipper::Generic_ordinal "

Generic ordinal gernerator.  

This is a generic fast ordinal calculator. It is supplied with a list of values,
from which it prepares a cumulative distribution function. This may the used to
return the approximate fracitonal ordinal (in the range 0...1) for any given
value from the distibution.  

The distibution may be initialised by providing a vector of values from the
distribution, or by adding the values and calling prep_ordinal().  

This distribution may also be inverted. Generation of a value from an ordinal
may be used for generating random values from a given distribution, or for
histogram matching.  

C++ includes: clipper_stats.h
";

%feature("docstring") clipper::Generic_ordinal::Generic_ordinal "

null constructor  
";

%feature("docstring") clipper::Generic_ordinal::Generic_ordinal "

constructor: from range and sampling  
";

%feature("docstring") clipper::Generic_ordinal::prep_ordinal "

generate the ordinal histogram  
";

%feature("docstring") clipper::Generic_ordinal::init "

initialiser: takes the source range and sampling  
";

%feature("docstring") clipper::Generic_ordinal::init "

initialiser: takes the source distibution and a number of bins  
";

%feature("docstring") clipper::Generic_ordinal::init "

DEPRECATED: initialiser: takes a number of bins for histogram.  
";

%feature("docstring") clipper::Generic_ordinal::invert "

invert distribution to get value from ordinal  
";

%feature("docstring") clipper::Generic_ordinal::ordinal "

return reflection ordinal  
";

%feature("docstring") clipper::Generic_ordinal::add_pass_2 "

DEPRECATED: add a value to the distribution (pass 2 of 2)  
";

%feature("docstring") clipper::Generic_ordinal::add_pass_1 "

DEPRECATED: add a value to the distribution (pass 1 of 2)  
";

%feature("docstring") clipper::Generic_ordinal::accumulate "

accumulate values to build the distribution  
";

%feature("docstring") clipper::Generic_ordinal::accumulate "

accumulate values to build the distribution  
";

// File: classclipper_1_1Grad__frac.xml


%feature("docstring") clipper::Grad_frac "

fractional (cell) gradient, with respect to fractional u,v,w  

C++ includes: derivs.h
";

%feature("docstring") clipper::Grad_frac::du "

get d/du  
";

%feature("docstring") clipper::Grad_frac::dv "

get d/dv  
";

%feature("docstring") clipper::Grad_frac::dw "

get d/dw  
";

%feature("docstring") clipper::Grad_frac::Grad_frac "

null constructor  
";

%feature("docstring") clipper::Grad_frac::Grad_frac "

constructor: copy/convert  
";

%feature("docstring") clipper::Grad_frac::Grad_frac "

constructor: from d/du,d/dv,d/dw  
";

%feature("docstring") clipper::Grad_frac::grad_orth "

fractional-orthogonal derivative conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed derivative.  
";

%feature("docstring") clipper::Grad_frac::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::Grad_frac::grad_map "

fractional-grid derivative conversion  

Parameters
----------
* `g` :  
    The grid concerned  

Returns
-------
The transformed derivative.  
";

// File: classclipper_1_1Grad__map.xml


%feature("docstring") clipper::Grad_map "

map coordinate gradient, with respect to grid u,v,w  

C++ includes: derivs.h
";

%feature("docstring") clipper::Grad_map::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::Grad_map::du "

get d/du  
";

%feature("docstring") clipper::Grad_map::dw "

get d/dw  
";

%feature("docstring") clipper::Grad_map::dv "

get d/dv  
";

%feature("docstring") clipper::Grad_map::Grad_map "

null constructor  
";

%feature("docstring") clipper::Grad_map::Grad_map "

constructor: copy/convert  
";

%feature("docstring") clipper::Grad_map::Grad_map "

constructor: from d/du,d/dv,d/dw  
";

%feature("docstring") clipper::Grad_map::grad_frac "

grid-fractional derivative conversion  

Parameters
----------
* `g` :  
    The grid concerned  

Returns
-------
The transformed derivative.  
";

// File: classclipper_1_1Grad__orth.xml


%feature("docstring") clipper::Grad_orth "

orthogonal (Angstom) gradient, with respect to orthogonal x,y,z  

C++ includes: derivs.h
";

%feature("docstring") clipper::Grad_orth::Grad_orth "

null constructor  
";

%feature("docstring") clipper::Grad_orth::Grad_orth "

constructor: copy/convert  
";

%feature("docstring") clipper::Grad_orth::Grad_orth "

constructor: from d/dx,d/dy,d/dz  
";

%feature("docstring") clipper::Grad_orth::dz "

get d/dz  
";

%feature("docstring") clipper::Grad_orth::dy "

get d/dy  
";

%feature("docstring") clipper::Grad_orth::dx "

get d/dx  
";

%feature("docstring") clipper::Grad_orth::grad_frac "

orthogonal-fractional derivative conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed derivative.  
";

%feature("docstring") clipper::Grad_orth::format "

return formatted String representation  

The result is an RT operator. This is a redudent representation, but is handy
for assembling compound operators.  

Returns
-------
The operator  

The formatted text string  
";

// File: classclipper_1_1Grid.xml


%feature("docstring") clipper::Grid "

generic grid  

This holds the dimensions of a 3D array, indexed from 0 along each dimension.  

C++ includes: coords.h
";

%feature("docstring") clipper::Grid::deindex "

grid deindexing operator  
";

%feature("docstring") clipper::Grid::debug "
";

%feature("docstring") clipper::Grid::in_grid "

determine if a point is in the grid  
";

%feature("docstring") clipper::Grid::size "

return size of grid array  
";

%feature("docstring") clipper::Grid::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::Grid::index "

grid indexing operator  
";

%feature("docstring") clipper::Grid::nu "

get nu  
";

%feature("docstring") clipper::Grid::nv "

get nv  
";

%feature("docstring") clipper::Grid::nw "

get nw  
";

%feature("docstring") clipper::Grid::Grid "

null constructor  
";

%feature("docstring") clipper::Grid::Grid "

constructor: from nu,nv,nw  
";

// File: classclipper_1_1Grid__range.xml


%feature("docstring") clipper::Grid_range "

Grid range class: defines array limits for a grid.  

This class is used for describing 3D grids covering an arbitrary part of the 3D
space, i.e. which do not start from (0,0,0).  

C++ includes: coords.h
";

%feature("docstring") clipper::Grid_range::Grid_range "

null constructor  
";

%feature("docstring") clipper::Grid_range::Grid_range "

constructor: takes grid limits  

Make a map grid with an oblong bounded by the coordinates min and max.  

Parameters
----------
* `min` :  
    The lower bound coordinate in u,v,w.  
* `max` :  
    The upper bound coordinate in u,v,w.  
";

%feature("docstring") clipper::Grid_range::Grid_range "

constructor: takes cell grid and fractional limits  

Make a map grid with an oblong bounded by the fractional coordinates min and
max, when the sampling of the cell is g  

Parameters
----------
* `g` :  
    The grid sampling of the whole unit cell.  
* `min` :  
    The lower bound coordinate in u,v,w.  
* `max` :  
    The upper bound coordinate in u,v,w.  
";

%feature("docstring") clipper::Grid_range::Grid_range "

constructor: make grid to hold a sphere from cell, grid, radius  

Make a map grid large enough to fully enclose a sphere about the origin of a
given radius with a given cell and grid sampling.  

Parameters
----------
* `cell` :  
    The cell parameters.  
* `grid` :  
    The grid sampling of the whole cell.  
* `radius` :  
    The radius of the sphere in Angstroms.  
";

%feature("docstring") clipper::Grid_range::in_grid "

determine if a point is in the grid  
";

%feature("docstring") clipper::Grid_range::deindex "

grid deindexing operator  
";

%feature("docstring") clipper::Grid_range::min "

access grid limits  
";

%feature("docstring") clipper::Grid_range::add_border "

border: increase grid to include given border  

Enlarge the grid by adding `b` cells in every direction. Will shrink the grid if
`b` is negative.  

Parameters
----------
* `b` :  
    The number of cells by which to enlarge/shrink.  
";

%feature("docstring") clipper::Grid_range::index "

grid indexing operator  
";

%feature("docstring") clipper::Grid_range::max "

access grid limits  
";

// File: classclipper_1_1Grid__sampling.xml


%feature("docstring") clipper::Grid_sampling "

Grid sampling of a unit cell.  

    This class represents the grid sampling of a unit cell. It is
 otherwise identical to its parent, clipper::Grid_cell, but has an additional
constructor which takes a spacegroup, cell and resolution and produces an
appropriate grid obeying all of the symmetry constraints, and using efficient
factors for the calculation of FFTs.  

note: The following methods are inherited from Grid and Grid_cell but are
    documented here for convenience: nu(), nv(), nw(), size(), index(),
    deindex(), format(), coord_frac(), coord_grid(), to_unit().  

C++ includes: coords.h
";

%feature("docstring") clipper::Grid_sampling::matrix_grid_frac "

return matrix which converts grid to fractional coordinates  

The result is an RT operator. This is a redudent representation, but is handy
for assembling compound operators.  

Returns
-------
The operator  
";

%feature("docstring") clipper::Grid_sampling::nu "
";

%feature("docstring") clipper::Grid_sampling::Grid_sampling "

null constructor  
";

%feature("docstring") clipper::Grid_sampling::Grid_sampling "

constructor: from nu, nv, nw  
";

%feature("docstring") clipper::Grid_sampling::Grid_sampling "

constructor: from Spacegroup, Cell, Resolution, Shannon rate  

A grid is chosen to represent the specified cell at the given resolution,
obeying any restrictions imposed by the spacegroup. A slightly finer grid may be
chosen if doing so is liable to significantly increase the speed of FFTs on that
grid.  

Parameters
----------
* `spacegroup` :  
    The spacegroup which the grid must obey.  
* `cell` :  
    The cell which the grid must contain.  
* `resol` :  
    The resolution to which the grid must sample.  
* `rate` :  
    The linear Shannon rate (oversampling) required. If rate = 1, the grid
    spaceing will be half the resolution (the the minimum required). For a grid
    spaceing of resol/3, use the default rate=1.5.  
";

%feature("docstring") clipper::Grid_sampling::deindex "
";

%feature("docstring") clipper::Grid_sampling::index "
";

%feature("docstring") clipper::Grid_sampling::init "

initialiser: from Spacegroup, Cell, Resolution, Shannon rate  

A grid is chosen to represent the specified cell at the given resolution,
obeying any restrictions imposed by the spacegroup. A slightly finer grid may be
chosen if doing so is liable to significantly increase the speed of FFTs on that
grid.  

Parameters
----------
* `spacegroup` :  
    The spacegroup which the grid must obey.  
* `cell` :  
    The cell which the grid must contain.  
* `resol` :  
    The resolution to which the grid must sample.  
* `rate` :  
    The linear Shannon rate (oversampling) required. If rate = 1, the grid
    spaceing will be half the resolution (the the minimum required). For a grid
    spaceing of resol/3, use the default rate=1.5.  
";

%feature("docstring") clipper::Grid_sampling::size "
";

%feature("docstring") clipper::Grid_sampling::is_null "

test if object has been initialised  

Returns
-------
true if the object has not been initalised.  
";

%feature("docstring") clipper::Grid_sampling::nw "
";

%feature("docstring") clipper::Grid_sampling::nv "
";

%feature("docstring") clipper::Grid_sampling::matrix_frac_grid "

return matrix which converts fractional to grid coordinates  

The result is an RT operator. This is a redudent representation, but is handy
for assembling compound operators.  

Returns
-------
The operator  
";

%feature("docstring") clipper::Grid_sampling::format "
";

// File: classclipper_1_1Histogram.xml


%feature("docstring") clipper::Histogram "

General histogram class.  

This class is used to accumulate and access a histogram of values spread over a
specified range. On storing data or retrieving by interpolation the range is
checked.  

C++ includes: clipper_stats.h
";

%feature("docstring") clipper::Histogram::size "
";

%feature("docstring") clipper::Histogram::x "
";

%feature("docstring") clipper::Histogram::y "

return value at index in histogram (Note: no bound check on i)  
";

%feature("docstring") clipper::Histogram::y "

return value at interpolated position in histogram  
";

%feature("docstring") clipper::Histogram::sum "

return sum of whole histogram  
";

%feature("docstring") clipper::Histogram::x_max "
";

%feature("docstring") clipper::Histogram::x_min "
";

%feature("docstring") clipper::Histogram::Histogram "

null constructor  
";

%feature("docstring") clipper::Histogram::Histogram "

constructor: from range and sampling  
";

%feature("docstring") clipper::Histogram::accumulate "

add value to histogram (if it is in range)  
";

%feature("docstring") clipper::Histogram::accumulate "

add specified value to histogram (if it is in range)  
";

// File: classclipper_1_1HKL.xml


%feature("docstring") clipper::HKL "

reflection 'Miller' index  

C++ includes: coords.h
";

%feature("docstring") clipper::HKL::HKL "

null constructor  
";

%feature("docstring") clipper::HKL::HKL "

constructor: copy/convert  
";

%feature("docstring") clipper::HKL::HKL "

constructor: from H,K,L  
";

%feature("docstring") clipper::HKL::invresolsq "

return inverse resolution squared for this reflection in given cell  

note: Normally you would get a value through clipper::HKL_info, unless you
    specifically want a value for a different cell.  
";

%feature("docstring") clipper::HKL::transform "

return transformed hkl  

Requires integer->ftype->integer transformation.  

Parameters
----------
* `op` :  
    The symmetry operator  

Returns
-------
The transformed coordinate  
";

%feature("docstring") clipper::HKL::transform "

return transformed hkl  

Optimal version.  

Parameters
----------
* `op` :  
    The symmetry operator  

Returns
-------
The transformed coordinate  
";

%feature("docstring") clipper::HKL::k "

get k  
";

%feature("docstring") clipper::HKL::k "

set k  
";

%feature("docstring") clipper::HKL::coord_reci_frac "

return fractional reciprocal coordinate (i.e. non-integer HKL)  

Returns
-------
The non-integer coordinate.  
";

%feature("docstring") clipper::HKL::h "

get h  
";

%feature("docstring") clipper::HKL::h "

set h  
";

%feature("docstring") clipper::HKL::l "

get l  
";

%feature("docstring") clipper::HKL::l "

set l  
";

%feature("docstring") clipper::HKL::format "

return formatted String representation  

Returns
-------
The formatted text string  
";

%feature("docstring") clipper::HKL::coord_reci_orth "

orthogonal-fractional reciprocal space coordinate conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed coordinate.  
";

%feature("docstring") clipper::HKL::sym_phase_shift "

return symmetry phase shift for this HKL under op  

Get the symmetry phase shift incurred when transforming a reflection by this
operator.  

Parameters
----------
* `hkl` :  
    The reflection HKL to transform.  

Returns
-------
The phase shift.  
";

// File: classclipper_1_1HKL__class.xml


%feature("docstring") clipper::HKL_class "

reflection class  

This describes the type of a reflection in a given spacegroup, including
centricity, systematic absence, phase restriction, and multiplicity.  

C++ includes: coords.h
";

%feature("docstring") clipper::HKL_class::centric "

is centric?  
";

%feature("docstring") clipper::HKL_class::sys_abs "

is sys abs?  
";

%feature("docstring") clipper::HKL_class::epsilon "

get epsilon  
";

%feature("docstring") clipper::HKL_class::epsilonc "

get epsilon for acentric, 2x epsilon for centric  
";

%feature("docstring") clipper::HKL_class::HKL_class "

null constructor  
";

%feature("docstring") clipper::HKL_class::HKL_class "

constructor - from spacegroup and HKL  

Determine the class of a reflection for a give spacegroup.  

Parameters
----------
* `spgr` :  
    The spacegroup.  
* `hkl` :  
    The reflection HKL  
";

%feature("docstring") clipper::HKL_class::allowed "

get allowed phase  
";

// File: classclipper_1_1HKL__data.xml


%feature("docstring") clipper::HKL_data "

HKL_data<>  

An actual hkl_data object, containing actual data of type T. This implements the
generic interface, and in addition provides type-specific access functions.  

note: The following methods are inherited from HKL_data_base but are documented
    here for convenience: base_hkl_info(), base_cell(), invresolsq(),
    invresolsq_range(), num_obs(), first(), first_data(), next_data().  

C++ includes: hkl_data.h
";

%feature("docstring") clipper::HKL_data::set_null "

set data entry in the list to its null value  
";

%feature("docstring") clipper::HKL_data::num_obs "
";

%feature("docstring") clipper::HKL_data::invresolsq "
";

%feature("docstring") clipper::HKL_data::set_data "

set data by HKL_reference_coord (returns false if no equivalent hkl)  

If a symmetry mate of the requested HKL exists in the list, then the correct
symmetry transformations are applied and data is set to the supplied values,
otherwise the function returns false.  

Parameters
----------
* `ih` :  
    The reference to the HKL.  
* `data` :  
    Value of the data to set.  

Returns
-------
true if the data was set.  
";

%feature("docstring") clipper::HKL_data::set_data "

set data by hkl (returns false if no equivalent hkl)  

If a symmetry mate of the requested HKL exists in the list, then the correct
symmetry transformations are applied and data is set to the supplied values,
otherwise the function returns false.  

Parameters
----------
* `hkl` :  
    The reflection HKL.  
* `data` :  
    Value of the data to set.  

Returns
-------
true if the data was set.  
";

%feature("docstring") clipper::HKL_data::get_data "

get data by HKL_reference_coord (returns false if no equivalent hkl)  

If a symmetry mate of the requested HKL exists in the list, then the correct
symmetry transformations are applied and the data is returned, otherwise the
value of 'missing' for the datatype is returned.  

Parameters
----------
* `ih` :  
    The reference to the HKL.  
* `data` :  
    Returned with the value of the data.  

Returns
-------
true if the data was returned.  
";

%feature("docstring") clipper::HKL_data::get_data "

get data by hkl (returns false if no equivalent hkl)  

If a symmetry mate of the requested HKL exists in the list, then the correct
symmetry transformations are applied and the supplied datatype is set, otherwise
the function returns false.  

Parameters
----------
* `hkl` :  
    The reflection HKL.  
* `data` :  
    Returned with the value of the data.  

Returns
-------
true if the data was returned.  
";

%feature("docstring") clipper::HKL_data::data_import "

conversion from array (for I/O)  
";

%feature("docstring") clipper::HKL_data::data_names "

return names of data elements in this type  
";

%feature("docstring") clipper::HKL_data::init "

initialiser: from parent hkl_info and cell  

Initialise the object using a given reflection list and cell.  

Parameters
----------
* `hkl_info` :  
    The reflection list object.  
* `cell` :  
    The unit cell for this datalist.  
";

%feature("docstring") clipper::HKL_data::init "

[CLIPPER2] initialiser: from spacegroup, cell, and HKL_sampling  

Initialise the object using a given spacegroup, cell, and sampling.  

Parameters
----------
* `spacegroup` :  
    The spacegroup for this datalist.  
* `cell` :  
    The unit cell for this datalist.  
* `hkl_sampling` :  
    The reflection list description.  
";

%feature("docstring") clipper::HKL_data::init "

[CLIPPER2] initialiser: from another HKL_data object  

    Initialise the object using a given HKL_data object. The
 properties of the object (spacegroup, cell, sampling) are the copied, but the
actual data is not.  

Parameters
----------
* `hkl_data` :  
    The HKL_data object to provide the data.  
";

%feature("docstring") clipper::HKL_data::debug "
";

%feature("docstring") clipper::HKL_data::data_export "

conversion to array (for I/O)  
";

%feature("docstring") clipper::HKL_data::compute "

Basic computation: fill this data list by function call.  
";

%feature("docstring") clipper::HKL_data::compute "

Unary computation: fill this data list by computation from another.  
";

%feature("docstring") clipper::HKL_data::compute "

Binary computation: fill this data list by computation from another.  
";

%feature("docstring") clipper::HKL_data::mask "

For each data element, if the corresponding element in `mask` is missing, then
that element in this list is also set to missing.  

Parameters
----------
* `mask` :  
    The list to provide the mask.  
";

%feature("docstring") clipper::HKL_data::next_data "
";

%feature("docstring") clipper::HKL_data::update "

update: synchornize info with parent HKL_info  

The datalist is resized if necessary to match the parent.  
";

%feature("docstring") clipper::HKL_data::invresolsq_range "
";

%feature("docstring") clipper::HKL_data::base_hkl_info "
";

%feature("docstring") clipper::HKL_data::first "
";

%feature("docstring") clipper::HKL_data::data_size "

return number of data elements in this type  
";

%feature("docstring") clipper::HKL_data::first_data "
";

%feature("docstring") clipper::HKL_data::missing "

check if a data entry in the list is marked as 'missing'  
";

%feature("docstring") clipper::HKL_data::type "

get data type (a list of names corresponding to the im/export values)  
";

%feature("docstring") clipper::HKL_data::base_cell "
";

%feature("docstring") clipper::HKL_data::HKL_data "

null constructor  
";

%feature("docstring") clipper::HKL_data::HKL_data "

constructor: from parent hkl_info  

Construct the object using a given reflection list and cell.  

Parameters
----------
* `hkl_info` :  
    The reflection list object.  
";

%feature("docstring") clipper::HKL_data::HKL_data "

constructor: from parent hkl_info and cell  

Construct the object using a given reflection list and cell.  

Parameters
----------
* `hkl_info` :  
    The reflection list object.  
* `cell` :  
    The unit cell for this datalist.  
";

%feature("docstring") clipper::HKL_data::HKL_data "

[CLIPPER2] constructor: from spacegroup, cell and hkl_sampling  

    Construct the object using a given spacegroup, cell, and sampling.
  

Parameters
----------
* `spacegroup` :  
    The spacegroup for this datalist.  
* `cell` :  
    The unit cell for this datalist.  
* `hkl_sampling` :  
    The reflection list description.  
";

%feature("docstring") clipper::HKL_data::HKL_data "

[CLIPPER2] constructor: from another HKL_data object  

    Construct the object using a given HKL_data object. The
 properties of the object (spacegroup, cell, sampling) are the copied, but the
actual data is not.  

Parameters
----------
* `hkl_data` :  
    The HKL_data object to provide the data.  
";

// File: classclipper_1_1HKL__data__base.xml


%feature("docstring") clipper::HKL_data_base "

HKL_data_base.  

This is the virtual base for the typed hkl_data objects. It exists to guarantee
and interface by which data can be managed without knowledge of the specific
data type.  

C++ includes: hkl_data.h
";

%feature("docstring") clipper::HKL_data_base::first "

return HKL_reference_index pointing to first reflection  

Returns
-------
HKL reference to the first data in this object.  
";

%feature("docstring") clipper::HKL_data_base::base_cell "

get the parent cell  
";

%feature("docstring") clipper::HKL_data_base::data_import "

conversion from array (for I/O)  
";

%feature("docstring") clipper::HKL_data_base::data_export "

conversion to array (for I/O)  
";

%feature("docstring") clipper::HKL_data_base::data_names "

return names of data elements in this type  
";

%feature("docstring") clipper::HKL_data_base::debug "
";

%feature("docstring") clipper::HKL_data_base::type "

get data type (a list of names corresponding to the im/export values)  
";

%feature("docstring") clipper::HKL_data_base::is_null "

test if object has been initialised  

Returns
-------
true if the object has not been initalised.  
";

%feature("docstring") clipper::HKL_data_base::init "

initialiser: from parent hkl_info and cell  

Initialise the object using a given reflection list and cell.  

Parameters
----------
* `hkl_info` :  
    The reflection list object.  
* `cell` :  
    The unit cell for this datalist.  
";

%feature("docstring") clipper::HKL_data_base::init "

initialiser: from another hkl_data  

Initialise the object using a given reflection list and cell.  

Parameters
----------
* `hkl_data` :  
    Object from which to inherit spacegrpoup, cell, sampling.  
";

%feature("docstring") clipper::HKL_data_base::init "

[CLIPPER2] initialiser: from spacegroup, cell, and HKL_sampling  

Initialise the object using a given spacegroup, cell, and sampling.  

Parameters
----------
* `spacegroup` :  
    The spacegroup for this datalist.  
* `cell` :  
    The unit cell for this datalist.  
* `hkl_sampling` :  
    The reflection list description.  
";

%feature("docstring") clipper::HKL_data_base::missing "

check if a data entry in the list is marked as 'missing'  
";

%feature("docstring") clipper::HKL_data_base::hkl_info "

[CLIPPER2] get HKL_info object  
";

%feature("docstring") clipper::HKL_data_base::data_size "

return number of data elements in this type  
";

%feature("docstring") clipper::HKL_data_base::invresolsq_range "

get resolution limits of the list (based on true cell and missing data)  

Returns
-------
The high and low resolution limits of the non-missing data.  
";

%feature("docstring") clipper::HKL_data_base::cell "

[CLIPPER2] get cell  
";

%feature("docstring") clipper::HKL_data_base::invresolsq "

get resolution by reflection index (based on true cell)  

Return the resolution of a particular reflection. If the cell of this list
closely matches (to within 0.5A) the cell of the parent list, this is a simple
lookup, otherwise a metric calculation is required.  
";

%feature("docstring") clipper::HKL_data_base::spacegroup "

[CLIPPER2] get spacegroup  
";

%feature("docstring") clipper::HKL_data_base::resolution "

[CLIPPER2] get resolution  
";

%feature("docstring") clipper::HKL_data_base::hkl_sampling "

[CLIPPER2] get HKL_sampling  
";

%feature("docstring") clipper::HKL_data_base::set_null "

set data entry in the list to its null value  
";

%feature("docstring") clipper::HKL_data_base::next_data "

increment HKL_reference_index to next non-missing data  

Parameters
----------
* `ih` :  
    The HKL reference to increment.  

Returns
-------
HKL reference to the next non-missing data in this object.  
";

%feature("docstring") clipper::HKL_data_base::base_hkl_info "

get the parent HKL_info object  
";

%feature("docstring") clipper::HKL_data_base::first_data "

return HKL_reference_index pointing to first non-missing data  

Returns
-------
HKL reference to the first non-missing data in this object.  
";

%feature("docstring") clipper::HKL_data_base::mask "

mask the data by marking any data missing in 'mask' as missing  
";

%feature("docstring") clipper::HKL_data_base::update "

update: synchornize info with parent HKL_info  
";

%feature("docstring") clipper::HKL_data_base::num_obs "

get number of observations in this list (based on missing data)  

Returns
-------
The number of non-missing data in the object.  
";

// File: classclipper_1_1HKL__data__cacheobj.xml


%feature("docstring") clipper::HKL_data_cacheobj "
";

%feature("docstring") clipper::HKL_data_cacheobj::format "
";

%feature("docstring") clipper::HKL_data_cacheobj::HKL_data_cacheobj "
";

%feature("docstring") clipper::HKL_data_cacheobj::matches "
";

// File: classclipper_1_1HKL__info.xml


%feature("docstring") clipper::HKL_info "

HKL list container and tree root.  

This object contains contains a reflection list, and all the properties on which
such a list depends, i.e. spacegroup, cell, resolution. It also keeps a fast
reflection lookup list and lookup lists for resolutions and reflection classes.  

C++ includes: hkl_info.h
";

%feature("docstring") clipper::HKL_info::debug "
";

%feature("docstring") clipper::HKL_info::add_hkl_list "

add new reflections to the list  

The new HKLs are transformed to the default reciprocal ASU, and added to the
reflection list. Duplicates and reflections outside the resoluution limit are
ignored. Then the fast lookup tables for HKL, invresolsq, and reflection class
are rebuilt.  

Parameters
----------
* `add` :  
    The list of new reflections to add.  
";

%feature("docstring") clipper::HKL_info::invresolsq_range "

get resolution limits of the list  
";

%feature("docstring") clipper::HKL_info::HKL_info::HKL_reference_index "
";

%feature("docstring") clipper::HKL_info::resolution "

get the resolution  
";

%feature("docstring") clipper::HKL_info::hkl_of "

reflection hkl from index  

Parameters
----------
* `index` :  
    The index.  

Returns
-------
The corresponding HKL.  
";

%feature("docstring") clipper::HKL_info::generate_hkl_list "

synthesize hkl list  

Using current cell, spacegroup, resolution.  
";

%feature("docstring") clipper::HKL_info::HKL_info::HKL_reference_coord "
";

%feature("docstring") clipper::HKL_info::cell "

get the cell  
";

%feature("docstring") clipper::HKL_info::find_sym "

find symop no and friedel to bring an HKL into ASU  

Returns the index of the reflection, the sym no. and Friedel flag.  
";

%feature("docstring") clipper::HKL_info::hkl_class "

get reflection class using lookup  
";

%feature("docstring") clipper::HKL_info::first "

return HKL_reference_index pointing to first reflection  
";

%feature("docstring") clipper::HKL_info::index_of "

reflection index from hkl  

This does not check symmetry equivalents (see find_sym).  

Parameters
----------
* `rfl` :  
    The HKL.  

Returns
-------
The index, or -1 if it does not exist.  
";

%feature("docstring") clipper::HKL_info::init "

initialiser: Takes spacegroup, cell, and resolution  

Initialise the HKL_info object. This updates the spacegroup and cell and clears
the reflection list. The resolution is used as a rejection criterion for
reflections - no HKL will be stored beyond the given limit. Initially there are
no reflections in the reflection list: see generate_hkl_list().  

If any of the parameters have null values, the existing values will be
unchanged. The object will only be fully initialised once all parameters are
available.  

Parameters
----------
* `spacegroup` :  
    The spacegroup.  
* `cell` :  
    The unit cell.  
* `resolution` :  
    The resolution limit.  
* `generate` :  
    If true, a reflection list will be generated for an ASU.  
";

%feature("docstring") clipper::HKL_info::init "

initialiser: Takes spacegroup, cell, and HKL_sampling  

Initialise the HKL_info object. This updates the spacegroup and cell and clears
the reflection list. The HKL_sampling determines the reflection list.  

If any of the parameters have null values, the existing values will be
unchanged. The object will only be fully initialised once all parameters are
available.  

Parameters
----------
* `spacegroup` :  
    The spacegroup.  
* `cell` :  
    The unit cell.  
* `hkl_sampling` :  
    The resolution limit.  
* `generate` :  
    If true, a reflection list will be generated for an ASU.  
";

%feature("docstring") clipper::HKL_info::spacegroup "

get the spacegroup  
";

%feature("docstring") clipper::HKL_info::num_reflections "

get number of reflections in the object  
";

%feature("docstring") clipper::HKL_info::HKL_info::HKL_reference_base "
";

%feature("docstring") clipper::HKL_info::hkl_sampling "

[CLIPPER2] get HKL_sampling  
";

%feature("docstring") clipper::HKL_info::HKL_info "

null constructor  
";

%feature("docstring") clipper::HKL_info::HKL_info "

constructor: Takes spacegroup, cell, and resolution  

Construct and initialise HKL_info object. This updates the spacegroup and cell
and clears the reflection list. The resolution is used as a rejection criterion
for reflections - no HKL will be stored beyond the given limit. Initially there
are no reflections in the reflection list: see generate_hkl_list().  

If any of the parameters have null values, the existing values will be
unchanged. The object will only be fully initialised once all parameters are
available.  

Parameters
----------
* `spacegroup` :  
    The spacegroup.  
* `cell` :  
    The unit cell.  
* `resolution` :  
    The resolution limit.  
";

%feature("docstring") clipper::HKL_info::invresolsq "

get reflection resolution using lookup  
";

%feature("docstring") clipper::HKL_info::is_null "

test if object has been initialised  

Returns
-------
true if the object has not been initalised.  
";

// File: classclipper_1_1HKL__lookup.xml


%feature("docstring") clipper::HKL_lookup "

Fast reflection lookup object.  

This version uses a ragged array index to get reflection indices from Miller
indices.  

C++ includes: hkl_lookup.h
";

%feature("docstring") clipper::HKL_lookup::debug "
";

%feature("docstring") clipper::HKL_lookup::index_of "

lookup function  
";

%feature("docstring") clipper::HKL_lookup::init "

initialise: make a reflection index for a list of HKLs  
";

// File: classclipper_1_1HKL__info_1_1HKL__reference__base.xml


%feature("docstring") clipper::HKL_info::HKL_reference_base "

HKL reference base class.  

This is a reference to an HKL. It forms a base class for index-like and
coordinate-like HKL references. If you write a method which will work with
either, then specify this instead of either of the derived classed.  

C++ includes: hkl_info.h
";

%feature("docstring") clipper::HKL_info::HKL_reference_base::index "

return the current index (-1 if invalid)  
";

%feature("docstring") clipper::HKL_info::HKL_reference_base::invresolsq "

return the inv resol sq for the reflection (assumes index valid)  
";

%feature("docstring") clipper::HKL_info::HKL_reference_base::invresolsq "

return the inv resol sq for the reflection (assumes index valid)  
";

%feature("docstring") clipper::HKL_info::HKL_reference_base::base_hkl_info "

return the parent HKL_info  
";

%feature("docstring") clipper::HKL_info::HKL_reference_base::last "

test if index has gone past last reflection  
";

// File: classclipper_1_1HKL__info_1_1HKL__reference__coord.xml


%feature("docstring") clipper::HKL_info::HKL_reference_coord "

HKL reference with coord-like behaviour.  

This is a reference to an HKL. It behaves like an HKL, but also stores the index
of the corresponding reflection in the reflection list, if it exists, and the
symmetry and friedel operators required to get there.  

note: The following methods are inherited from HKL_reference_base but are
    documented here for convenience: base_hkl_info(), index(), invresolsq(),
    last().  

C++ includes: hkl_info.h
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::prev_h "

decrement h  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::prev_k "

decrement k  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::prev_l "

decrement l  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::sym "

get current symop number  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::next "

increment to next reflection  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::index "
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::invresolsq "
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::last "
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::base_hkl_info "
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::HKL_reference_coord "

Null constructor.  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::HKL_reference_coord "

Constructor: takes parent HKL_info and initial HKL.  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::hkl "

return the current HKL  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::next_l "

increment l  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::next_k "

increment k  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::next_h "

increment h  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::friedel "

get current friedel flag  
";

%feature("docstring") clipper::HKL_info::HKL_reference_coord::set_hkl "

assign from HKL  

The new HKL must exist in the reflection list. The calculation is optimised for
the case when the new HKL is near the old one.  
";

// File: classclipper_1_1HKL__info_1_1HKL__reference__index.xml


%feature("docstring") clipper::HKL_info::HKL_reference_index "

HKL reference with index-like behaviour.  

This is a reference to an HKL. It behaves like a simple index into the
reflection list, but can be easily converted into an HKL as and when required.
It also implements methods for iterating through a reflection list.  

note: The following methods are inherited from HKL_reference_base but are
    documented here for convenience: base_hkl_info(), index(), invresolsq(),
    last().  

C++ includes: hkl_info.h
";

%feature("docstring") clipper::HKL_info::HKL_reference_index::invresolsq "
";

%feature("docstring") clipper::HKL_info::HKL_reference_index::hkl "

return the current HKL  
";

%feature("docstring") clipper::HKL_info::HKL_reference_index::next "

increment to next reflection  
";

%feature("docstring") clipper::HKL_info::HKL_reference_index::HKL_reference_index "

Null constructor.  
";

%feature("docstring") clipper::HKL_info::HKL_reference_index::HKL_reference_index "

Constructor: takes parent HKL_info and initial index.  
";

%feature("docstring") clipper::HKL_info::HKL_reference_index::last "
";

%feature("docstring") clipper::HKL_info::HKL_reference_index::base_hkl_info "
";

%feature("docstring") clipper::HKL_info::HKL_reference_index::index "
";

%feature("docstring") clipper::HKL_info::HKL_reference_index::hkl_class "

return the reflection class for the reflection  
";

// File: classclipper_1_1HKL__sampling.xml


%feature("docstring") clipper::HKL_sampling "

HKL sampling of reciprocal space.  

The HKL_sampling class uniquely describes a P0 reflection list bounded by some
resolution limit in reciprocal space. It is described in terms of large
integers, and so immune from rounding errors once the object is constructed.  

C++ includes: coords.h
";

%feature("docstring") clipper::HKL_sampling::hkl_limit "

return limiting values of H, K, L  

Returned HKL contains maximum possible values of H, K, L respectively.  

Returns
-------
Limiting h,k,l.  
";

%feature("docstring") clipper::HKL_sampling::HKL_sampling "

null constructor  

Null constructor  

Initialise to 'null'  
";

%feature("docstring") clipper::HKL_sampling::HKL_sampling "

constructor: takes parameters of normal or inverse cell  

Initialise using cell and resolution.  
";

%feature("docstring") clipper::HKL_sampling::format "

return formatted String representation  
";

%feature("docstring") clipper::HKL_sampling::resolution "

return approximate resolution given cell  

Returned resolution is an estimate based on highest reflection in list.  

Returns
-------
The resolution.  
";

%feature("docstring") clipper::HKL_sampling::is_null "

test if object has been initialised  
";

%feature("docstring") clipper::HKL_sampling::in_resolution "

test if a reflection is within the resolution limit  
";

// File: structclipper_1_1HKL__lookup_1_1hlookup.xml


%feature("docstring") clipper::HKL_lookup::hlookup "

lookup on h  

C++ includes: hkl_lookup.h
";

%feature("docstring") clipper::HKL_lookup::hlookup::hlookup "
";

// File: classclipper_1_1datatypes_1_1I__sigI.xml


%feature("docstring") clipper::datatypes::I_sigI "

Reflection data type: I + sigI.  

Note that I_sigI also has methods for returning I_pl(), sigI_pl(), I_mi,
sigI_mi(), so you can use this type in any template type where you would use
I_sigI_ano.  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::I_sigI::I_mi "
";

%feature("docstring") clipper::datatypes::I_sigI::cov "
";

%feature("docstring") clipper::datatypes::I_sigI::sigI_pl "
";

%feature("docstring") clipper::datatypes::I_sigI::data_export "
";

%feature("docstring") clipper::datatypes::I_sigI::scale "

this type is scalable - apply magnitude scale factor  
";

%feature("docstring") clipper::datatypes::I_sigI::sigI_mi "
";

%feature("docstring") clipper::datatypes::I_sigI::I "
";

%feature("docstring") clipper::datatypes::I_sigI::I "
";

%feature("docstring") clipper::datatypes::I_sigI::shift_phase "
";

%feature("docstring") clipper::datatypes::I_sigI::sigI "
";

%feature("docstring") clipper::datatypes::I_sigI::sigI "
";

%feature("docstring") clipper::datatypes::I_sigI::data_names "
";

%feature("docstring") clipper::datatypes::I_sigI::data_size "
";

%feature("docstring") clipper::datatypes::I_sigI::missing "
";

%feature("docstring") clipper::datatypes::I_sigI::I_pl "
";

%feature("docstring") clipper::datatypes::I_sigI::data_import "
";

%feature("docstring") clipper::datatypes::I_sigI::friedel "
";

%feature("docstring") clipper::datatypes::I_sigI::I_sigI "
";

%feature("docstring") clipper::datatypes::I_sigI::I_sigI "
";

%feature("docstring") clipper::datatypes::I_sigI::type "
";

%feature("docstring") clipper::datatypes::I_sigI::set_null "
";

// File: classclipper_1_1datatypes_1_1I__sigI__ano.xml


%feature("docstring") clipper::datatypes::I_sigI_ano "

Reflection data type: I(+) I(+) sigI(+) sigI(-) cov+-.  

Note that I_sigI_ano also has methods for returning I(), sigI(), so you can use
this type in any template type where you would use I_sigI.  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::I_sigI_ano::I_sigI_ano "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::data_names "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::friedel "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::sigI_mi "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::sigI_mi "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::set_null "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::scale "

this type is scalable - apply magnitude scale factor  
";

%feature("docstring") clipper::datatypes::I_sigI_ano::sigI "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::I_mi "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::I_mi "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::data_size "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::I "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::shift_phase "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::type "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::data_import "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::cov "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::cov "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::sigI_pl "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::sigI_pl "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::data_export "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::missing "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::I_pl "
";

%feature("docstring") clipper::datatypes::I_sigI_ano::I_pl "
";

// File: classclipper_1_1Interp__cubic.xml


%feature("docstring") clipper::Interp_cubic "

Wrapper class for third-order (cubic) interpolation fns.  

These can be used through the built-in methods in Xmap/NXmap, or passed to
methods to allow a choice of interpolation methods, or directly by providing the
map as an argument. For example:  

C++ includes: map_interp.h
";

%feature("docstring") clipper::Interp_cubic::interp "

Interpolate map M using type T at coord.  

The value of the map at the supplied map coordinate is calculated by third order
(cubic) interpolation based on the surrounding 64 points.  

Parameters
----------
* `pos` :  
    The fractional coord at which the density is to be calcuated.  

Returns
-------
The value of the density at that point.  
";

%feature("docstring") clipper::Interp_cubic::can_interp "

Test if we can interpolate in map M at coord.  

The map is queried to see if interpolation is possible at the given coord. For a
crystallographic map, this is always true. For a non-crystallographic map, this
depends if the point and enough neighbours are in the grid.  

Parameters
----------
* `map` :  
    The map on which to perform the calculation.  
* `pos` :  
    The map coord at which the density is to be calcuated.  
";

%feature("docstring") clipper::Interp_cubic::interp_grad "

The value of the map at the supplied map coordinate and its gradient are
calculated by third order (cubic) interpolation based on the surrounding 64
points.  

Parameters
----------
* `pos` :  
    The fractional coord at which the density is to be calcuated.  
* `val` :  
    The value of the density at that point.  
* `grad` :  
    The interpolated value as a gradient vector with respect to the fractional
    coordinates (see Cell::coord_orth).  
";

%feature("docstring") clipper::Interp_cubic::order "

Order of interpolant.  
";

%feature("docstring") clipper::Interp_cubic::interp_curv "

The value of the map at the supplied map coordinate and its gradient are
calculated by third order (cubic) interpolation based on the surrounding 64
points.  

Parameters
----------
* `pos` :  
    The fractional coord at which the density is to be calcuated.  
* `val` :  
    The value of the density at that point.  
* `grad` :  
    The interpolated value as a gradient vector with respect to the fractional
    coordinates (see Cell::coord_orth).  
";

// File: classclipper_1_1Interp__linear.xml


%feature("docstring") clipper::Interp_linear "

Wrapper class for first-order (linear) interpolation fns.  

These can be used through the built-in methods in Xmap/NXmap, or passed to
methods to allow a choice of interpolation methods, or directly by providing the
map as an argument. For example:  

C++ includes: map_interp.h
";

%feature("docstring") clipper::Interp_linear::order "

Order of interpolant.  
";

%feature("docstring") clipper::Interp_linear::can_interp "

Test if we can interpolate in map M at coord.  

The map is queried to see if interpolation is possible at the given coord. For a
crystallographic map, this is always true. For a non-crystallographic map, this
depends if the point and enough neighbours are in the grid.  

Parameters
----------
* `map` :  
    The map on which to perform the calculation.  
* `pos` :  
    The map coord at which the density is to be calcuated.  
";

%feature("docstring") clipper::Interp_linear::interp "

Interpolate map M using type T at coord.  

The value of the map at the supplied map coordinate is calculated by first order
(linear) interpolation based on 8 neighbouring points.  

Parameters
----------
* `map` :  
    The map on which to perform the calculation.  
* `pos` :  
    The map coord at which the density is to be calcuated.  

Returns
-------
The value of the density at that point.  
";

// File: classclipper_1_1Interp__nearest.xml


%feature("docstring") clipper::Interp_nearest "

Wrapper class for zeroth-order (nearest neighbour) interpolation fns.  

These can be used through the built-in methods in Xmap/NXmap, or passed to
methods to allow a choice of interpolation methods, or directly by providing the
map as an argument. For example:  

C++ includes: map_interp.h
";

%feature("docstring") clipper::Interp_nearest::order "

Order of interpolant.  
";

%feature("docstring") clipper::Interp_nearest::can_interp "

Test if we can interpolate in map M at coord.  

The map is queried to see if interpolation is possible at the given coord. For a
crystallographic map, this is always true. For a non-crystallographic map, this
depends if the point and enough neighbours are in the grid.  

Parameters
----------
* `map` :  
    The map on which to perform the calculation.  
* `pos` :  
    The map coord at which the density is to be calcuated.  
";

%feature("docstring") clipper::Interp_nearest::interp "

Interpolate map M using type T at coord.  

The value of the map at the supplied map coordinate is calculated by zeroth
order (nearest neighbour) interpolation based on 1 point.  

Parameters
----------
* `map` :  
    The map on which to perform the calculation.  
* `pos` :  
    The map coord at which the density is to be calcuated.  

Returns
-------
The value of the density at that point.  
";

// File: classclipper_1_1Isymop.xml


%feature("docstring") clipper::Isymop "

Integerised symmetry matrix.  

This is used for optimised calculations in real and reciprocal space  

C++ includes: symop.h
";

%feature("docstring") clipper::Isymop::Isymop "

null constructor  
";

%feature("docstring") clipper::Isymop::Isymop "

constructor: RTop  
";

%feature("docstring") clipper::Isymop::Isymop "

constructor  

Integerised symops are more efficient when handling integer coordinate types,
e.g. HKL, Coord_grid. The rotation parts of the integerised symop are general
and can be used for any recirpocal space data. The translation part is specific
to an individual grid.  

Parameters
----------
* `symop` :  
    The conventional symop.  
* `grid` :  
    The specific grid.  
";

// File: classclipper_1_1HKL__data__cacheobj_1_1Key.xml


%feature("docstring") clipper::HKL_data_cacheobj::Key "
";

%feature("docstring") clipper::HKL_data_cacheobj::Key::cell_descr "
";

%feature("docstring") clipper::HKL_data_cacheobj::Key::Key "
";

%feature("docstring") clipper::HKL_data_cacheobj::Key::hkl_sampling "
";

%feature("docstring") clipper::HKL_data_cacheobj::Key::spgr_descr "
";

// File: classclipper_1_1Xmap__cacheobj_1_1Key.xml


%feature("docstring") clipper::Xmap_cacheobj::Key "
";

%feature("docstring") clipper::Xmap_cacheobj::Key::grid_sampling "
";

%feature("docstring") clipper::Xmap_cacheobj::Key::Key "
";

%feature("docstring") clipper::Xmap_cacheobj::Key::spgr_descr "
";

// File: structclipper_1_1HKL__lookup_1_1klookup.xml


%feature("docstring") clipper::HKL_lookup::klookup "

lookup on k  

C++ includes: hkl_lookup.h
";

%feature("docstring") clipper::HKL_lookup::klookup::klookup "
";

// File: structclipper_1_1data_1_1LGdata.xml


%feature("docstring") clipper::data::LGdata "
";

// File: structclipper_1_1HKL__lookup_1_1llookup.xml


%feature("docstring") clipper::HKL_lookup::llookup "

lookup on l  

C++ includes: hkl_lookup.h
";

%feature("docstring") clipper::HKL_lookup::llookup::llookup "
";

// File: classclipper_1_1LogPhaseProb.xml


%feature("docstring") clipper::LogPhaseProb "

Log phase probability distribution object.  

This object is used to store and manipulate phase log-probability distributions.
Centrics are handled by two values on the phase circle, acentrics by a list of
values. The values can be indexed like and array. The phase() function returns
the phase corresponding to the given array index. Conversion to and from
Hendrickson-Lattman coefficients is provided. The object is templatised on the
sampling of the phase circle.  

C++ includes: hkl_operators.h
";

%feature("docstring") clipper::LogPhaseProb::size "

return num. of phases  
";

%feature("docstring") clipper::LogPhaseProb::sampling "

return phase sampling  
";

%feature("docstring") clipper::LogPhaseProb::LogPhaseProb "

constructor: from HKL class  
";

%feature("docstring") clipper::LogPhaseProb::phase "

return phase associated with index  
";

%feature("docstring") clipper::LogPhaseProb::get_abcd "

get HL coeffs  
";

%feature("docstring") clipper::LogPhaseProb::set_abcd "

set HL coeffs  
";

%feature("docstring") clipper::LogPhaseProb::set_phi_fom "

set phi/fom  
";

%feature("docstring") clipper::LogPhaseProb::get_phi_fom "

get phi/fom  
";

// File: classclipper_1_1Map__index__sort.xml


%feature("docstring") clipper::Map_index_sort "

Generic map sorting class.  

This class is used to sort a vector of integer indices into a map. This includes
sorting the whole map to get highest or lowest density first, or sorting some
subset, e.g. a list of peak indices. Integer indices are used because they are
the most compact way of referencing a unique map location. e.g.  

C++ includes: map_utils.h
";

%feature("docstring") clipper::Map_index_sort::sort_decreasing "

Sort a list into decreasing order.  

The index is sorted in place  map The map to be sorted.  index The list of
indices to sort.  
";

%feature("docstring") clipper::Map_index_sort::sort_increasing "

Sort a list into increasing order.  

The index is sorted in place  map The map to be sorted.  index The list of
indices to sort.  
";

// File: classclipper_1_1Xmap__base_1_1Map__reference__base.xml


%feature("docstring") clipper::Xmap_base::Map_reference_base "

Map reference base class.  

This is a reference to an Map. It forms a base class for index-like and
coordinate-like Map references. If you write a method which will work with
either, then specify this instead of either of the derived classed.  

C++ includes: xmap.h
";

%feature("docstring") clipper::Xmap_base::Map_reference_base::index "

Get the index into the map data array.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_base::base_xmap "

return the parent Xmap  
";

%feature("docstring") clipper::Xmap_base::Map_reference_base::last "

Check for end of map.  
";

// File: classclipper_1_1NXmap__base_1_1Map__reference__base.xml


%feature("docstring") clipper::NXmap_base::Map_reference_base "

Map reference base class.  

This is a reference to an Map. It forms a base class for index-like and
coordinate-like Map references. If you write a method which will work with
either, then specify this instead of either of the derived classed.  

C++ includes: nxmap.h
";

%feature("docstring") clipper::NXmap_base::Map_reference_base::last "

Check for end of map.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_base::index "

Get the index into the map data array.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_base::base_nxmap "

return the parent NXmap  
";

// File: classclipper_1_1NXmap__base_1_1Map__reference__coord.xml


%feature("docstring") clipper::NXmap_base::Map_reference_coord "

Map reference with coordinate-like behaviour.  

This is a reference to a map coordinate. It behaves like a coordinate, but also
stores the index of the corresponding point in the map. It also implements
methods for iterating through the a map. Since the current coordinate is stored,
coord() is fast. However it required 5 words of storage.  

note: The following methods are inherited from Map_reference_base but are
    documented here for convenience: base_nxmap(), index(), last().  

C++ includes: nxmap.h
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::base_nxmap "
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::coord_orth "

Get current value of orthogonal coordinate.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::set_coord "

Set current value of coordinate - optimised for nearby coords.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::next "

Simple increment.  

Use of this function resets the stored coordinate and sym  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::prev_u "

decrement u  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::prev_w "

decrement w  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::prev_v "

decrement v  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::index "
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::next_u "

increment u  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::next_w "

increment w  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::next_v "

increment v  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::coord "

Get current value of coordinate.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::Map_reference_coord "

Null constructor.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::Map_reference_coord "

Constructor: need parent map.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::Map_reference_coord "

Constructor: need parent map and coord.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_coord::last "
";

// File: classclipper_1_1Xmap__base_1_1Map__reference__coord.xml


%feature("docstring") clipper::Xmap_base::Map_reference_coord "

Map reference with coordinate-like behaviour.  

This is a reference to a map coordinate. It behaves like a coordinate, but also
stores the index of the corresponding point in the map, and the symmetry
operator required to get there. It also implements methods for iterating through
the a map. Since the current coordinate and symmetry are stored, coord() is
fast. However, it requires 1 pointer and 5 words of storage.  

note: The following methods are inherited from Map_reference_base but are
    documented here for convenience: base_xmap(), index(), last().  

C++ includes: xmap.h
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::next "

Simple increment.  

Use of this function resets the stored coordinate and sym  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::base_xmap "
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::sym "

Get current symmetry operator.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::next_v "

increment v  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::next_w "

increment w  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::next_u "

increment u  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::prev_v "

decrement v  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::last "
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::Map_reference_coord "

Null constructor.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::Map_reference_coord "

Constructor: takes parent map.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::Map_reference_coord "

Constructor: takes parent map and coord.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::coord_orth "

Get current value of orthogonal coordinate.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::prev_w "

decrement w  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::coord "

Get current value of coordinate.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::set_coord "

Set current value of coordinate - optimised for nearby coords.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::prev_u "

increment u  
";

%feature("docstring") clipper::Xmap_base::Map_reference_coord::index "
";

// File: classclipper_1_1Xmap__base_1_1Map__reference__index.xml


%feature("docstring") clipper::Xmap_base::Map_reference_index "

Map reference with index-like behaviour.  

This is a reference to a map coordinate. It behaves like a simple index into the
map, but can be easily converted into a coordinate as and when required. It also
implements methods for iterating through the unique portion of a map. It is very
compact, but coord() involves some overhead and loses any information concerning
symmetry equivelents.  

note: The following methods are inherited from Map_reference_base but are
    documented here for convenience: base_xmap(), index(), last().  

C++ includes: xmap.h
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::coord_orth "

Get current value of orthogonal coordinate.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::coord "

Get current grid coordinate.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::next "

Simple increment.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::index_offset "

Index of neighbouring point.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::index "
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::Map_reference_index "

Null constructor.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::Map_reference_index "

Constructor: takes parent map.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::Map_reference_index "

Constructor: takes parent map and coord.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::base_xmap "
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::set_coord "

Set current value of coordinate - optimised for nearby coords.  
";

%feature("docstring") clipper::Xmap_base::Map_reference_index::last "
";

// File: classclipper_1_1NXmap__base_1_1Map__reference__index.xml


%feature("docstring") clipper::NXmap_base::Map_reference_index "

Map reference with index-like behaviour.  

This is a reference to a map coordinate. It behaves like a simple index into the
map, but can be easily converted into a coordinate as and when required. It also
implements methods for iterating through the map. It is very compact, but
coord() involves some overhead.  

note: The following methods are inherited from Map_reference_base but are
    documented here for convenience: base_nxmap(), index(), last().  

C++ includes: nxmap.h
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::index_offset "

Index of neighbouring point.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::base_nxmap "
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::coord "

Get current grid coordinate.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::last "
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::index "
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::next "

Simple increment.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::set_coord "

Set current value of coordinate - optimised for nearby coords.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::Map_reference_index "

Null constructor.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::Map_reference_index "

Constructor: need parent map.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::Map_reference_index "

Constructor: need parent map and coord.  
";

%feature("docstring") clipper::NXmap_base::Map_reference_index::coord_orth "

Get current value of orthogonal coordinate.  
";

// File: classclipper_1_1Map__stats.xml


%feature("docstring") clipper::Map_stats "

Generic map statistics class.  

This class is used to calculate and store the mean and standard deviation of a
generic map object of scalar types (e.g. Xmap, NXmap). If the map contains NaN
values, those points are excluded for the calculation. In the case of an Xmap,
the proper multiplicty corrections are applied to give statistics for a whole
unit cell  

C++ includes: map_utils.h
";

%feature("docstring") clipper::Map_stats::min "

Minimum of map.  
";

%feature("docstring") clipper::Map_stats::mean "

Mean of map.  
";

%feature("docstring") clipper::Map_stats::std_dev "

Std deviation of map.  
";

%feature("docstring") clipper::Map_stats::max "

Maximum of map.  
";

%feature("docstring") clipper::Map_stats::Map_stats "

null constructor  
";

%feature("docstring") clipper::Map_stats::Map_stats "

Constructor: from Xmap.  

For float and double maps  map The map for which moments are to be calculated.  
";

%feature("docstring") clipper::Map_stats::range "

Range.  
";

// File: classclipper_1_1Mat33.xml


%feature("docstring") clipper::Mat33 "

3x3-matrix class  

C++ includes: clipper_types.h
";

%feature("docstring") clipper::Mat33::inverse "

inverse  
";

%feature("docstring") clipper::Mat33::Mat33 "

null constructor  
";

%feature("docstring") clipper::Mat33::Mat33 "

constructor  
";

%feature("docstring") clipper::Mat33::Mat33 "

constructor: copy/convert  
";

%feature("docstring") clipper::Mat33::Mat33 "

constructor: copy/convert from symmetric matrix  
";

%feature("docstring") clipper::Mat33::transpose "

transpose  
";

%feature("docstring") clipper::Mat33::null "

return null matrix (only valid for floating point types)  
";

%feature("docstring") clipper::Mat33::det "

determinant  
";

%feature("docstring") clipper::Mat33::equals "

test equality  
";

%feature("docstring") clipper::Mat33::is_null "

test for null matrix (only valid for floating point types)  
";

%feature("docstring") clipper::Mat33::format "

return formatted String representation  
";

%feature("docstring") clipper::Mat33::identity "

return identity matrix  
";

// File: classclipper_1_1Mat33sym.xml


%feature("docstring") clipper::Mat33sym "

Compressed form for 3x3 symmetric matrix class.  

C++ includes: clipper_types.h
";

%feature("docstring") clipper::Mat33sym::inverse "

inverse  
";

%feature("docstring") clipper::Mat33sym::null "

return null matrix (only valid for floating point types)  
";

%feature("docstring") clipper::Mat33sym::sqrt "

square root  
";

%feature("docstring") clipper::Mat33sym::is_null "

test for null matrix (only valid for floating point types)  
";

%feature("docstring") clipper::Mat33sym::mat02 "

element (0,2)  
";

%feature("docstring") clipper::Mat33sym::format "

return formatted String representation  
";

%feature("docstring") clipper::Mat33sym::Mat33sym "

null constructor  
";

%feature("docstring") clipper::Mat33sym::Mat33sym "

constructor: from Mat33 (does not check for symmetry)  
";

%feature("docstring") clipper::Mat33sym::Mat33sym "

constructor: from Mat33sym  
";

%feature("docstring") clipper::Mat33sym::Mat33sym "

constructor: from coefficients  
";

%feature("docstring") clipper::Mat33sym::mat22 "

element (2,2)  
";

%feature("docstring") clipper::Mat33sym::det "

determinant  
";

%feature("docstring") clipper::Mat33sym::mat00 "

element (0,0)  
";

%feature("docstring") clipper::Mat33sym::mat01 "

element (0,1)  
";

%feature("docstring") clipper::Mat33sym::mat11 "

element (1,1)  
";

%feature("docstring") clipper::Mat33sym::identity "

return identity matrix  
";

%feature("docstring") clipper::Mat33sym::quad_form "

return quadratic form with vector  
";

%feature("docstring") clipper::Mat33sym::mat12 "

element (1,2)  
";

// File: classclipper_1_1Matrix.xml


%feature("docstring") clipper::Matrix "

General matrix class: like Array2d but with numerical methods.  

C++ includes: clipper_types.h
";

%feature("docstring") clipper::Matrix::Matrix "

null constructor  
";

%feature("docstring") clipper::Matrix::Matrix "

constructor  
";

%feature("docstring") clipper::Matrix::Matrix "

constructor  
";

%feature("docstring") clipper::Matrix::eigen "

eigenvalue calculation (square symmetric matrices only)  

Find the Eigenvalues and Eigenvectors of the matrix. Uses the Jacobi method.
Only suitable for small systems (dimension<20). The matrix is replaced by the
matrix of eigenvectors (as columns).  

Parameters
----------
* `sort` :  
    Sort the eigenvalues and vectors, smallest first. (default=true)  

Returns
-------
Vector of eigenvalues.  
";

%feature("docstring") clipper::Matrix::solve "

equation solver (square matrices only)  

Solve the system of linear equations Ax=b for x Uses elimination. Only suitable
for small systems.  
";

// File: classclipper_1_1Message.xml


%feature("docstring") clipper::Message "

Message handler class.  

The message handler is a static class which handles messages and errors. It has
3 properties:  

*   the output stream: to which messages will be directed (default stderr)  
*   message level: messages with a level >= this will be output (default 5)  
*   fatal level: messages with a level >= this will be fatal (default 9)  

Levels may be in the range 1-9. They are priorotised as follows:  

*   1-4 messages, never fatal  
*   5-8 warning, may be fatal  
*   9: always fatal. The fatal level must be greater than or equal to the
    message level, and greater than or equal to 5.  

A message is any object which implements the following methods:  The level
method may be static. Messages are usually derived from Message_base.  

C++ includes: clipper_message.h
";

%feature("docstring") clipper::Message::set_stream "

set the output stream  
";

%feature("docstring") clipper::Message::stream "

return the current stream  
";

%feature("docstring") clipper::Message::set_message_level "

set the current message level  
";

%feature("docstring") clipper::Message::Message "

null constuctor  
";

%feature("docstring") clipper::Message::message_level "

return the current message level  
";

%feature("docstring") clipper::Message::fatal_level "

return the current fatal error level  
";

%feature("docstring") clipper::Message::set_fatal_level "

set the current fatal error level  
";

%feature("docstring") clipper::Message::message "

pass a message  
";

// File: classclipper_1_1Message__base.xml


%feature("docstring") clipper::Message_base "

Base type for messages.  

C++ includes: clipper_message.h
";

%feature("docstring") clipper::Message_base::level "
";

%feature("docstring") clipper::Message_base::text "
";

// File: classclipper_1_1Message__ctor.xml


%feature("docstring") clipper::Message_ctor "

Constructor message (level = 2)  

C++ includes: clipper_message.h
";

%feature("docstring") clipper::Message_ctor::Message_ctor "
";

%feature("docstring") clipper::Message_ctor::level "
";

%feature("docstring") clipper::Message_ctor::text "
";

// File: classclipper_1_1Message__dtor.xml


%feature("docstring") clipper::Message_dtor "

Destructor message (level = 2)  

C++ includes: clipper_message.h
";

%feature("docstring") clipper::Message_dtor::Message_dtor "
";

%feature("docstring") clipper::Message_dtor::text "
";

%feature("docstring") clipper::Message_dtor::level "
";

// File: classclipper_1_1Message__fatal.xml


%feature("docstring") clipper::Message_fatal "

Fatal message (level = 9)  

C++ includes: clipper_message.h
";

%feature("docstring") clipper::Message_fatal::Message_fatal "
";

%feature("docstring") clipper::Message_fatal::level "
";

%feature("docstring") clipper::Message_fatal::text "
";

// File: classclipper_1_1Message__generic.xml


%feature("docstring") clipper::Message_generic "

Generic message.  

C++ includes: clipper_message.h
";

%feature("docstring") clipper::Message_generic::Message_generic "
";

%feature("docstring") clipper::Message_generic::text "
";

%feature("docstring") clipper::Message_generic::level "
";

// File: classclipper_1_1Message__info.xml


%feature("docstring") clipper::Message_info "

Info message (level = 1)  

C++ includes: clipper_message.h
";

%feature("docstring") clipper::Message_info::Message_info "
";

%feature("docstring") clipper::Message_info::text "
";

%feature("docstring") clipper::Message_info::level "
";

// File: classclipper_1_1Message__warn.xml


%feature("docstring") clipper::Message_warn "

Warning message (level = 5)  

C++ includes: clipper_message.h
";

%feature("docstring") clipper::Message_warn::text "
";

%feature("docstring") clipper::Message_warn::level "
";

%feature("docstring") clipper::Message_warn::Message_warn "
";

// File: classclipper_1_1Metric__tensor.xml


%feature("docstring") clipper::Metric_tensor "

Metric tensor.  

The metric tensor is used to determine a distance in real or reciprocal space
using fraction coordinates or Miller indices. It is symmetrical, so only the
upper triangle is stored with the off-diagonal elements doubled.  

C++ includes: cell.h
";

%feature("docstring") clipper::Metric_tensor::lengthsq "

apply metric to vector  
";

%feature("docstring") clipper::Metric_tensor::lengthsq "

apply metric to int vector  
";

%feature("docstring") clipper::Metric_tensor::format "

return formatted String representation  
";

%feature("docstring") clipper::Metric_tensor::Metric_tensor "

null constructor  

Null constructor  
";

%feature("docstring") clipper::Metric_tensor::Metric_tensor "

constructor: takes parameters of normal or inverse cell  

Construct and initialise a metric tensor, given a set of real or reciprocal cell
parameters.  

Parameters
----------
* `a` :  
    Length of **a** axis in Angstroms or reciprocal Angstroms.  
* `b` :  
    Length of **b** axis in Angstroms or reciprocal Angstroms.  
* `c` :  
    Length of **c** axis in Angstroms or reciprocal Angstroms.  
* `alph` :  
    Angle between **b** and **c** in radians.  
* `beta` :  
    Angle between **a** and **c** in radians.  
* `gamm` :  
    Angle between **a** and **b** in radians.  
";

// File: classclipper_1_1Mutex.xml


%feature("docstring") clipper::Mutex "

Mutex class: used for locking and unlocking shared resources.  

Create a mutex for any sharted resource, i.e. non-stack object used by a multi-
threaded program. The lock and unlock methods lock that resource. Recursive
locks are not allowed.  

C++ includes: clipper_thread.h
";

%feature("docstring") clipper::Mutex::lock "

lock the mutex  
";

%feature("docstring") clipper::Mutex::unlock "

unlock the mutex  
";

%feature("docstring") clipper::Mutex::Mutex "

constructor: create the mutex  
";

%feature("docstring") clipper::Mutex::~Mutex "

destructor: destroy the mutex  
";

// File: classclipper_1_1NX__operator.xml


%feature("docstring") clipper::NX_operator "

NX_operator: non-crystal map operator.  

This class holds a reference to a non-crystal map frame from somewhere within a
crystallographic map frame. In the general case, an orthogonal rotation-
translation operator is provided which maps the orthogonal frame of the crystal
space onto the orthogonal frame of the NXmap space.  

The object calculates and stores optimised transformations between the
crystallgoraphic frame (described either in fractional or grid coordinates), and
the NXmap grid. Fast paths are generated automatically if the grids are related.  

C++ includes: nxmap_operator.h
";

%feature("docstring") clipper::NX_operator::NX_operator "

null constructor  

The object is not initialised, and will return is_null().  
";

%feature("docstring") clipper::NX_operator::NX_operator "

constructor: from Xmap, NXmap, and operator  

The operator and inverse operator, together with any possible optimisations, are
constructed to relate the give crystallographic and non-crystallographic grid
frames, using the supplied orthogonal operator.  

Parameters
----------
* `xmap` :  
    An Xmap defining the crystal grid frame.  
* `nxmap` :  
    An NXmap defining the non-crystal grid frame.  
* `rtop` :  
    The operator relating the orthogonal frame of the NXmap onto the orthogonal
    frame of the Xmap.  
";

%feature("docstring") clipper::NX_operator::NX_operator "

constructor: from cell, grid sampling, NXmap, and operator  

The operator and inverse operator, together with any possible optimisations, are
constructed to relate the give crystallographic and non-crystallographic grid
frames, using the supplied orthogonal operator.  

Parameters
----------
* `cell` :  
    The cell defining the crystal grid frame.  
* `grid` :  
    The grid defining the crystal grid frame.  
* `nxmap` :  
    An NXmap defining the non-crystal grid frame.  
* `rtop` :  
    The operator relating the orthogonal frame of the NXmap onto the orthogonal
    frame of the Xmap.  
";

%feature("docstring") clipper::NX_operator::coord_frac "

convert nxmap map coord to xtal frac coord  
";

%feature("docstring") clipper::NX_operator::init "

initialiser:: from Xmap, NXmap, and operator  

The operator and inverse operator, together with any possible optimisations, are
constructed to relate the give crystallographic and non-crystallographic grid
frames, using the supplied orthogonal operator.  

Parameters
----------
* `xmap` :  
    An Xmap defining the crystal grid frame.  
* `nxmap` :  
    An NXmap defining the non-crystal grid frame.  
* `rtop` :  
    The operator relating the orthogonal frame of the NXmap onto the orthogonal
    frame of the Xmap.  
";

%feature("docstring") clipper::NX_operator::init "

initialiser:: from cell, grid sampling, NXmap, and operator  

The operator and inverse operator, together with any possible optimisations, are
constructed to relate the give crystallographic and non-crystallographic grid
frames, using the supplied orthogonal operator.  

Parameters
----------
* `cell` :  
    The cell defining the crystal grid frame.  
* `grid` :  
    The grid defining the crystal grid frame.  
* `nxmap` :  
    An NXmap defining the non-crystal grid frame.  
* `rtop` :  
    The operator relating the orthogonal frame of the NXmap onto the orthogonal
    frame of the Xmap.  
";

%feature("docstring") clipper::NX_operator::debug "
";

%feature("docstring") clipper::NX_operator::nxmap_data "

get value of nxmap at xmap grid coord using fastest appropriate method  

The density of the non-crystal map at the position corresponding to a
crystallographic map grid coordinate is returned. If the grids match exactly
either by pure translation or by rotation+translation, then fast paths are used
to return the requested density directly. Otherwise the supplied interpolation
template is used. No checking is performed for coordinates outside the NXmap.  

Parameters
----------
* `nxmap` :  
    The non-crystal map (NXmap) to be queried.  
* `c` :  
    The grid coordinate in the crystallographic coordinate frame.  

Returns
-------
The value of the NXmap at the requested position.  
";

%feature("docstring") clipper::NX_operator::is_null "

test if object has been initialised  
";

%feature("docstring") clipper::NX_operator::coord_map "

convert xtal frac coord to nxmap map coord  
";

%feature("docstring") clipper::NX_operator::xmap_data "

get value of xmap at nxmap grid coord using fastest appropriate method  

The density of the crystal map at the position corresponding to a non-
crystallographic map grid coordinate is returned. If the grids match exactly
either by pure translation or by rotation+translation, then fast paths are used
to return the requested density directly. Otherwise the supplied interpolation
template is used.  

Parameters
----------
* `xmap` :  
    The crystal map (Xmap) to be queried.  
* `c` :  
    The grid coordinate in the crystallographic coordinate frame.  

Returns
-------
The value of the Xmap at the requested position.  
";

// File: classclipper_1_1NXmap.xml


%feature("docstring") clipper::NXmap "

NXmap<T>: actual non-crystallographic map class.  

The non-crystallographic map class stores a map of arbitrary data type. Unlike
an Xmap it is finite in extent and has no symmetry. An RT operator provides
mapping onto an arbitrary orthogonal coordinate frame. Iterators provide
efficient access to data.  

This is derived from NXmap_base, and adds the templatised data itself and the
methods which deal with it.  

note: The following methods are inherited from NXmap_base but are documented
    here for convenience: grid(), coord_orth(), coord_grid(), first(),
    first_coord().  

C++ includes: nxmap.h
";

%feature("docstring") clipper::NXmap::first_coord "
";

%feature("docstring") clipper::NXmap::get_data "

get a density value for an arbitrary position  
";

%feature("docstring") clipper::NXmap::set_data "

set a density value for an arbitrary position  
";

%feature("docstring") clipper::NXmap::interp_curv "

get map value and curv for map coord using supplied interpolator  

The value of the map at the desired non-grid map coordinate and its gradient and
curvature are calculated using the supplied interpolator template.  

Parameters
----------
* `pos` :  
    The map coord at which the density is to be calcuated.  
* `val` :  
    The value of the density at that point.  
* `grad` :  
    The interpolated gradient vector with respect to the map coordinates (see
    Cell::coord_orth).  
* `curv` :  
    The interpolated curvature matrix with respect to the map coordinates (see
    Cell::coord_orth).  
";

%feature("docstring") clipper::NXmap::first "
";

%feature("docstring") clipper::NXmap::grid "
";

%feature("docstring") clipper::NXmap::coord_map "
";

%feature("docstring") clipper::NXmap::init "

initialiser: takes grid and orthogonal->grid coordinate operator  

Initialise an NXmap to some rhomboid chosen from within a crystal coordinate
space, specified by the grid and a transformation from orthogonal to grid
coordinates.  

Parameters
----------
* `grid` :  
    The grid dimensions of the desired map.  
* `rt` :  
    The rotation/transln op from orthogonal to grid coordinates.  
";

%feature("docstring") clipper::NXmap::init "

initialiser: takes grid, cell, and fraction limits  

Initialise an NXmap to some rhomboid chosen from within a crystal grid
coordinate space, specified by a cell, sampling and box within that grid. This
is useful for creating an NXmap which exactly matches some subregion of a
crystallographic map.  

Parameters
----------
* `cell` :  
    Unit cell defining the crystal space.  
* `grid` :  
    The grid sampling of the given unit cell.  
* `grid_extent` :  
    The map extent within that cell.  
";

%feature("docstring") clipper::NXmap::NXmap "

Null constructor, for later initialisation.  
";

%feature("docstring") clipper::NXmap::NXmap "

Constructor: takes grid and orthogonal->grid coordinate operator.  

Initialise an NXmap to some rhomboid chosen from within a crystal coordinate
space, specified by the grid and a transformation from orthogonal to grid
coordinates.  

Parameters
----------
* `grid` :  
    The grid dimensions of the desired map.  
* `rt` :  
    The rotation/transln op from orthogonal to grid coordinates.  
";

%feature("docstring") clipper::NXmap::NXmap "

Constructor: takes grid, cell, and extent.  

Initialise an NXmap to some rhomboid chosen from within a crystal grid
coordinate space, specified by a cell, sampling and box within that grid. This
is useful for creating an NXmap which exactly matches some subregion of a
crystallographic map.  

Parameters
----------
* `cell` :  
    Unit cell defining the crystal space.  
* `grid` :  
    The grid sampling of the given unit cell.  
* `grid_extent` :  
    The map extent within that cell.  
";

%feature("docstring") clipper::NXmap::coord_orth "
";

%feature("docstring") clipper::NXmap::interp_grad "

get map value and grad for map coord using supplied interpolator  

The value of the map at the desired non-grid map coordinate and its gradient are
calculated using the supplied interpolator template.  

Parameters
----------
* `pos` :  
    The map coord at which the density is to be calcuated.  
* `val` :  
    The value of the density at that point.  
* `grad` :  
    The interpolated gradient vector with respect to the map coordinates (see
    Cell::coord_orth).  
* `curv` :  
    The interpolated curvature matrix with respect to the map coordinates (see
    Cell::coord_orth).  
";

%feature("docstring") clipper::NXmap::interp "

get map value for map coord using supplied interpolator  

The value of the map at the desired non-grid map coordinate are calculated using
the supplied interpolator template.  

Parameters
----------
* `pos` :  
    The map coord at which the density is to be calcuated.  

Returns
-------
The value of the density at that point. map coordinates (see Cell::coord_orth).  
";

// File: classclipper_1_1NXmap__base.xml


%feature("docstring") clipper::NXmap_base "

NXmap_base: base for non-crystallographic map class.  

The non-crystallographic map class stores a map of arbitrary data type. Unlike
an Xmap it is finite in extent and has no symmetry. An RT operator provides
mapping onto an arbitrary orthogonal coordinate frame. Iterators provide
efficient access to data.  

This base contains everything except the data, which is templated in the derived
type clipper::NXmap<T>.  

C++ includes: nxmap.h
";

%feature("docstring") clipper::NXmap_base::multiplicity "

get multiplicity of a map grid point (always 1 for NXmap)  
";

%feature("docstring") clipper::NXmap_base::coord_orth "

convert map coordinate to orthogonal  

Parameters
----------
* `cm` :  
    The grid coordinate to be converted.  

Returns
-------
The equivalent orthogonal coordinate.  
";

%feature("docstring") clipper::NXmap_base::NXmap_base::Map_reference_base "
";

%feature("docstring") clipper::NXmap_base::is_null "

test if object has been initialised  

Returns
-------
true if the object has not been initalised.  
";

%feature("docstring") clipper::NXmap_base::grid "

return the grid dimensions for this map  
";

%feature("docstring") clipper::NXmap_base::coord_map "

convert orthogonal coordinate to map  

Parameters
----------
* `co` :  
    The orthogonal coordinate to be converted.  

Returns
-------
The equivalent grid coordinate.  
";

%feature("docstring") clipper::NXmap_base::NXmap_base::Map_reference_coord "
";

%feature("docstring") clipper::NXmap_base::NXmap_base::Map_reference_index "
";

%feature("docstring") clipper::NXmap_base::first_coord "

return a coord Map_reference_index for this map  
";

%feature("docstring") clipper::NXmap_base::in_map "

is the given coord available in the map?  
";

%feature("docstring") clipper::NXmap_base::in_map "

is the given coord available in the map using the given interpolant?  

Note that the higher the order of the interpolant, the more of the boundary of
the map becomes inaccessible.  

Parameters
----------
* `cm` :  
    The coord_map to test.  

Returns
-------
true if interpolation can be performed at that coordinate.  
";

%feature("docstring") clipper::NXmap_base::first "

return a basic Map_reference_index for this map  
";

// File: classclipper_1_1NXmap__operator.xml


%feature("docstring") clipper::NXmap_operator "

NXmap_operator: non-crystal map operator referencing a particular NXmap.  

This class holds a reference to a non-crystal map object from somewhere within a
crystallographic map frame. In the general case, an orthogonal rotation-
translation operator is provided which maps the orthogonal frame of the crystal
space onto the orthogonal frame of the NXmap space.  

The object calculates and stores optimised transformations between the
crystallgoraphic frame (described either in fractional or grid coordinates), and
the NXmap grid. Fast paths are generated automatically if the grids are related.  

note: This object differes from NX_operator in that it keeps a reference to an
    individual NXmap, which may be used to access that object directly.  

C++ includes: nxmap_operator.h
";

%feature("docstring") clipper::NXmap_operator::NXmap_operator "

null constructor  
";

%feature("docstring") clipper::NXmap_operator::NXmap_operator "

constructor: from Xmap, NXmap, and operator  
";

%feature("docstring") clipper::NXmap_operator::NXmap_operator "

constructor: from cell, grid sampling, NXmap, and operator  
";

%feature("docstring") clipper::NXmap_operator::init "

initialiser:: from Xmap, NXmap, and operator  
";

%feature("docstring") clipper::NXmap_operator::init "

initialiser:: from cell, grid sampling, NXmap, and operator  
";

%feature("docstring") clipper::NXmap_operator::nxmap_data "

access NXmap directly from xmap grid coord using fastest method  
";

%feature("docstring") clipper::NXmap_operator::nxmap "

get the target NXmap of this operator  
";

// File: classclipper_1_1ObjectCache.xml


%feature("docstring") clipper::ObjectCache "

Object Cache manager.  

The object cache is a tool for storing information which may appear several
times. Examples include tables of information for spacegroups or
crystallographic maps. When a new object is created, a check is first done to
see if such an object already exists in the cache, in which case that copy is
used. Otherwise a new copy is added to the cache.  

A cached object must implement:  

*   a constructor from a type T  
*   a method 'matches(T)' which tests if it constructed from that object  
*   a 'format()' method, returning a string description of the contents The type
    T should be a compact unique description of the object.  

Referring to the cache returns an ObjectCache<T>::Reference to a cache object.
This object performs reference counting, which is used for garbage collection.  

To retrieve the actual cached data, use the ObjectCache<T>::Reference::data()
method. The data is held at a fixed memory location, therefore pointers to the
data may be safely kept, as long as they are discarded as soon as the reference
is discarded (at which point garbage collection may occur).  

Ideally this would be a class with static members only, but some compilers have
trouble with static members of template classes.  

Garbage collection modes include:  

*   ObjectCache<T>::NORMAL : Remove an old object only when a new object is
    required and an old object is no longer in use. (default)  
*   ObjectCache<T>::MINMEM : Remove an old object as soon as it is no longer in
    use.  
*   ObjectCache<T>::MAXMEM : Never remove old objects. The more memory hungry
    modes may improve performance for some problems where a new object may be
    created which was already used and destroyed before.  

C++ includes: clipper_memory.h
";

%feature("docstring") clipper::ObjectCache::ObjectCache "

constructor  
";

%feature("docstring") clipper::ObjectCache::set_mode "

set garbage collection mode  
";

%feature("docstring") clipper::ObjectCache::~ObjectCache "

destructor, can message on contents  
";

%feature("docstring") clipper::ObjectCache::cache "

cache or return data by key  
";

%feature("docstring") clipper::ObjectCache::debug "
";

%feature("docstring") clipper::ObjectCache::purge "

purge unreferenced objects from cache  
";

%feature("docstring") clipper::ObjectCache::destroy "

VERY DANGEROUS, DO NOT USE.  
";

// File: classclipper_1_1datatypes_1_1Phi__fom.xml


%feature("docstring") clipper::datatypes::Phi_fom "

Reflection data type: best phi + fom.  

C++ includes: hkl_datatypes.h
";

%feature("docstring") clipper::datatypes::Phi_fom::shift_phase "
";

%feature("docstring") clipper::datatypes::Phi_fom::missing "
";

%feature("docstring") clipper::datatypes::Phi_fom::data_size "
";

%feature("docstring") clipper::datatypes::Phi_fom::Phi_fom "
";

%feature("docstring") clipper::datatypes::Phi_fom::Phi_fom "
";

%feature("docstring") clipper::datatypes::Phi_fom::phi "
";

%feature("docstring") clipper::datatypes::Phi_fom::phi "
";

%feature("docstring") clipper::datatypes::Phi_fom::fom "
";

%feature("docstring") clipper::datatypes::Phi_fom::fom "
";

%feature("docstring") clipper::datatypes::Phi_fom::data_import "
";

%feature("docstring") clipper::datatypes::Phi_fom::friedel "
";

%feature("docstring") clipper::datatypes::Phi_fom::set_null "
";

%feature("docstring") clipper::datatypes::Phi_fom::type "
";

%feature("docstring") clipper::datatypes::Phi_fom::data_names "
";

%feature("docstring") clipper::datatypes::Phi_fom::data_export "
";

// File: classclipper_1_1Polar__ccp4.xml


%feature("docstring") clipper::Polar_ccp4 "

Polar_ccp4 angle class.  

C++ includes: rotation.h
";

%feature("docstring") clipper::Polar_ccp4::Polar_ccp4 "

null constructor  
";

%feature("docstring") clipper::Polar_ccp4::Polar_ccp4 "

constructor: from specified angles  
";

%feature("docstring") clipper::Polar_ccp4::kappa "

return kappa  
";

%feature("docstring") clipper::Polar_ccp4::format "

return formatted String representation  
";

%feature("docstring") clipper::Polar_ccp4::psi "

return omega  
";

%feature("docstring") clipper::Polar_ccp4::omega "

return omega  
";

%feature("docstring") clipper::Polar_ccp4::phi "

return phi  
";

// File: classclipper_1_1Prob__phi__2d.xml


%feature("docstring") clipper::Prob_phi_2d "

2-d angular probability distibution class  

Base for Ramachandran class (and other similar classes, such as a pseudo-
ramachandran plot or the JPD of two phases ).  

C++ includes: ramachandran.h
";

%feature("docstring") clipper::Prob_phi_2d::probability "

get probability for a particular pair of angles  

linear interpolation off of grid  
";

%feature("docstring") clipper::Prob_phi_2d::format "

formatted string representation (as C++ code)  
";

%feature("docstring") clipper::Prob_phi_2d::normalise "

normalise to integrate to 1/(2pi)^2  

normalise mean value to 1/(2pi)^2  
";

%feature("docstring") clipper::Prob_phi_2d::accumulate "

accumulate new table of samples to probability  
";

%feature("docstring") clipper::Prob_phi_2d::accumulate "

accumulate new sample to probability  

linear interpolation onto grid  
";

%feature("docstring") clipper::Prob_phi_2d::init "

initialise: with sampling  
";

%feature("docstring") clipper::Prob_phi_2d::data "

2d read access  
";

%feature("docstring") clipper::Prob_phi_2d::data "

2d write access  
";

// File: classclipper_1_1Property.xml


%feature("docstring") clipper::Property "

Template for a property holding an arbitrary type.  

C++ includes: clipper_memory.h
";

%feature("docstring") clipper::Property::clone "

factory copy method  
";

%feature("docstring") clipper::Property::Property "

constructor: takes contents  
";

%feature("docstring") clipper::Property::~Property "
";

%feature("docstring") clipper::Property::value "

return value of contents  
";

// File: classclipper_1_1Property__base.xml


%feature("docstring") clipper::Property_base "

Base class for properties of arbitrary types.  

C++ includes: clipper_memory.h
";

%feature("docstring") clipper::Property_base::~Property_base "
";

%feature("docstring") clipper::Property_base::clone "

factory copy method  
";

// File: classclipper_1_1PropertyManager.xml


%feature("docstring") clipper::PropertyManager "

Class for holding a list of labelled properties of arbitrary types.  

To add a property list to an object, derive it from this class, or include a
member and mirror the methods. To add a property, simply call
insert_property(label,property). Properties must be objects derived from
clipper::Propert_base. Usually, you can just use the template form,
clipper::Property<T>.  

To read a property which you know exists and is of a particular type, use:  If
you are unsure if a property is present, use the exists_property(label) method.
If you are unsure of a property's type, dynamic cast a pointer and test for
null. e.g.  

C++ includes: clipper_memory.h
";

%feature("docstring") clipper::PropertyManager::exists_property "

test for property  

Parameters
----------
* `label` :  
    The label of the property to be tested.  

Returns
-------
true on success.  
";

%feature("docstring") clipper::PropertyManager::get_property "

get a labelled property from the list  

Parameters
----------
* `label` :  
    The label of the property to be returned.  

Returns
-------
the property object.  
";

%feature("docstring") clipper::PropertyManager::set_property "

add a labelled property to the list  
";

%feature("docstring") clipper::PropertyManager::copy "

copy manager  

This function is used by the copy constructor and assignement operator and is
also useful for derived classes.  
";

%feature("docstring") clipper::PropertyManager::~PropertyManager "

destructor  

Deletes all stored properties.  
";

%feature("docstring") clipper::PropertyManager::delete_property "

delete property  
";

%feature("docstring") clipper::PropertyManager::PropertyManager "

null constructor  
";

%feature("docstring") clipper::PropertyManager::PropertyManager "

copy constructor  

Makes copies of all property objects.  
";

// File: classclipper_1_1Ramachandran.xml


%feature("docstring") clipper::Ramachandran "

Ramachandran plot class.  

This class provides a reference Ramachandran plot for Gly, Pro, other, and
combinations of those types of residues. The source data comes from the best
residues from the 'top500' best-determined structures list of D. C. and J. S.
Richardson, http://kinemage.biochem.duke.edu/index.html  

The Ramachandran plot is normalised in inverse radians squared, so the mean
value of a probability is 1/(2 pi)2.  

C++ includes: ramachandran.h
";

%feature("docstring") clipper::Ramachandran::favored "

test if a pair of angles are in the favored region  
";

%feature("docstring") clipper::Ramachandran::Ramachandran "

null constructor  
";

%feature("docstring") clipper::Ramachandran::Ramachandran "

constructor: from standard plot  

Construct a Ramachandran plot of a given type.  

Parameters
----------
* `type` :  
    The residue type of the plot. Options include: Ramachandran::Gly,
    Ramachandran::Pro, Ramachandran::NonGlyPro, Ramachandran::NonGly,
    Ramachandran::All  
";

%feature("docstring") clipper::Ramachandran::set_thresholds "

change threshholds to different values  

Set thresholds for favorable and allowed regions of the Ramachandran plot. The
US spelling is used because it is the same length as 'allowed'. I should get out
more. Sorry.  

Parameters
----------
* `prob_favored` :  
    The probability threshold for the favored region.  
* `prob_allowed` :  
    The probability threshold for the allowed region.  
";

%feature("docstring") clipper::Ramachandran::probability "

get probability for a particular pair of angles  
";

%feature("docstring") clipper::Ramachandran::allowed "

test if a pair of angles are in the allowed region  
";

%feature("docstring") clipper::Ramachandran::init "

initialise: from standard plot  

Construct a Ramachandran plot of a given type.  

Parameters
----------
* `type` :  
    The residue type of the plot. Options include: Ramachandran::Gly,
    Ramachandran::Pro, Ramachandran::NonGlyPro, Ramachandran::NonGly,
    Ramachandran::All  
";

// File: classclipper_1_1Range.xml


%feature("docstring") clipper::Range "

Range - upper and lower bounds of some type.  

C++ includes: clipper_stats.h
";

%feature("docstring") clipper::Range::min "

minimum value  
";

%feature("docstring") clipper::Range::Range "

null constructor  
";

%feature("docstring") clipper::Range::Range "

constructor  
";

%feature("docstring") clipper::Range::max "

maximum value  
";

%feature("docstring") clipper::Range::truncate "

truncate data to be within range  
";

%feature("docstring") clipper::Range::range "

range = max - min  
";

%feature("docstring") clipper::Range::include "

update limits to include a new datum  
";

%feature("docstring") clipper::Range::contains "

test if data is within limits ( min <= datum <= max )  
";

// File: classclipper_1_1Range__sampling.xml


%feature("docstring") clipper::Range_sampling "

Range sampling: discrete sampling of a real range.  

C++ includes: clipper_stats.h
";

%feature("docstring") clipper::Range_sampling::Range_sampling "

null constructor  
";

%feature("docstring") clipper::Range_sampling::Range_sampling "

constructor: from number of samplings  
";

%feature("docstring") clipper::Range_sampling::Range_sampling "

constructor: from range and number of samplings  
";

%feature("docstring") clipper::Range_sampling::x "

return x-value (0..n) from fractional posn in counting range  
";

%feature("docstring") clipper::Range_sampling::x "

return x-value corresponding to centre of i'th range  
";

%feature("docstring") clipper::Range_sampling::size "

return number of samplings in range  
";

%feature("docstring") clipper::Range_sampling::x_min "

return x-value corresponding to bottom of i'th range  
";

%feature("docstring") clipper::Range_sampling::x_max "

return x-value corresponding to top of i'th range  
";

%feature("docstring") clipper::Range_sampling::index "

return nearest index to particular x-value  
";

%feature("docstring") clipper::Range_sampling::index_bounded "

return nearest index to particular x-value (bounded 0...n-1)  
";

%feature("docstring") clipper::Range_sampling::indexf "

return fractional posn in counting range from x-value (0..n)  
";

// File: classclipper_1_1TargetFn__base_1_1Rderiv.xml


%feature("docstring") clipper::TargetFn_base::Rderiv "

object holding the residual function and first two derivatives  

C++ includes: resol_fn.h
";

// File: classclipper_1_1ObjectCache_1_1Reference.xml


%feature("docstring") clipper::ObjectCache::Reference "

ObjectCache reference class.  

C++ includes: clipper_memory.h
";

%feature("docstring") clipper::ObjectCache::Reference::Reference "
";

%feature("docstring") clipper::ObjectCache::Reference::Reference "
";

%feature("docstring") clipper::ObjectCache::Reference::is_null "
";

%feature("docstring") clipper::ObjectCache::Reference::~Reference "
";

%feature("docstring") clipper::ObjectCache::Reference::data "
";

// File: classclipper_1_1Resolution.xml


%feature("docstring") clipper::Resolution "

Resolution in angstroms.  

This object represents a resolution limit which will be used for all aspects of
a calculation. This is a base for a donor type.  

C++ includes: coords.h
";

%feature("docstring") clipper::Resolution::Resolution "

null constructor  
";

%feature("docstring") clipper::Resolution::Resolution "

constructor: from ftype  

Parameters
----------
* `resol_` :  
    The resolution limit in Angstroms.  
";

%feature("docstring") clipper::Resolution::invresolsq_limit "

get invresolsq limit  

Returns
-------
The resolution limit in inverse squared Angstroms.  
";

%feature("docstring") clipper::Resolution::init "

initialiser: from ftype  

Parameters
----------
* `resol_` :  
    The resolution limit in Angstroms.  
";

%feature("docstring") clipper::Resolution::limit "

get resolution limit  

Returns
-------
The resolution limit in Angstroms.  
";

%feature("docstring") clipper::Resolution::is_null "

test if value has been initialised  

Returns
-------
true if the object has not been initalised.  
";

// File: classclipper_1_1Resolution__ordinal.xml


%feature("docstring") clipper::Resolution_ordinal "

Resolution ordinal gernerator.  

This class is a helper class for functions which need to divide reflections up
by resolution whilst guaranteeing a certain distribution of number of
reflections per range. It takes a list of reflections, one at a time, and
calculates a function to get the approximate ordinal number of a reflection in a
list sorted by resolution.  

C++ includes: resol_basisfn.h
";

%feature("docstring") clipper::Resolution_ordinal::init "

initialiser: takes an HKL_info and uses all reflections.  
";

%feature("docstring") clipper::Resolution_ordinal::init "

initialiser: takes an HKL_data & uses non-missing reflections.  
";

%feature("docstring") clipper::Resolution_ordinal::init "

initialiser: takes an HKL_data + Cell & uses non-missing reflections.  
";

// File: classclipper_1_1ResolutionFn.xml


%feature("docstring") clipper::ResolutionFn "

2nd order resolution function evaluator  

This is an automatic evaluator for arbitrary functions of HKL, most commonly
used for evaluating a function of resolution (such a mean F^2 or sigmaa),
although more general tasks including local scaling of reflections and
anisotropic functions can also be handled. This form is for target functions
which approach zero quadratically, e.g. least-squares targets.  

This version implements a naive Newton-Raphson minimiser, which only uses the
gradient and curvature of the target function, ignoring its value. It is ideal
for quadratic targets with linear basis functions.  

To evaluate a resolution function, this class must be provided with two objects:  

*   The basis function (and gradients), which describes the value of the
    function for any reflection given a set of paramters.  
*   The target function (and derivatives), which is used to determine the values
    of the basis function parameters.  

For example, the following code may be used to calculate a smooth scaling
function to scale one set of data to another using an anisotropic Gaussian
scaling function:  

The most useful isotropic resolution function is the BasisFn_spline, since it is
linear and provides a good fit to most data.  

C++ includes: resol_fn.h
";

%feature("docstring") clipper::ResolutionFn::f "

return the value of the basis function with the current paramters  
";

%feature("docstring") clipper::ResolutionFn::params "

return the values of the parameters  

Returns
-------
The refined basis function parameters  
";

%feature("docstring") clipper::ResolutionFn::ResolutionFn "

constructor: need reflections, basis fn and target fn.  

The constructor performs the full minimisation calculation.  

Parameters
----------
* `hkl_info` :  
    HKL_info object which provides the reflection list.  
* `basisfn` :  
    The basis function used to describe the desired property.  
* `targetfn` :  
    The target function to be minimised.  
* `params` :  
    Initial values for the function parameters.  
* `damp_` :  
    If > 0.0, shifts are fdamped during early cycles to help convergence with
    difficult bases/target conbinations.  
";

%feature("docstring") clipper::ResolutionFn::debug "

print the target, gradient, and curvatures with respect to the params  
";

// File: classclipper_1_1ResolutionFn__nonlinear.xml


%feature("docstring") clipper::ResolutionFn_nonlinear "

2nd order resolution function evaluator  

This is an automatic evaluator for arbitrary functions of HKL, most commonly
used for evaluating a function of resolution (such a mean F^2 or sigmaa),
although more general tasks including local scaling of reflections and
anisotropic functions can also be handled. This form is for target functions
which approach zero quadratically, e.g. least-squares targets.  

note: This version implements a minimiser which uses both Newton-Raphson and
    gradient steps depending on the situation. It can be used for non-quadratic
    targets or non-linear basis functions.  

To evaluate a resolution function, this class must be provided with two objects:  

*   The basis function (and gradients), which describes the value of the
    function for any reflection given a set of paramters.  
*   The target function (and derivatives), which is used to determine the values
    of the basis function parameters.  

C++ includes: resol_fn.h
";

%feature("docstring") clipper::ResolutionFn_nonlinear::ResolutionFn_nonlinear "

constructor: need reflections, basis fn and target fn.  

The constructor performs the full minimisation calculation.  

Parameters
----------
* `hkl_info` :  
    HKL_info object which provides the reflection list.  
* `basisfn` :  
    The basis function used to describe the desired property.  
* `targetfn` :  
    The target function to be minimised.  
* `damp_` :  
    If > 0.0, shifts are fdamped during early cycles to help convergence with
    difficult bases/target conbinations  
";

// File: classclipper_1_1Rotation.xml


%feature("docstring") clipper::Rotation "

Rotation class.  

This class represents a rotation. The internal representation is as a unit
quaternion, which is easily combined, inverted, or converted to or from other
commonly used forms.  

C++ includes: rotation.h
";

%feature("docstring") clipper::Rotation::z "

return z component  
";

%feature("docstring") clipper::Rotation::y "

return y component  
";

%feature("docstring") clipper::Rotation::x "

return x component  
";

%feature("docstring") clipper::Rotation::zero "

return zero rotation  
";

%feature("docstring") clipper::Rotation::is_null "

test for null (uninitialised) rotation  
";

%feature("docstring") clipper::Rotation::null "

return null rotation  
";

%feature("docstring") clipper::Rotation::euler "

< return Euler angles  
";

%feature("docstring") clipper::Rotation::w "

return w component  
";

%feature("docstring") clipper::Rotation::matrix "

return 3x3 matrix  

The resulting rotation matrix would commonly be used to construct a
clipper::RTop_orth.  

Returns
-------
The rotation matrix.  
";

%feature("docstring") clipper::Rotation::Rotation "

null constructor  
";

%feature("docstring") clipper::Rotation::Rotation "

constructor: from generic Euler  
";

%feature("docstring") clipper::Rotation::Rotation "

constructor: from Euler_ccp4  
";

%feature("docstring") clipper::Rotation::Rotation "

constructor: from Polar_ccp4  
";

%feature("docstring") clipper::Rotation::Rotation "

constructor: from Matrix  
";

%feature("docstring") clipper::Rotation::Rotation "

constructor: from components  
";

%feature("docstring") clipper::Rotation::format "

return formatted String representation  
";

%feature("docstring") clipper::Rotation::polar_ccp4 "

return Polar_ccp4 angles  

If omega ~= 0, then phi is set to zero.  

Returns
-------
The Polar_ccp4 angles.  
";

%feature("docstring") clipper::Rotation::norm "

normalise this quaternion  

The normalisation is performed in-place. If a rotation becomes significantly
denormalised, the conversion methods will fail. Therefore it may be safer to
call this before a conversion.  
";

%feature("docstring") clipper::Rotation::inverse "

return inverse rotation  
";

%feature("docstring") clipper::Rotation::euler_ccp4 "

return Euler_ccp4 angles  

If beta ~= 0, then alpha is set to zero.  

Returns
-------
The Euler_ccp4 angles.  
";

%feature("docstring") clipper::Rotation::abs_angle "

return absolute rotation angle  

Positive magnitude of the angle of rotation.  

Returns
-------
The angle in radians.  
";

// File: classclipper_1_1RTop.xml


%feature("docstring") clipper::RTop "

Rotation-translation operator.  

C++ includes: clipper_types.h
";

%feature("docstring") clipper::RTop::RTop "

null constructor  
";

%feature("docstring") clipper::RTop::RTop "

constructor: from rotation  
";

%feature("docstring") clipper::RTop::RTop "

constructor: from rotation and translation  
";

%feature("docstring") clipper::RTop::rot "

get rotation  
";

%feature("docstring") clipper::RTop::rot "

set rotation  
";

%feature("docstring") clipper::RTop::equals "

test equality with some tolerance  
";

%feature("docstring") clipper::RTop::format "

return formatted String representation  
";

%feature("docstring") clipper::RTop::trn "

get translation  
";

%feature("docstring") clipper::RTop::trn "

set translation  
";

%feature("docstring") clipper::RTop::is_null "

test for null operator  
";

%feature("docstring") clipper::RTop::null "

return identity operator  
";

%feature("docstring") clipper::RTop::identity "

return identity operator  
";

%feature("docstring") clipper::RTop::inverse "

inverse  
";

// File: classclipper_1_1RTop__frac.xml


%feature("docstring") clipper::RTop_frac "

Fractional operator class.  

This class is used for any RT-operator which operates on fractional coordinates.
For a full list of methods, see clipper::RTop  

C++ includes: symop.h
";

%feature("docstring") clipper::RTop_frac::rtop_orth "

fractional-orthogonal conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed coordinate.  
";

%feature("docstring") clipper::RTop_frac::RTop_frac "

null constructor  
";

%feature("docstring") clipper::RTop_frac::RTop_frac "

constructor: copy/convert  
";

%feature("docstring") clipper::RTop_frac::RTop_frac "

constructor: from rotation  
";

%feature("docstring") clipper::RTop_frac::RTop_frac "

constructor: from string description  

Construct an RT operator from a string description, e.g. 1/2x,z-y+2/3,x '*' is
optional for multiplication, commas are compulsory.  
";

%feature("docstring") clipper::RTop_frac::RTop_frac "

constructor: from rotation and translation  
";

%feature("docstring") clipper::RTop_frac::null "

return null (uninitialised) operator  

Returns
-------
The null (uninitialised) operator.  
";

%feature("docstring") clipper::RTop_frac::inverse "

inverse operator  

Returns
-------
The inverse of the operator.  
";

%feature("docstring") clipper::RTop_frac::identity "

return identity operator  

Returns
-------
The identity operator.  
";

// File: classclipper_1_1RTop__orth.xml


%feature("docstring") clipper::RTop_orth "

Orthogonal operator class.  

This class is used for any RT-operator which operates on orthogonal coordinates.
For a full list of methods, see clipper::RTop  

C++ includes: coords.h
";

%feature("docstring") clipper::RTop_orth::RTop_orth "

null constructor  
";

%feature("docstring") clipper::RTop_orth::RTop_orth "

constructor: copy/convert  
";

%feature("docstring") clipper::RTop_orth::RTop_orth "

constructor: from rotation  
";

%feature("docstring") clipper::RTop_orth::RTop_orth "

constructor: from rotation and translation  
";

%feature("docstring") clipper::RTop_orth::RTop_orth "

constructor: from two vectors of Coord_orth  

Construct the operator which give the least-squares fit of one set of
coordinates onto another. The coodinates are stored as STL vectors of
Coord_orth. The lists must be the same size, and each atom in the source list
must correspond to the same atom in the target list. The algorithm employed is
that of Kearsley, S.K. (1989) 'On the orthogonal transformation used for
structural comparisons'. Acta Cryst. A45, 208-210.  

Parameters
----------
* `src` :  
    The source list (i.e. the atoms to be transformed).  
* `tgt` :  
    The target list (i.e. the fixed atoms).  
";

%feature("docstring") clipper::RTop_orth::RTop_orth "

constructor: from two vectors of Coord_orth  

Construct the operator which give the least-squares fit of one set of
coordinates onto another. The coodinates are stored as STL vectors of
Coord_orth. The lists must be the same size, and each atom in the source list
must correspond to the same atom in the target list. The algorithm employed is
that of Kearsley, S.K. (1989) 'On the orthogonal transformation used for
structural comparisons'. Acta Cryst. A45, 208-210.  

Parameters
----------
* `src` :  
    The source list (i.e. the atoms to be transformed).  
* `tgt` :  
    The target list (i.e. the fixed atoms).  
* `wgt` :  
    The weight to apply to each atom.  
";

%feature("docstring") clipper::RTop_orth::RTop_orth "

constructor: from two atom-list type objects  

Construct the operator which relates one atom-list like object onto another. The
lists must be the same size, and have the following properties:  

*   a size() method.  
*   a [int] operator, with int ranging from 0 to size()-1.  
*   the object returned by the [] operator must have a coord_orth() method.
    Suitable objects include a vector of Atom, or an Atom_list.  
";

%feature("docstring") clipper::RTop_orth::screw_translation "

return screw translation  

Returns
-------
screw translation, 000 if rotation is zero  
";

%feature("docstring") clipper::RTop_orth::axis_coordinate_near "

return point on axis near the specified coordinate  

Parameters
----------
* `centre` :  
    An arbitrary point.  

Returns
-------
point on axis near the specified coordinate, 000 if rotation is zero  
";

%feature("docstring") clipper::RTop_orth::identity "

return identity operator  

Returns
-------
The identity operator.  
";

%feature("docstring") clipper::RTop_orth::inverse "

inverse operator  

Returns
-------
The inverse of the operator.  
";

%feature("docstring") clipper::RTop_orth::null "

return null (uninitialised) operator  

Returns
-------
The null (uninitialised) operator.  
";

%feature("docstring") clipper::RTop_orth::rtop_frac "

orthogonal-fractional conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed coordinate.  
";

// File: structclipper_1_1data_1_1SFdata.xml


%feature("docstring") clipper::data::SFdata "
";

// File: structclipper_1_1data_1_1SGdata.xml


%feature("docstring") clipper::data::SGdata "
";

// File: classclipper_1_1Spacegroup.xml


%feature("docstring") clipper::Spacegroup "

Spacegroup object.  

The spacegroup object is a full description of a spacegroup, including all the
most regularly used information in an efficient form. It may be initialised from
a clipper::Spgr_descr. This object.  

For more details of spacegroup symbols, see Sydney R. Hall & Ralf W. Grosse-
Kunstleve 'Concise Space-Group Symbols',
http://www.kristall.ethz.ch/LFK/software/sginfo/hall_symbols.html  

C++ includes: spacegroup.h
";

%feature("docstring") clipper::Spacegroup::descr "

get spacegroup description  
";

%feature("docstring") clipper::Spacegroup::symop "

get n'th symop  
";

%feature("docstring") clipper::Spacegroup::centering_symop "

get n'th centering symop (0...3 max)  
";

%feature("docstring") clipper::Spacegroup::null "

Return null spacegroup.  
";

%feature("docstring") clipper::Spacegroup::recip_asu "

test if hkl is in default reciprocal ASU  

The reciprocal ASU is chosen from one of 47 optimised functions.  

Parameters
----------
* `hkl` :  
    The HKL to test.  

Returns
-------
true if the HKL is in the ASU.  
";

%feature("docstring") clipper::Spacegroup::primitive_symop "

get n'th primitive symop (identical to symop(sym_no))  
";

%feature("docstring") clipper::Spacegroup::num_primops "

get number of primitive symops (identical to num_primitive_symops())  
";

%feature("docstring") clipper::Spacegroup::hkl_class "

get 'class' of reflection: multiplicity, allowed phase, absence  

The reflection class describes the type of a reflection in a given spacegroup,
including centricity, systematic absence, phase restriction, and multiplicity.  

This is a shortcut to constructing an HKL_class from the spacegroup and HKL.  

Parameters
----------
* `hkl` :  
    The reflection HKL  
";

%feature("docstring") clipper::Spacegroup::order_of_symmetry_about_axis "

get the order of rotational symmetry about a given axis  

The number of rotational operators parallel to the specified axis is returned.  

Parameters
----------
* `axis` :  
    The axis, A, B or C.  

Returns
-------
The order of the axis.  
";

%feature("docstring") clipper::Spacegroup::inverse_op "

get symop number corresponding to the inverse of a symop  
";

%feature("docstring") clipper::Spacegroup::inversion_symop "

get n'th inversion symop (0...1 max)  
";

%feature("docstring") clipper::Spacegroup::num_primitive_symops "

get number of primitive symops (inc identity and inversion)  
";

%feature("docstring") clipper::Spacegroup::num_centering_symops "

get number of centering symops (inc identity)  
";

%feature("docstring") clipper::Spacegroup::symbol_hm "
";

%feature("docstring") clipper::Spacegroup::invariant_under_change_of_hand "

test if change of hand preserves spacegroup  

Test if hand-change is possible.  

Returns
-------
true if a change of hand preserves the spacegroup.  
";

%feature("docstring") clipper::Spacegroup::num_primitive_noninversion_symops "

get number of primitive non-inversion symops (inc identity)  
";

%feature("docstring") clipper::Spacegroup::product_op "

get symop number corresponding to the product of two symops  
";

%feature("docstring") clipper::Spacegroup::num_symops "

get number of symops  
";

%feature("docstring") clipper::Spacegroup::symbol_hall "
";

%feature("docstring") clipper::Spacegroup::Spacegroup "

null constructor  
";

%feature("docstring") clipper::Spacegroup::Spacegroup "

constructor: fast constructor for Null or P1 spacegroup  

Construct null or P1 spacegroup. This is faster than the normal constructor.  

Parameters
----------
* `type` :  
    Spacegroup::Null or Spacegroup::P1  
";

%feature("docstring") clipper::Spacegroup::Spacegroup "

constructor: from spacegroup description  

Construct a spacegroup and initialise with a spacegroup description.  

Parameters
----------
* `spgr_descr` :  
    The spacegroup description.  
";

%feature("docstring") clipper::Spacegroup::p1 "

Return P1 spacegroup.  
";

%feature("docstring") clipper::Spacegroup::spacegroup_number "
";

%feature("docstring") clipper::Spacegroup::is_null "

test if object has been initialised  

Returns
-------
true if the object has not been initalised.  
";

%feature("docstring") clipper::Spacegroup::num_inversion_symops "

get number of inversion symops (inc identity)  
";

%feature("docstring") clipper::Spacegroup::debug "
";

%feature("docstring") clipper::Spacegroup::asu_min "

get map ASU, lower bound  

The map ASU is an oblong which contains at least one assymetric unit. It is
guaranteed to be contained withing the unit box. The lower limit is always
0,0,0.  

Returns
-------
Fractional coordinate of the lower bound of the ASU.  
";

%feature("docstring") clipper::Spacegroup::init "

initialiser: from spacegroup description  

Initialise the spacegroup.  

Parameters
----------
* `spgr_descr` :  
    The spacegroup description.  
";

%feature("docstring") clipper::Spacegroup::symbol_laue "

return the Laue group symbol  

Returns
-------
The Laue group symbol. i.e. one of -1, 2/m, 2/mmm, -3, -3m, 4/m, 4/mmm, 6/m,
6/mmm, m-3, m-3m  
";

%feature("docstring") clipper::Spacegroup::asu_max "

get map ASU, upper bound  

The map ASU is an oblong which contains at least one assymetric unit. It is
guaranteed to be contained withing the unit box. The lower limit is always
0,0,0.  

Returns
-------
Fractional coordinate of the upper bound of the ASU.  
";

// File: classclipper_1_1Spgr__cacheobj.xml


%feature("docstring") clipper::Spgr_cacheobj "
";

%feature("docstring") clipper::Spgr_cacheobj::format "

string description  
";

%feature("docstring") clipper::Spgr_cacheobj::Spgr_cacheobj "

construct entry  
";

%feature("docstring") clipper::Spgr_cacheobj::matches "

compare entry  
";

// File: classclipper_1_1Spgr__descr.xml


%feature("docstring") clipper::Spgr_descr "

spacegroup description  

The spacegroup description is a compact description of a spacegroup. It may be
initialised from Hall or H-M symbols, a string of symops or a number. Internally
a hash code is used to refer to the spacegroup, so this object is only 32 bits
in size.  

For more details of spacegroup symbols, see Sydney R. Hall & Ralf W. Grosse-
Kunstleve 'Concise Space-Group Symbols',
http://www.kristall.ethz.ch/LFK/software/sginfo/hall_symbols.html  

C++ includes: spacegroup.h
";

%feature("docstring") clipper::Spgr_descr::set_preferred "

set preferred default spacegroup choice  

Sets the preferred origin or setting for initialising all Spgr_descr objects
using H-M symbols or Spacegroup numbers. cctbx uses origin choice '1' by
default, CCP4 uses '2'. Both packages use 'H' in preference to 'R'. Preferred
values are stored for both. Defaults are '1' and 'H'.  

CCP4 users may wish to add the following before using H-M codes or numbers.  

Parameters
----------
* `c` :  
    Either '1' or '2', 'H' or 'R'.  
";

%feature("docstring") clipper::Spgr_descr::symbol_hm "

return the H-M symbol  

The H-M symbol is only available if the spacegroup exists in the internal table,
see Hall & Grosse-Kunstleve.  

Returns
-------
The H-M symbol, or \"Unknown\" if unavailable.  
";

%feature("docstring") clipper::Spgr_descr::symbol_hm_ext "

return the extension H-M symbol  

The extension H-M symbol is only available if the spacegroup exists in the
internal table, see Hall & Grosse-Kunstleve.  

Returns
-------
The extension H-M symbol, or \"\"  
";

%feature("docstring") clipper::Spgr_descr::spacegroup_number "

return the spacegroup number  

The spacegroup number is only available if the spacegroup exists in the internal
table, see Hall & Grosse-Kunstleve.  

Returns
-------
The spacegroup number, or 0 if unavailable.  
";

%feature("docstring") clipper::Spgr_descr::symbol_xhm "

return the extended H-M symbol  

The extended H-M symbol is only available if the spacegroup exists in the
internal table, see Hall & Grosse-Kunstleve.  

Returns
-------
The extended H-M symbol, or \"Unknown\" if unavailable.  
";

%feature("docstring") clipper::Spgr_descr::hash "

return the hash code for the spacegroup  
";

%feature("docstring") clipper::Spgr_descr::symbol_hall "

return the Hall symbol  

The Hall symbol is only available if the spacegroup exists in the internal
table, see Hall & Grosse-Kunstleve.  

Returns
-------
The Hall symbol, or \"Unknown\" if unavailable.  
";

%feature("docstring") clipper::Spgr_descr::Spgr_descr "

null constructor  

Construct a null description spacegroup. The result is initialised to an invalid
spacegroup code.  
";

%feature("docstring") clipper::Spgr_descr::Spgr_descr "

constructor: from symbol or operators.  

Construct a spacegroup description from a text description, i.e. a symbol or
operators. This may be one of the following:  

*   Hall symbol, e.g. \" P 2ac 2ab\"  
*   H-M symbol, e.g. \"P 21 21 21\"  
*   Number, e.g. \"19\"  
*   List of symmetry operators separated by semicolons, e.g.
    \"x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2\"  

It is best to specify the type of symbol being used, however if this parameter
is omitted a guess will be made. Unfortunately, Hall and H-M symbols may be
ambiguous. Any ambiguity may be removed by appending parentheses \"()\" to the
end of the Hall symbol, otherwise the symbol will be interpreted as an H-M
symbol, and a Hall symbol if that fails.  

H-M symbols and spacegroup numbers may correspond to 2 different entries in
international tables. The choice between 2 origin settings or
hexagonal/rhomohedral settings is made using the set_preferred() method.  

Parameters
----------
* `name` :  
    The spacegroup symbol or operators.  
* `type` :  
    The type of symbol: Spacegroup::Symops, Spacegroup::Hall, Spacegroup::HM,
    Spacegroup::XHM, Spacegroup::Number  
";

%feature("docstring") clipper::Spgr_descr::Spgr_descr "

constructor: from number.  

See previous constuctor.  

Parameters
----------
* `num` :  
    The spacegroup number.  
";

%feature("docstring") clipper::Spgr_descr::Spgr_descr "

constructor: from symop list.  

This is not normally used, except in conjunction with Spgr_desc::generator_ops()
to derive one group from another.  
";

%feature("docstring") clipper::Spgr_descr::generator_ops "

return the generators for the spacegroup  
";

// File: classclipper_1_1String.xml


%feature("docstring") clipper::String "

String extension with simple parsing methods.  

String extension with primitive 'split' operation for parsing and pathname
processing operations.  

C++ includes: clipper_types.h
";

%feature("docstring") clipper::String::f64 "

convert to double  
";

%feature("docstring") clipper::String::f "

convert to ftype  
";

%feature("docstring") clipper::String::split "

String splitter - a very simple parser component.  
";

%feature("docstring") clipper::String::trim "

Return copy of string without leading and trailing blanks.  
";

%feature("docstring") clipper::String::head "

remove trailing path element  
";

%feature("docstring") clipper::String::f32 "

convert to float  
";

%feature("docstring") clipper::String::l "

convert to long  
";

%feature("docstring") clipper::String::tail "

get trailing path element  
";

%feature("docstring") clipper::String::i "

convert to int  
";

%feature("docstring") clipper::String::String "

null constructor  
";

%feature("docstring") clipper::String::String "

constructor: from string  
";

%feature("docstring") clipper::String::String "

constructor: from char*  
";

%feature("docstring") clipper::String::String "

constructor: from char*  
";

%feature("docstring") clipper::String::String "

constructor: from int  
";

%feature("docstring") clipper::String::String "

constructor: from long  
";

%feature("docstring") clipper::String::String "

constructor: from float  
";

%feature("docstring") clipper::String::String "

constructor: from double  
";

%feature("docstring") clipper::String::rational "

convert from rational to ftype  
";

%feature("docstring") clipper::String::rational "

construct string from rational f using base b  
";

%feature("docstring") clipper::String::nohead "

get leading path element  
";

%feature("docstring") clipper::String::notail "

remove leading path element  
";

// File: classclipper_1_1Symop.xml


%feature("docstring") clipper::Symop "

Crystallographic symmetry operator.  

This is identical to a fractional RTop, but has its own class since not all
fractional RTops are symops. For a full list of methods, see clipper::RTop and
clipper::RTop_frac  

C++ includes: symop.h
";

%feature("docstring") clipper::Symop::Symop "

null constructor  
";

%feature("docstring") clipper::Symop::Symop "

constructor: RTop  

Construct a symmetry operator and initialise it to the supplied RTop.
Translations are rounded to a basis of 48, and put on the range 0..1  

Parameters
----------
* `mat` :  
    The RTop to use.  
";

%feature("docstring") clipper::Symop::Symop "

constructor: from 4x4 matrix  

Construct a symmetry operator and initialise it to the supplied matrix.
Translations are rounded to a basis of 48, and put on the range 0..1  

Parameters
----------
* `mat` :  
    The 4x4 matrix to use. The [i][3] elements contain the translation.  
";

%feature("docstring") clipper::Symop::format "

return formatted String representation  

Return formatted representation of the symmetry operator.  

Returns
-------
The formatted text string, e.g. -x, -y+1/2, z.  
";

// File: classclipper_1_1Symop__code.xml


%feature("docstring") clipper::Symop_code "

Compressed encoded symmetry operator.  

This is a compresses representation of a crystallographic symmetry operator,
stored as a single 32-bit integer. It may be converted to or from a symop or an
int and compared, sorted, etc. The following guarantees are made concerning the
code:  

*   The identity operator has a code of zero.  
*   Operators with non-identity rotations will have higher codes than operators
    with identity rotations, for the same translation.  
*   Operators with non-zero translations will have higher codes than operators
    with zero translations, for the same rotation.  

C++ includes: symop.h
";

%feature("docstring") clipper::Symop_code::Symop_code "

null constructor  
";

%feature("docstring") clipper::Symop_code::Symop_code "

constructor: from int  
";

%feature("docstring") clipper::Symop_code::Symop_code "

constructor: from Symop  
";

%feature("docstring") clipper::Symop_code::Symop_code "

constructor: from Isymop  
";

%feature("docstring") clipper::Symop_code::symop "

convert to symop  

Construct a symmetry operator and initialise it to the matrix encoded in the
given int.  

Parameters
----------
* `code` :  
    The integer code.  
";

%feature("docstring") clipper::Symop_code::isymop "

convert to integerised symop  

Construct an integerised symmetry operator and initialise it to the matrix
encoded in the given int, with a grid (base) of (24,24,24).  

Parameters
----------
* `code` :  
    The integer code.  
";

%feature("docstring") clipper::Symop_code::init "

initialiser: from Isymop  
";

%feature("docstring") clipper::Symop_code::code_trn "

return code for translation part  
";

%feature("docstring") clipper::Symop_code::code_rot "

return code for rotation part  
";

%feature("docstring") clipper::Symop_code::identity "

identity code  
";

// File: classclipper_1_1Spgr__descr_1_1Symop__codes.xml


%feature("docstring") clipper::Spgr_descr::Symop_codes "

Vector of symop codes and associated methods.  

C++ includes: spacegroup.h
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::primitive_noninversion_ops "

return primitive non-inversion ops (by computation)  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::expand "

expand (incomplete) list of symops  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::inversion_ops "

return inversion ops (by computation)  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::pgrp_ops "

return point group ops  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::laue_ops "

return Laue ops  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::hash "

return hash code of symop list  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::primitive_ops "

return primitive incl inversion ops (by computation)  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::patterson_ops "

return Patterson ops  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::init_hall "

initialise from Hall symbol  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::centering_ops "

return lattice centering ops (by computation)  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::generator_ops "

return minimal list of generator ops  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::product "

return product of this (expanded) list by another (expanded) list  
";

%feature("docstring") clipper::Spgr_descr::Symop_codes::init_symops "

initialise from symops  
";

// File: classclipper_1_1TargetFn__base.xml


%feature("docstring") clipper::TargetFn_base "

abstract base class for least-squares resolution function target functions  

A target function must be able to return its value given the value of the basis
function for all HKL, and its derivative with respect the values of the basis
function for all HKL.  

Optionally, performance can be improved by returning a flag to indicate if the
target function is quadratic.  

C++ includes: resol_fn.h
";

%feature("docstring") clipper::TargetFn_base::rderiv "

return the value and derivatives of the target function  

If the value of f(h) is invalid, rderiv.r should be set to NaN  
";

%feature("docstring") clipper::TargetFn_base::debug "

test that the residuals, gradients, and curvatures are consistent  
";

%feature("docstring") clipper::TargetFn_base::type "

the type of the function: optionally used to improve convergence  
";

%feature("docstring") clipper::TargetFn_base::~TargetFn_base "

destructor  
";

// File: classclipper_1_1TargetFn__meanEnth.xml


%feature("docstring") clipper::TargetFn_meanEnth "

simple mean |E|n target  

This class implements the target function for calculating mean |E|n as a
function of position in reciprocal space. It includes the appropriate
multiplicity correction, and so can be applied to any type with an 'E' member
with the same dimensions as an |E| (or corrected |F| or |U|).  

This function should not be used to scale F's to E's. See TargetFn_scaleEsq.  

C++ includes: resol_targetfn.h
";

%feature("docstring") clipper::TargetFn_meanEnth::rderiv "

return the value and derivatives of the target function  
";

%feature("docstring") clipper::TargetFn_meanEnth::TargetFn_meanEnth "

constructor: takes the datalist against which to calc target, and power  
";

%feature("docstring") clipper::TargetFn_meanEnth::type "

the type of the function: optionally used to improve convergence  
";

// File: classclipper_1_1TargetFn__meanFnth.xml


%feature("docstring") clipper::TargetFn_meanFnth "

simple mean |F|n target  

This class implements the target function for calculating mean |F|n as a
function of position in reciprocal space. It includes the appropriate
multiplicity correction, and so can be applied to any type with an 'f' member
with the same dimensions as an |F| or |U| (or an uncorrected |E|).  

This function should not be used to scale F's to E's. See TargetFn_scaleEsq.  

C++ includes: resol_targetfn.h
";

%feature("docstring") clipper::TargetFn_meanFnth::TargetFn_meanFnth "

constructor: takes the datalist against which to calc target, and power  
";

%feature("docstring") clipper::TargetFn_meanFnth::rderiv "

return the value and derivatives of the target function  
";

%feature("docstring") clipper::TargetFn_meanFnth::type "

the type of the function: optionally used to improve convergence  
";

// File: classclipper_1_1TargetFn__scaleEsq.xml


%feature("docstring") clipper::TargetFn_scaleEsq "

|E|2 scaling target  

This class implements the target function for calculating the scale factor to
normalise to <|E|2> = 1. Note that this is not the same as dividing by <|E|2>,
except in a few special cases, e.g. a simple resolution bins calculation. The
resulting targen function is the square of the value by which |E| should be
multiplied to acheive the correct normalisation.  

C++ includes: resol_targetfn.h
";

%feature("docstring") clipper::TargetFn_scaleEsq::rderiv "

return the value and derivatives of the target function  
";

%feature("docstring") clipper::TargetFn_scaleEsq::type "

the type of the function: optionally used to improve convergence  
";

%feature("docstring") clipper::TargetFn_scaleEsq::TargetFn_scaleEsq "

constructor: takes the datalist against which to calc target  
";

// File: classclipper_1_1TargetFn__scaleF1F2.xml


%feature("docstring") clipper::TargetFn_scaleF1F2 "

|F|2 scaling target  

This class implements the target function for calculating the scale factor to
scale one set of F's to another. The resulting scale is the square of the factor
that scales the first set of data to match the second.  

C++ includes: resol_targetfn.h
";

%feature("docstring") clipper::TargetFn_scaleF1F2::rderiv "

return the value and derivatives of the target function  
";

%feature("docstring") clipper::TargetFn_scaleF1F2::TargetFn_scaleF1F2 "

constructor: takes the datalist against which to calc target  
";

%feature("docstring") clipper::TargetFn_scaleF1F2::type "

the type of the function: optionally used to improve convergence  
";

// File: classclipper_1_1TargetFn__scaleI1I2.xml


%feature("docstring") clipper::TargetFn_scaleI1I2 "

This class implements the target function for calculating the scale factor to
scale one set of I's to another. The resulting scale is the square of the factor
that scales the first set of data to match the second.  

C++ includes: resol_targetfn.h
";

%feature("docstring") clipper::TargetFn_scaleI1I2::rderiv "

return the value and derivatives of the target function  
";

%feature("docstring") clipper::TargetFn_scaleI1I2::TargetFn_scaleI1I2 "

constructor: takes the datalist against which to calc target  
";

%feature("docstring") clipper::TargetFn_scaleI1I2::type "

the type of the function: optionally used to improve convergence  
";

// File: classclipper_1_1TargetFn__scaleLogF1F2.xml


%feature("docstring") clipper::TargetFn_scaleLogF1F2 "

log |F|2 scaling target  

This class implements the target function for calculating the scale factor to
scale the weighted log of one set of F's to another. The resulting scale is the
square of the factor that scales the first set of data to match the second. The
log scaling target is used in conjunction with the log-Gaussian basis functions
for a fast and robust approximation to iso/aniso Gaussian scaling.  

C++ includes: resol_targetfn.h
";

%feature("docstring") clipper::TargetFn_scaleLogF1F2::rderiv "

return the value and derivatives of the target function  
";

%feature("docstring") clipper::TargetFn_scaleLogF1F2::type "

the type of the function: optionally used to improve convergence  
";

%feature("docstring") clipper::TargetFn_scaleLogF1F2::TargetFn_scaleLogF1F2 "

constructor: takes the datalist against which to calc target  
";

// File: classclipper_1_1TargetFn__scaleLogI1I2.xml


%feature("docstring") clipper::TargetFn_scaleLogI1I2 "

log |I| scaling target  

This class implements the target function for calculating the scale factor to
scale the weighted log of one set of I's to another. The resulting scale is the
square of the factor that scales the first set of data to match the second. The
log scaling target is used in conjunction with the log-Gaussian basis functions
for a fast and robust approximation to iso/aniso Gaussian scaling.  

C++ includes: resol_targetfn.h
";

%feature("docstring") clipper::TargetFn_scaleLogI1I2::rderiv "

return the value and derivatives of the target function  
";

%feature("docstring") clipper::TargetFn_scaleLogI1I2::type "

the type of the function: optionally used to improve convergence  
";

%feature("docstring") clipper::TargetFn_scaleLogI1I2::TargetFn_scaleLogI1I2 "

constructor: takes the datalist against which to calc target  
";

// File: classclipper_1_1TargetFn__sigmaa.xml


%feature("docstring") clipper::TargetFn_sigmaa "

Deprecated
simple sigma_a target function  

par: Warning: Convergence of this basis-function can be
    unreliable under some circumstances. Use clipper::TargetFn_sigmaa_omegaa
    instead, except for development purposes.  

This class implements the target function for calculating sigma_a. Required is a
datalist containing Eo, Ec.  

This version simplifies terms in |Eo|^2 and |Ec|^2 which should average out to 1
if the normalisation scheme is consistent with the sigmaa calc.  

C++ includes: resol_targetfn.h
";

%feature("docstring") clipper::TargetFn_sigmaa::TargetFn_sigmaa "

constructor: takes the datalist against which to calc target  
";

%feature("docstring") clipper::TargetFn_sigmaa::rderiv "

return the value and derivatives of the target function  
";

%feature("docstring") clipper::TargetFn_sigmaa::sigmaa "

convert function to sigmaa  
";

// File: classclipper_1_1TargetFn__sigmaa__omegaa.xml


%feature("docstring") clipper::TargetFn_sigmaa_omegaa "

Deprecated
simple sigma_a target function  

This class implements the target function for calculating sigma_a. Required is a
datalist containing Eo, Ec.  

It actually refines omegaa = sigmaa/(1-sigmaa^2). This has better proerties for
refinement. To get sigmaa use  This is available as a static function:  

This version simplifies terms in |Eo|^2 and |Ec|^2 which should average out to 1
if the normalisation scheme is consistent with the sigmaa calc.  

Convergence is good for calculations using the 'binner' basis function, however
the smooth basis function have convergence problems. This is still under
investigation.  

C++ includes: resol_targetfn.h
";

%feature("docstring") clipper::TargetFn_sigmaa_omegaa::rderiv "

return the value and derivatives of the target function  
";

%feature("docstring") clipper::TargetFn_sigmaa_omegaa::TargetFn_sigmaa_omegaa "

constructor: takes the datalist against which to calc target  
";

%feature("docstring") clipper::TargetFn_sigmaa_omegaa::sigmaa "

convert omegaa to sigmaa  
";

// File: classclipper_1_1Test__base.xml


%feature("docstring") clipper::Test_base "

Base class for clipper self-test classes.  

Override this class with test classes for each clipper package.  

C++ includes: clipper_test.h
";

%feature("docstring") clipper::Test_base::set_stream "
";

%feature("docstring") clipper::Test_base::Test_base "
";

// File: classclipper_1_1Test__core.xml


%feature("docstring") clipper::Test_core "
";

// File: classclipper_1_1data_1_1Test__data.xml


%feature("docstring") clipper::data::Test_data "

Class to return test data.  

C++ includes: test_data.h
";

%feature("docstring") clipper::data::Test_data::hkl_data_abcd "

Return HKL_data class.  
";

%feature("docstring") clipper::data::Test_data::atom_list "

Return atom list.  
";

%feature("docstring") clipper::data::Test_data::hkl_data_f_sigf "

Return HKL_data class.  
";

%feature("docstring") clipper::data::Test_data::Test_data "

Null constructor: fills the arrays.  
";

// File: structclipper_1_1data_1_1TESThkldata.xml


%feature("docstring") clipper::data::TESThkldata "
";

// File: structclipper_1_1data_1_1TESTxyzdata.xml


%feature("docstring") clipper::data::TESTxyzdata "
";

// File: classclipper_1_1Thread__base.xml


%feature("docstring") clipper::Thread_base "

Thread base class: Override this to create new threads.  

      To create a thread, override this class. Store data as members
    with accessors to set input and read output. Override the Run()
    method to do the actual work. e.g. the following class implements
    a thread which can sum a list of numbers.
  

C++ includes: clipper_thread.h
";

%feature("docstring") clipper::Thread_base::id "
";

%feature("docstring") clipper::Thread_base::~Thread_base "
";

%feature("docstring") clipper::Thread_base::Thread_base "
";

%feature("docstring") clipper::Thread_base::unlock "
";

%feature("docstring") clipper::Thread_base::join "
";

%feature("docstring") clipper::Thread_base::run "
";

%feature("docstring") clipper::Thread_base::lock "
";

// File: unionclipper_1_1Util_1_1U32.xml

// File: unionclipper_1_1Util_1_1U64.xml

// File: classclipper_1_1U__aniso__frac.xml


%feature("docstring") clipper::U_aniso_frac "

Anisotropic fractional atomic displacement parameters.  

These are defined on fractional atomic coordinates in A-2, i.e. they are
anisotropic U values.  

C++ includes: coords.h
";

%feature("docstring") clipper::U_aniso_frac::U_aniso_frac "

null constructor  
";

%feature("docstring") clipper::U_aniso_frac::U_aniso_frac "

constructor: from Mat33sym  
";

%feature("docstring") clipper::U_aniso_frac::U_aniso_frac "

constructor: from Uij  
";

%feature("docstring") clipper::U_aniso_frac::transform "

return transformed U_aniso  

The aniso U is transformed by the given RT op.  

Parameters
----------
* `u` :  
    The aniso U.  
";

%feature("docstring") clipper::U_aniso_frac::u_aniso_orth "

fractional-orthogonal conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed coordinate.  
";

// File: classclipper_1_1U__aniso__orth.xml


%feature("docstring") clipper::U_aniso_orth "

Anisotropic orthogonal atomic displacement parameters.  

These are defined on orthogonal atomic coordinates in A-2, i.e. they are
anisotropic U values.  

C++ includes: coords.h
";

%feature("docstring") clipper::U_aniso_orth::u_iso "

return nearest isotropic U  

The best isotropic U is the cube root of the determinant of the matrix of
anisotropic coefficients. NOTE: This is not the conventional definition, but the
mathematically correct one, and gives a better approximation to the anisotropic
U (i.e. lower R-factors).  

Returns
-------
The nearest isotropic U.  
";

%feature("docstring") clipper::U_aniso_orth::U_aniso_orth "

null constructor  
";

%feature("docstring") clipper::U_aniso_orth::U_aniso_orth "

constructor: from Mat33sym  
";

%feature("docstring") clipper::U_aniso_orth::U_aniso_orth "

constructor: from isotropic U  
";

%feature("docstring") clipper::U_aniso_orth::U_aniso_orth "

constructor: from Uij  
";

%feature("docstring") clipper::U_aniso_orth::transform "

return transformed U_aniso  

The aniso U is transformed by the given RT op.  

Parameters
----------
* `u` :  
    The aniso U.  
";

%feature("docstring") clipper::U_aniso_orth::u_aniso_frac "

orthogonal-fractional conversion  

Parameters
----------
* `cell` :  
    The cell concerned  

Returns
-------
The transformed coordinate.  
";

// File: classclipper_1_1Util.xml


%feature("docstring") clipper::Util "

Utility class.  

This class holds a set of useful static functions and values. You should never
need to instantiate this class: Refer to members using Util::  

C++ includes: clipper_util.h
";

%feature("docstring") clipper::Util::set_null "

set null generic value  
";

%feature("docstring") clipper::Util::set_null "

set null floating value - a specific value of NaN used for missings  
";

%feature("docstring") clipper::Util::set_null "

set null floating value - a specific value of NaN used for missings  
";

%feature("docstring") clipper::Util::invsim "

Inverse Sim function: I1(X)/I0(X)  

Parameters
----------
* `x` :  
    I1(y)/I0(y)  

Returns
-------
y  
";

%feature("docstring") clipper::Util::sim_integ "

Integral of Sim function: log(I0(X))  
";

%feature("docstring") clipper::Util::min "

min  
";

%feature("docstring") clipper::Util::is_nan "

fast Util::nan() test  

Used for missing entries: THIS DOES NOT DISTINGUISH BETWEEN NAN & INF  
";

%feature("docstring") clipper::Util::is_nan "

fast Util::nan() test  

Used for missing entries: THIS DOES NOT DISTINGUISH BETWEEN NAN & INF  
";

%feature("docstring") clipper::Util::sqr "

square  
";

%feature("docstring") clipper::Util::mod "

Corrected mod.  
";

%feature("docstring") clipper::Util::mod "

Corrected mod.  
";

%feature("docstring") clipper::Util::nand "

fast Util::nan() value  
";

%feature("docstring") clipper::Util::d2rad "

degree-to-radian conversion  

Parameters
----------
* `x` :  
    Angle in degrees  

Returns
-------
Angle in radians  
";

%feature("docstring") clipper::Util::nanf "

fast Util::nan() value  
";

%feature("docstring") clipper::Util::bound "

bound a value by limits  
";

%feature("docstring") clipper::Util::twopi "

2 pi  
";

%feature("docstring") clipper::Util::twopi2 "

2 pi squared  
";

%feature("docstring") clipper::Util::mean "

Convert F+/F- to mean F, with NaN checks.  
";

%feature("docstring") clipper::Util::sim_deriv "

Derivative of Sim function: d/dx( I1(X)/I0(x) )  
";

%feature("docstring") clipper::Util::isnan "

slow general NaN test for compatibility  

Works for all architectures with IEEE arithmetic only  
";

%feature("docstring") clipper::Util::isnan "

slow general NaN test for compatibility  

Works for all architectures with IEEE arithmetic only  
";

%feature("docstring") clipper::Util::intc "

Truncate-to-integer above: int(ceil(a))  
";

%feature("docstring") clipper::Util::isqrt "

Integer square root (returns floor of sqrt)  
";

%feature("docstring") clipper::Util::sim_deriv_recur "

Derivative of Sim function using recurrance: -sim(x)/x + (1 - sim(x)^2)  
";

%feature("docstring") clipper::Util::u2b "

Convert isotropic U-value to B-factor.  
";

%feature("docstring") clipper::Util::Util "

null constructor  
";

%feature("docstring") clipper::Util::sig_mean "

Convert sigF+/sigF-/cov to sig F, with NaN checks.  
";

%feature("docstring") clipper::Util::sim "

Sim function: I1(X)/I0(X)  

Parameters
----------
* `x` :  
    The argument  

Returns
-------
I1(x)/I0(x)  
";

%feature("docstring") clipper::Util::intf "

Truncate-to-integer: int(floor(a))  
";

%feature("docstring") clipper::Util::eightpi2 "

8 pi squared  
";

%feature("docstring") clipper::Util::intr "

Round-to-integer: int(round(a))  
";

%feature("docstring") clipper::Util::b2u "

Convert isotropic B-factor to U-value.  
";

%feature("docstring") clipper::Util::bessel_i0 "

Modified Bessel function of the first kind.  
";

%feature("docstring") clipper::Util::nan "

fast Util::nan() value  
";

%feature("docstring") clipper::Util::max "

max  
";

%feature("docstring") clipper::Util::pi "

pi  
";

%feature("docstring") clipper::Util::is_null "

fast test for null floating value - only works if set from Util::null()  
";

%feature("docstring") clipper::Util::is_null "

fast test for null floating value - only works if set from Util::null()  
";

%feature("docstring") clipper::Util::rad2d "

degree-to-radian conversion  

Parameters
----------
* `x` :  
    Angle in radians  

Returns
-------
Angle in degrees  
";

%feature("docstring") clipper::Util::swap "

swap the contents of two objects  
";

%feature("docstring") clipper::Util::swap "

swap the contents of two objects, using third as store (for speed)  
";

%feature("docstring") clipper::Util::atanh "

Arc hyperbolic tangent.  
";

// File: classclipper_1_1Vec3.xml


%feature("docstring") clipper::Vec3 "

3-vector class  

C++ includes: clipper_types.h
";

%feature("docstring") clipper::Vec3::is_null "

test for null vector  
";

%feature("docstring") clipper::Vec3::format "

return formatted String representation  
";

%feature("docstring") clipper::Vec3::Vec3 "

null constructor  
";

%feature("docstring") clipper::Vec3::Vec3 "

constructor: from individual values  
";

%feature("docstring") clipper::Vec3::Vec3 "

constructor: copy/convert  
";

%feature("docstring") clipper::Vec3::cross "

Vector cross product.  
";

%feature("docstring") clipper::Vec3::dot "

Vector dot product (equivalent to *)  
";

%feature("docstring") clipper::Vec3::null "

return null vector (only valid for floating point types)  
";

%feature("docstring") clipper::Vec3::zero "

return zero vector  
";

%feature("docstring") clipper::Vec3::unit "

return unit vector with same direction as this vector  
";

%feature("docstring") clipper::Vec3::equals "

test equality  
";

// File: classclipper_1_1Xmap.xml


%feature("docstring") clipper::Xmap "

Xmap<T>: actual crystallographic map class.  

The crystallographic map class stores a map of arbitrary data type. Its main
difference from a 3-d array is that the data extent appears to be infinite, and
yet internally only a unique ASU is stored. Iterators provide efficient access
to data.  

This is derived from Xmap_base, and adds the templatised data itself and the
methods which deal with it.  

note: The following methods are inherited from Xmap_base but are documented here
    for convenience: cell(), spacegroup(), grid_sampling(), grid_asu(),
    in_asu(), multiplicity(), to_map_unit(), first(), first_coord().  

C++ includes: xmap.h
";

%feature("docstring") clipper::Xmap::Xmap "

Null constructor, for later initialisation.  
";

%feature("docstring") clipper::Xmap::Xmap "

constructor: from spacegroup, cell, and grid  
";

%feature("docstring") clipper::Xmap::coord_of "
";

%feature("docstring") clipper::Xmap::set_data "

set a density value for an arbitrary position  

Accessing the data by coordinate, rather than by Map_reference_index or
Map_reference_coord, involves a symmetry lookup and is therefore slow. Avoid
using these methods when you can.  
";

%feature("docstring") clipper::Xmap::set_data "

set data by index (not recommended)  

Accessing the data by index, rather than by Map_reference_index or
Map_reference_coord, is generally to be avoided since the indices do not start
at zero and do not increase contiguously. These methods are only useful when a
large number of references into a map must be stored, e.g. for sorting into
density order.  

Returns
-------
true if data was set, i.e. index is valid.  
";

%feature("docstring") clipper::Xmap::spacegroup "
";

%feature("docstring") clipper::Xmap::interp_grad "

get map value and grad for fractional coord using supplied interpolator  

The value of the map at the desired non-grid fractional coordinate and its
gradient are calculated using the supplied interpolator template.  

Parameters
----------
* `pos` :  
    The fractional coord at which the density is to be calcuated.  
* `val` :  
    The value of the density at that point.  
* `grad` :  
    The interpolated gradient vector with respect to the fractional coordinates
    (see Cell::coord_orth).  
";

%feature("docstring") clipper::Xmap::interp_grad "

get map value and grad for map coord using supplied interpolator  

The value of the map at the desired non-grid map coordinate and its gradient are
calculated using the supplied interpolator template.  

Parameters
----------
* `pos` :  
    The map coord at which the density is to be calcuated.  
* `val` :  
    The value of the density at that point.  
* `grad` :  
    The interpolated gradient vector with respect to the map coordinates (see
    Cell::coord_orth).  
";

%feature("docstring") clipper::Xmap::interp "

get map value for fractional coord using supplied interpolator  

The value of the map at the desired non-grid fractional coordinate are
calculated using the supplied interpolator template.  

Parameters
----------
* `pos` :  
    The fractional coord at which the density is to be calcuated.  

Returns
-------
The value of the density at that point.  
";

%feature("docstring") clipper::Xmap::interp "

get map value for map coord using supplied interpolator  

The value of the map at the desired non-grid map coordinate are calculated using
the supplied interpolator template.  

Parameters
----------
* `pos` :  
    The map coord at which the density is to be calcuated.  

Returns
-------
The value of the density at that point.  
";

%feature("docstring") clipper::Xmap::grid_sampling "
";

%feature("docstring") clipper::Xmap::fft_from "

FFT from reflection list to map.  

An FFT is calculated using the provided reflection list of F_phi, and used to
fill this map. The reflection list is unchanged.  

Parameters
----------
* `fphidata` :  
    The reflection data list to use  
";

%feature("docstring") clipper::Xmap::first_coord "
";

%feature("docstring") clipper::Xmap::fft_to "

FFT from map to reflection list.  

The Fourier transform of this map is calculated and used to fill a reflection
list of F_phi. The map is unchanged.  

Arguably this should be part of hkl_data<F_phi<T>>. But that requires writing a
specialisation of hkl_data for F_phi. This is simpler and imposes less demands
on the compiler.  

Parameters
----------
* `fphidata` :  
    The reflection data list to set.  
";

%feature("docstring") clipper::Xmap::interp_curv "

get map value and curv for fractional coord using supplied interpolator  

The value of the map at the desired non-grid fractional coordinate and its
gradient and curvature are calculated using the supplied interpolator template.
e.g.  

Parameters
----------
* `pos` :  
    The fractional coord at which the density is to be calcuated.  
* `val` :  
    The value of the density at that point.  
* `grad` :  
    The interpolated gradient vector with respect to the fractional coordinates
    (see Cell::coord_orth).  
* `curv` :  
    The interpolated curvature matrix with respect to the fractional coordinates
    (see Cell::coord_orth).  
";

%feature("docstring") clipper::Xmap::interp_curv "

get map value and curv for map coord using supplied interpolator  

The value of the map at the desired non-grid map coordinate and its gradient and
curvature are calculated using the supplied interpolator template. e.g.  

Parameters
----------
* `pos` :  
    The map coord at which the density is to be calcuated.  
* `val` :  
    The value of the density at that point.  
* `grad` :  
    The interpolated gradient vector with respect to the map coordinates (see
    Cell::coord_orth).  
* `curv` :  
    The interpolated curvature matrix with respect to the map coordinates (see
    Cell::coord_orth).  
";

%feature("docstring") clipper::Xmap::index_of "
";

%feature("docstring") clipper::Xmap::multiplicity "
";

%feature("docstring") clipper::Xmap::to_map_unit "
";

%feature("docstring") clipper::Xmap::first "
";

%feature("docstring") clipper::Xmap::init "

initialiser: from spacegroup, cell, and grid  
";

%feature("docstring") clipper::Xmap::get_data "

get a density value for an arbitrary position  

Accessing the data by coordinate, rather than by Map_reference_index or
Map_reference_coord, involves a symmetry lookup and is therefore slow. Avoid
using these methods when you can.  
";

%feature("docstring") clipper::Xmap::get_data "

get data by index (not recommended)  

Accessing the data by index, rather than by Map_reference_index or
Map_reference_coord, is generally to be avoided since the indices do not start
at zero and do not increase contiguously. These methods are only useful when a
large number of references into a map must be stored, e.g. for sorting into
density order.  
";

%feature("docstring") clipper::Xmap::cell "
";

%feature("docstring") clipper::Xmap::grid_asu "
";

// File: classclipper_1_1Xmap__base.xml


%feature("docstring") clipper::Xmap_base "

Xmap_base: base for crystallographic map class.  

The crystallographic map class stores a map of arbitrary data type. Its main
difference from a 3-d array is that the data extent appears to be infinite, and
yet internally only a unique ASU is stored. Iterators provide efficient access
to data.  

This base contains everything except the data, which is templated in the derived
type Xmap<T>  

C++ includes: xmap.h
";

%feature("docstring") clipper::Xmap_base::grid_sampling "

get the cell grid  
";

%feature("docstring") clipper::Xmap_base::spacegroup "

get the spacegroup  
";

%feature("docstring") clipper::Xmap_base::coord_map "

convert orthogonal coordinate to map  

Parameters
----------
* `co` :  
    The orthogonal coordinate to be converted.  

Returns
-------
The equivalent grid coordinate.  
";

%feature("docstring") clipper::Xmap_base::Xmap_base::Map_reference_index "
";

%feature("docstring") clipper::Xmap_base::Xmap_base::Map_reference_base "
";

%feature("docstring") clipper::Xmap_base::grid_asu "

get the ASU grid  
";

%feature("docstring") clipper::Xmap_base::default_type "

set/get default backend type  
";

%feature("docstring") clipper::Xmap_base::first "

return a Map_reference_index for this map  
";

%feature("docstring") clipper::Xmap_base::to_map_unit "

function to pick right cell repeat for any grid coord  
";

%feature("docstring") clipper::Xmap_base::coord_orth "

convert map coordinate to orthogonal  

Parameters
----------
* `cm` :  
    The grid coordinate to be converted.  

Returns
-------
The equivalent orthogonal coordinate.  
";

%feature("docstring") clipper::Xmap_base::multiplicity "

get multiplicity of a map grid point  

The multiplicity is the number of times the spacegroup operators map a
particular grid point onto itself. This is required in order to properly weight
map statistics so as to get the same result from just an ASU as using the whole
cell.  

Parameters
----------
* `pos` :  
    The coordinate of the grid point.  

Returns
-------
The multiplicty of the point.  
";

%feature("docstring") clipper::Xmap_base::coord_of "

map coordinate from index  

Parameters
----------
* `index` :  
    The index.  

Returns
-------
The corresponding grid coordinate.  
";

%feature("docstring") clipper::Xmap_base::first_coord "

return a Map_reference_coord for this map  
";

%feature("docstring") clipper::Xmap_base::in_map "

(This method is for compatibility with NXmap - it always returns true)  
";

%feature("docstring") clipper::Xmap_base::in_map "

(This method is for compatibility with NXmap - it always returns true)  
";

%feature("docstring") clipper::Xmap_base::is_null "

test if object has been initialised  

Returns
-------
true if the object has not been initalised.  
";

%feature("docstring") clipper::Xmap_base::cell "

get the cell  
";

%feature("docstring") clipper::Xmap_base::Xmap_base::Map_reference_coord "
";

%feature("docstring") clipper::Xmap_base::index_of "

map index from coordinate  

This does not check symmetry equivalents.  

Parameters
----------
* `coord` :  
    The coordinate.  

Returns
-------
The index, or -1 if it does not exist.  
";

// File: classclipper_1_1Xmap__cacheobj.xml


%feature("docstring") clipper::Xmap_cacheobj "
";

%feature("docstring") clipper::Xmap_cacheobj::Xmap_cacheobj "

construct entry  
";

%feature("docstring") clipper::Xmap_cacheobj::matches "

compare entry  
";

%feature("docstring") clipper::Xmap_cacheobj::format "

string description  
";

// File: namespaceclipper.xml

%feature("docstring") clipper::datatypes::reci_asu "
";

%feature("docstring") clipper::datatypes::real_asu "
";

// File: namespaceclipper_1_1data.xml

%feature("docstring") clipper::data::ASU_223 "
";

%feature("docstring") clipper::data::ASU_32D "
";

%feature("docstring") clipper::data::ASU_32B "
";

%feature("docstring") clipper::data::ASU_32C "
";

%feature("docstring") clipper::data::ASU_232 "
";

%feature("docstring") clipper::data::ASU_T11 "
";

%feature("docstring") clipper::data::ASU_121 "
";

%feature("docstring") clipper::data::ASU_32V "
";

%feature("docstring") clipper::data::ASU_32W "
";

%feature("docstring") clipper::data::ASU_32U "
";

%feature("docstring") clipper::data::ASU_32Z "
";

%feature("docstring") clipper::data::ASU_32X "
";

%feature("docstring") clipper::data::ASU_1T1 "
";

%feature("docstring") clipper::data::ASU_22W "
";

%feature("docstring") clipper::data::ASU_22V "
";

%feature("docstring") clipper::data::ASU_22U "
";

%feature("docstring") clipper::data::ASU_422 "
";

%feature("docstring") clipper::data::ASU_11T "
";

%feature("docstring") clipper::data::ASU_32A "
";

%feature("docstring") clipper::data::ASU_M3B "
";

%feature("docstring") clipper::data::ASU_M3M "
";

%feature("docstring") clipper::data::ASU_211 "
";

%feature("docstring") clipper::data::ASU_311 "
";

%feature("docstring") clipper::data::ASU_141 "
";

%feature("docstring") clipper::data::ASU_111 "
";

%feature("docstring") clipper::data::ASU_222 "
";

%feature("docstring") clipper::data::ASU_113 "
";

%feature("docstring") clipper::data::ASU_112 "
";

%feature("docstring") clipper::data::ASU_114 "
";

%feature("docstring") clipper::data::ASU_224 "
";

%feature("docstring") clipper::data::ASU_131 "
";

%feature("docstring") clipper::data::ASU_31D "
";

%feature("docstring") clipper::data::ASU_31C "
";

%feature("docstring") clipper::data::ASU_31B "
";

%feature("docstring") clipper::data::ASU_31A "
";

%feature("docstring") clipper::data::ASU_21Z "
";

%feature("docstring") clipper::data::ASU_21X "
";

%feature("docstring") clipper::data::ASU_21Y "
";

%feature("docstring") clipper::data::ASU_21V "
";

%feature("docstring") clipper::data::ASU_21W "
";

%feature("docstring") clipper::data::ASU_21U "
";

%feature("docstring") clipper::data::ASU_32Y "
";

%feature("docstring") clipper::data::ASU_411 "
";

%feature("docstring") clipper::data::ASU_322 "
";

%feature("docstring") clipper::data::ASU_242 "
";

// File: namespaceclipper_1_1data32.xml

// File: namespaceclipper_1_1data64.xml

// File: namespaceclipper_1_1datatypes.xml

// File: atomsf_8cpp.xml

// File: atomsf_8h.xml

// File: cell_8cpp.xml

// File: cell_8h.xml

// File: class__overview_8dox.xml

// File: clipper_8dox.xml

// File: clipper__contrib_8dox.xml

// File: clipper__instance_8cpp.xml

// File: clipper__instance_8h.xml

// File: clipper__memory_8cpp.xml

// File: clipper__memory_8h.xml

// File: clipper__message_8cpp.xml

// File: clipper__message_8h.xml

// File: clipper__minimol_8dox.xml

// File: clipper__mmdb_8dox.xml

// File: clipper__mmdbold_8dox.xml

// File: clipper__precision_8h.xml

// File: clipper__stats_8cpp.xml

// File: clipper__stats_8h.xml

// File: clipper__sysdep_8h.xml

// File: clipper__test_8cpp.xml

// File: clipper__test_8h.xml

// File: clipper__thread_8cpp.xml

// File: clipper__thread_8h.xml

// File: clipper__types_8cpp.xml

// File: clipper__types_8h.xml

// File: clipper__util_8cpp.xml

// File: clipper__util_8h.xml

// File: container_8cpp.xml

// File: container_8h.xml

// File: container__hkl_8cpp.xml

// File: container__hkl_8h.xml

// File: container__map_8cpp.xml

// File: container__map_8h.xml

// File: container__types_8cpp.xml

// File: container__types_8h.xml

// File: conventions_8dox.xml

// File: coords_8cpp.xml

// File: coords_8h.xml

// File: coordtypes_8dox.xml

// File: derivs_8cpp.xml

// File: derivs_8h.xml

// File: develop_8dox.xml

// File: develop__hkl_8dox.xml

// File: develop__map_8dox.xml

// File: develop__model_8dox.xml

// File: fftmap_8cpp.xml

// File: fftmap_8h.xml

// File: fftmap__sparse_8cpp.xml

// File: fftmap__sparse_8h.xml

// File: hkl__compute_8cpp.xml

// File: hkl__compute_8h.xml

// File: hkl__data_8cpp.xml

// File: hkl__data_8h.xml

// File: hkl__datatypes_8cpp.xml

// File: hkl__datatypes_8h.xml

// File: hkl__info_8cpp.xml

// File: hkl__info_8h.xml

// File: hkl__lookup_8cpp.xml

// File: hkl__lookup_8h.xml

// File: hkl__operators_8cpp.xml

// File: hkl__operators_8h.xml

// File: installation_8dox.xml

// File: map__interp_8cpp.xml

// File: map__interp_8h.xml

// File: map__utils_8cpp.xml

// File: map__utils_8h.xml

// File: nxmap_8cpp.xml

// File: nxmap_8h.xml

// File: nxmap__operator_8cpp.xml

// File: nxmap__operator_8h.xml

// File: ramachandran_8cpp.xml

// File: ramachandran_8h.xml

// File: resol__basisfn_8cpp.xml

// File: resol__basisfn_8h.xml

// File: resol__fn_8cpp.xml

// File: resol__fn_8h.xml

// File: resol__targetfn_8cpp.xml

// File: resol__targetfn_8h.xml

// File: rotation_8cpp.xml

// File: rotation_8h.xml

// File: spacegroup_8cpp.xml

// File: spacegroup_8h.xml

// File: spacegroup__data_8cpp.xml

// File: spacegroup__data_8h.xml

// File: symop_8cpp.xml

// File: symop_8h.xml

// File: test__core_8cpp.xml

// File: test__core_8h.xml

// File: test__data_8cpp.xml

// File: test__data_8h.xml

// File: wheretolook_8dox.xml

// File: xmap_8cpp.xml

// File: xmap_8h.xml

// File: p_class_overview.xml

// File: p_conventions.xml

// File: p_coords.xml

// File: p_develop.xml

// File: p_develop_hkl.xml

// File: p_develop_map.xml

// File: p_develop_model.xml

// File: p_installation.xml

// File: p_wheretolook.xml

// File: deprecated.xml

// File: dir_a214a1b60c6c26183dabe05fa6aecf4a.xml

// File: dir_85cb9eaae0d3200e6e62b752294ff8ec.xml

// File: indexpage.xml

