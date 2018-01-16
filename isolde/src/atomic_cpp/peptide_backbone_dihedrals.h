#ifndef ISOLDE_BACKBONE_DIHEDRALS
#define ISOLDE_BACKBONE_DIHEDRALS

#define PYINSTANCE_EXPORT
#include <atomstruct/destruct.h>
#include <atomstruct/AtomicStructure.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Coord.h>
#include <atomstruct/Residue.h>
#include <pyinstance/PythonInstance.declare.h>

#include "dihedral.h"

using namespace atomstruct;

namespace isolde 
{
//! Encapsulates phi, psi and omega dihedrals for an amino acid residue
/*! Each amino acid residue inside a polypeptide chain has the following
 * three backbone dihedrals:
 * omega: Ca(-1)-C(-1)-N-Ca
 * phi: C(-1)-N-Ca-C
 * psi: N-Ca-C-N(+1)
 * The omega torsion defines the conformation (cis or trans) of the 
 * peptide bond connecting this residue to the previous one. Phi and psi
 * together define the "twist" of the backbone, and are commonly plotted
 * against each other on the Ramachandran plot. N-terminal residues have
 * no omega or phi dihedrals, while C-terminal residues lack psi.
 */ 
class Peptide_Backbone_Dihedrals: public DestructionObserver, public pyinstance::PythonInstance<Peptide_Backbone_Dihedrals>
{
public:
    //! Constructor for an internal residue
    Peptide_Backbone_Dihedrals(const Residue &N_residue, const Residue &this_residue, const Residue &C_residue);
    //! Constructor for a N-terminal residue
    Peptide_Backbone_Dihedrals(const Residue &this_residue, const Residue &C_residue);
    //! Constructor for a C-terminal residue
    Peptide_Backbone_Dihedrals(const Residue &N_residue, const Residue &this_residue);

private:
    bool _has_omega;
    bool _has_phi;
    bool _has_psi;
    Proper_Dihedral* _omega;
    Proper_Dihedral* _phi;
    Proper_Dihedral* _psi;

}; // class Peptide_Backbone_Dihedrals



] //namespace isolde
#endif //ISOLDE_BACKBONE_DIHEDRALS
