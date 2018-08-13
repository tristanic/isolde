#pragma once
#include <unordered_map>
#include <clipper/clipper.h>

//! Van der Waals radii for use in map calculations. Mostly as per CCTBX,
//  except that metals use the radius of the ion corresponding to their PDB
//  residue name.
namespace clipper_cx
{
namespace data
{
    static const std::unordered_map<std::string, clipper::ftype32> vdw_radii = {
        {"H", 1.20},
        {"He", 1.40},
        {"Li", 0.9}, // Li+
        {"Be", 0.63},
        {"B", 1.75},
        {"C", 1.77},
        {"N", 1.50},
        {"O", 1.45},
        {"F", 1.47},
        {"Ne", 1.54},
        {"Na", 1.16}, // Na+
        {"Mg", 1.03}, // Mg2+
        {"Al", 0.675},
        {"Si", 2.10},
        {"P", 1.90},
        {"S", 1.80},
        {"Cl", 1.67}, // Cl-
        {"Ar", 1.88},
        {"K", 1.65}, // K+
        {"Ca", 1.26}, // Ca2+
        {"Sc", 1.32},
        {"Ti", 1.95},
        {"V", 1.06},
        {"Cr", 1.13},
        {"Mn", 1.19},
        {"Fe", 0.92}, // Fe3+
        {"Co", 0.79}, // Co2+
        {"Ni", 0.83}, // Ni2+
        {"Cu", 0.87}, // Cu2+
        {"Zn", 1.04}, // Zn2+
        {"Ga", 0.76}, // Ga3+
        {"Ge", 1.48},
        {"As", 0.83},
        {"Se", 1.90},
        {"Br", 1.85},
        {"Kr", 2.02},
        {"Rb", 1.75}, // Rb1+
        {"Sr", 1.40}, // Sr2+
        {"Y", 1.61},
        {"Zr", 0.98}, // Zr4+
        {"Nb", 1.33},
        {"Mo", 1.75},
        {"Tc", 2.00},
        {"Ru", 1.20},
        {"Rh", 1.22},
        {"Pd", 1.63},
        {"Ag", 1.42}, // Ag1+
        {"Cd", 1.24}, // Cd2+
        {"In", 1.93},
        {"Sn", 2.17},
        {"Sb", 0.9}, // Sb3+
        {"Te", 1.26},
        {"I", 2.06}, // I-
        {"Xe", 2.16},
        {"Cs", 1.88}, // Cs1+
        {"Ba", 1.56}, // Ba2+
        {"La", 1.83},
        {"Ce", 1.86},
        {"Pr", 1.62},
        {"Nd", 1.79},
        {"Pm", 1.76},
        {"Sm", 1.74},
        {"Eu", 1.96},
        {"Gd", 1.69},
        {"Tb", 1.66},
        {"Dy", 1.63},
        {"Ho", 1.61},
        {"Er", 1.59},
        {"Tm", 1.57},
        {"Yb", 1.54},
        {"Lu", 1.53},
        {"Hf", 1.40},
        {"Ta", 1.22},
        {"W", 1.26},
        {"Re", 1.30},
        {"Os", 1.58},
        {"Ir", 1.22},
        {"Pt", 1.72},
        {"Au", 1.66},
        {"Hg", 1.55},
        {"Tl", 1.96},
        {"Pb", 2.02},
        {"Bi", 1.73},
        {"Po", 1.21},
        {"At", 1.12},
        {"Rn", 2.30},
        {"Fr", 3.24},
        {"Ra", 2.57},
        {"Ac", 2.12},
        {"Th", 1.84},
        {"Pa", 1.60},
        {"U", 1.75},
        {"Am", 1.66},
        {"Cm", 1.65},
        {"Bk", 1.64},
        {"Cf", 1.63},
    };
} // namespace clipper_cx::data
} // namespace clipper_cx
