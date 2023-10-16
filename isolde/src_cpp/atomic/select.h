#include "../molc.h"

using namespace atomstruct;

std::set<Residue *> residue_neighbors(Residue *r)
{
    std::set<Residue *> neighbors;
    for (auto a: r->atoms())
        for (auto n: a-> neighbors())
            if (n->residue() != r)
                neighbors.insert(n->residue());
    return neighbors;
}

void recursive_select(Residue* r, int i, int nsteps) {
    for (auto a: r->atoms())
        a->set_selected(true);
    if (i==nsteps) return;
    for (auto n: residue_neighbors(r))
    {
        if (n->atoms()[0]->selected()) continue;
        recursive_select(n, i+1, nsteps);
    }
    
}
extern "C" EXPORT void
expand_selection(void *residues, size_t nres, int nsteps)
{
    Residue **rr = static_cast<Residue **>(residues);
    try {
        for (size_t i=0; i<nres; ++i)
        {
            recursive_select(rr[i], 0, nsteps);
        }
    } catch (...) {
        molc_error();
    }
}

