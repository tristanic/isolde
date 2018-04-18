
.. contents::
    :local:

Getting Started
===============

Assuming you already have ISOLDE installed (if not, you can do so via
`Tools/More Tools...` in the ChimeraX menu), then you can start it up via
`Tools/General/ISOLDE`. This should yield a new floating panel looking something
like this:

.. image:: images/isolde_initial_screen.png
    :alt: Initial GUI appearance

If you already have a model loaded, you'll probably notice that it's gained a
few (but hopefully not *too* many) new features:

.. figure:: images/basic_display.jpg
    :alt: Basic rota/rama/omega markup

    Basics of ISOLDE validation markup (*sans* map)

    +---+--------------------------------------------------------------+
    | a | Crosshairs denoting the pivot point of the display. Red,     |
    |   | green and blue point along the x, y and z axes respectively. |
    +---+--------------------------------------------------------------+
    | b | This exclamation mark/spiral motif denotes a rotamer outlier |
    |   | (that is, a sidechain in an unlikely conformation). The more |
    |   | unlikely the conformation, the larger and redder the         |
    |   | indicator becomes. Below and to the left you can see a less  |
    |   | severe "iffy" rotamer.                                       |
    +---+--------------------------------------------------------------+
    | c | The red trapezoids you see here are highlighting non-proline |
    |   | *cis* peptide bonds (where the amide hydrogen and carbonyl   |
    |   | oxygen are pointing in the same direction). In the real world|
    |   | these are vanishingly rare (around 3 per 10,000 amino acid   |
    |   | residues), and real ones tend to be heavily stabilised by    |
    |   | surrounding packing/H-bond interactions (and hence are       |
    |   | usually among the better-resolved sites in the molecule). A  |
    |   | string of non-proline *cis* bonds on a flexible loop as seen |
    |   | here is essentially impossible. When all atoms are shown, the|
    |   | trapezoids can be seen to fill in the "cup" formed by the    |
    |   | C, O, N and CA atoms. The less-rare proline *cis* bonds are  |
    |   | similarly shown in green, and peptide bonds twisted more than|
    |   | 30 degrees from planar in yellow.                            |
    +---+--------------------------------------------------------------+
    |d,e| The protein backbone is not infinitely free to move, but has |
    |   | clearly preferred conformations that have been well          |
    |   | characterised by studying high-resolution structures. The    |
    |   | best-established method for characterising backbone          |
    |   | conformation is via the Ramachandran plot, a plot of the phi |
    |   | (C-N-CA-C) and psi (N-CA-C-N) dihedral angles against each   |
    |   | other. The probabilities of finding different (phi, psi)     |
    |   | combinations have been mapped out in high detail for         |
    |   | various groups of amino acids [MolProbity]_. While ISOLDE   |
    |   | also provides a Ramachandran plot, the current probability   |
    |   | score for each protein residue is mapped in real time to the |
    |   | colour of its alpha carbon (CA) atom as shown. Green denotes |
    |   | a "happy" residue, yellow is marginal (possible, but somewhat|
    |   | strained), and red is a serious outlier (still possible, but |
    |   | needs very strong support to justify it).                    |
    +---+--------------------------------------------------------------+

.. [MolProbity] https://doi.org/10.1107/S0907444909042073
