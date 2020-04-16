# @Author: Tristan Croll <tic20>
# @Date:   16-Apr-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 16-Apr-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def register_isolde_add(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        AtomSpecArg, ModelArg,
        FloatArg, IntArg, BoolArg, StringArg, NoArg,
        ListOf, EnumOf, RepeatOf
    )
    from chimerax.atomic import AtomsArg, ResiduesArg, StructureArg

    def register_isolde_add_water():
        from .place_ligand import place_water
        desc = CmdDesc(
            optional = [
                ('model', StructureArg),
            ],
            keyword = [
                ('position', ListOf(FloatArg,3,3)),
                ('bfactor', FloatArg),
                ('chain', StringArg),
                ('distance_cutoff', FloatArg),
                ('sim_settle', BoolArg)
            ],
            synopsis = 'Add a water molecule'
        )
        register('isolde add water', desc, place_water, logger=logger)
    register_isolde_add_water()

    def register_isolde_add_ligand():
        from .place_ligand import place_ligand
        desc = CmdDesc(
            required = [
                ('ligand_id', StringArg),
            ],
            optional = [
                ('model', StructureArg),
            ],
            keyword = [
                ('position', ListOf(FloatArg,3,3)),
                ('bfactor', FloatArg),
                ('chain', StringArg),
                ('distance_cutoff', FloatArg),
                ('sim_settle', BoolArg)
            ],
            synopsis = 'Add a ligand from the Chemical Components Dictionary'
        )
        register('isolde add ligand', desc, place_ligand, logger=logger)
    register_isolde_add_ligand()
