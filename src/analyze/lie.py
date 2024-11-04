from openmm import *
from openmm.app import *
from openmm.unit import *
from pdbfixer import PDBFixer

def pdbfixit(inpdb, outpdb):
    fixer = PDBFixer(filename=inpdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()
    
    # Save the fixed structure to a temporary file
    fixed_pdb_file = outpdb
    with open(fixed_pdb_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

def coulomb_energy(solute_scale, solvent_scale):
    context.setParameter("solute_scale", solute_scale)
    context.setParameter("solvent_scale", solvent_scale)
    return context.getState(getEnergy=True,
                            groups={0}).getPotentialEnergy()


def energy(context, solute_coulomb_scale, solute_lj_scale, solvent_coulomb_scale, solvent_lj_scale):
    context.setParameter("solute_coulomb_scale", solute_coulomb_scale)
    context.setParameter("solute_lj_scale", solute_lj_scale)
    context.setParameter("solvent_coulomb_scale", solvent_coulomb_scale)
    context.setParameter("solvent_lj_scale", solvent_lj_scale)
    return context.getState(getEnergy=True, groups={0}).getPotentialEnergy()

def lie(inppdb,
        fixedpdb,
        sel_chain,
        sel_rst,
        sel_rend):

    #pdb = PDBFile('fixed_structure.pdb')
    pdbfixit(inppdb, fixedpdb)
    pdb = PDBFile(fixedpdb)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic)

    solvent = set([a.index for a in pdb.topology.atoms() if a.residue.chain.id == sel_chain and int(a.residue.id)>sel_rst and int(a.residue.id)<sel_rend])
    protein = set([a.index for a in pdb.topology.atoms() if a.index not in solvent])
        
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            force.setForceGroup(0)
            force.addGlobalParameter("solute_coulomb_scale", 1)
            force.addGlobalParameter("solute_lj_scale", 1)
            force.addGlobalParameter("solvent_coulomb_scale", 1)
            force.addGlobalParameter("solvent_lj_scale", 1)
            for i in range(force.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(i)
                force.setParticleParameters(i, 0, 0, 0)
                if i in protein:
                    force.addParticleParameterOffset("solute_coulomb_scale", i, charge, 0, 0)
                    force.addParticleParameterOffset("solute_lj_scale", i, 0, sigma, epsilon)
                else:
                    force.addParticleParameterOffset("solvent_coulomb_scale", i, charge, 0, 0)
                    force.addParticleParameterOffset("solvent_lj_scale", i, 0, sigma, epsilon)
            for i in range(force.getNumExceptions()):
                p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
                force.setExceptionParameters(i, p1, p2, 0, 0, 0)
        else:
            force.setForceGroup(2)
    
    integrator = VerletIntegrator(0.001*picosecond)
    context = Context(system, integrator)
    context.setPositions(pdb.positions)
    
    
    total_coulomb = energy(context, 1, 0, 1, 0)
    solute_coulomb = energy(context, 1, 0, 0, 0)
    solvent_coulomb = energy(context, 0, 0, 1, 0)
    total_lj = energy(context, 0, 1, 0, 1)
    solute_lj = energy(context, 0, 1, 0, 0)
    solvent_lj = energy(context, 0, 0, 0, 1)
    coul_final = total_coulomb - solute_coulomb - solvent_coulomb
    lj_final = total_lj - solute_lj - solvent_lj
    return coul_final.value_in_unit(kilocalories_per_mole), lj_final.value_in_unit(kilocalories_per_mole)


if __name__ == "__main__": 
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-p',
                        '--pdbinp',
                        type=str,
                        help='input pdb')

    parser.add_argument('-pf',
                        '--pdbfixed',
                        type=str,
                        help='fixed pdb out name')

    parser.add_argument('-c',
                        '--chainsel',
                        type=str,
                        help='chain id for sel (A, B)')

    parser.add_argument('-rs',
                        '--rstrtsel',
                        type=int,
                        help='starting resid for sel (int)')

    parser.add_argument('-re',
                        '--rendsel',
                        type=int,
                        help='ending resid for sel (int)')


    args = parser.parse_args()

    coul, lj = lie(args.pdbinp,
                    args.pdbfixed,
                    args.chainsel,
                    args.rstrtsel,
                    args.rendsel)


    print('Coulomb interaction energy:', coul)
    print('LJ interaction energy:', lj)
