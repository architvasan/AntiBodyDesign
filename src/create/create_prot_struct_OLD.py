import os
from Bio.PDB import PDBList
import subprocess
import propka
from threading import Timer

def obtain_nions(logfile, ionconc):
    with open(logfile, 'r') as f:
        for line in f.readlines():
            if 'residues' in line:
                n_wat = line.split()[1]
    return int(0.0187* 0.15 * int(n_wat))#int(N_A*ionconc*vol)

def protonate_lig(lig):
   os.system(f'obabel -ipdb {lig}.pdb -opdb -O {lig}.H.pdb -h') 
 

class Watchdog(BaseException):
    def __init__(self, timeout, userHandler=None):  # timeout in seconds
        self.timeout = timeout
        self.handler = userHandler if userHandler is not None else self.defaultHandler
        self.timer = Timer(self.timeout, self.handler)

    def reset(self):
        self.timer.cancel()
        self.timer = Timer(self.timeout, self.handler)

    def stop(self):
        self.timer.cancel()

    def defaultHandler(self):
        raise self

class CreationSteps():
    def __init__(self, pdb_id, ligid, ambertools_env, leapin, leapin_ion):
        self.pdb_id = pdb_id
        self.ligid = ligid
        self.ambertools_env = ambertools_env
        self.leapin = leapin
        self.leapin_ion = leapin_ion
        os.makedirs('struct', exist_ok=True)

    def step0(self):
        # Create an instance of the PDBList class
        pdb_list = PDBList()
        
        # Download the MMCIF file using the retrieve_pdb_file method
        pdb_filename = pdb_list.retrieve_pdb_file(self.pdb_id, pdir=f"struct/{self.pdb_id}", file_format="pdb")
        
        subprocess.run(["mv", pdb_filename, f"{pdb_filename}.pdb"])

        os.system(f'grep -v -e "{self.ligid}" -e "GWS" -e "B" -e "DMS" -e "CONECT" -e "HOH" {pdb_filename}.pdb > struct/protein.pdb')
        os.system(f'grep "HOH" {pdb_filename}.pdb > struct/waters.pdb')
        os.system(f'grep "{self.ligid}" {pdb_filename}.pdb > struct/ligid.pdb')
        os.system('sed -i 1,4d struct/ligid.pdb')
        return
    
    #using pdbfixer code
    def step1(self, inp_pdb, out_pdb):
        from openmm import unit
        import tempfile
        pH = 7.0
    
        outfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
        output_pdb_filename = outfile.name
    
        timeout_seconds = 30
        watchdog = Watchdog(timeout_seconds)
        build_successful = False
        # Keep track of list of failures.
        failures = list()
        try:        
            from pdbfixer.pdbfixer import PDBFixer
            from openmm import app
            stage = "Creating PDBFixer..."
            fixer = PDBFixer(f'struct/{inp_pdb}')
            stage = "Finding missing residues..."
            fixer.findMissingResidues()
            chains = list(fixer.topology.chains())
            keys = fixer.missingResidues.keys()
            for key in list(keys):
                chain = chains[key[0]]
                print(chain)
                if key[1] == 0 or key[1] == len(list(chain.residues())):
                    del fixer.missingResidues[key]
            fixer.findNonstandardResidues()
            #fixer.missingResidues={}
            stage = "Finding nonstandard residues..."
            fixer.findNonstandardResidues()
            stage = "Replacing nonstandard residues..."
            fixer.replaceNonstandardResidues()
            stage = "Removing heterogens..."
            fixer.removeHeterogens(False)
            stage = "Finding missing atoms..."
            fixer.findMissingAtoms()
            stage = "Adding missing atoms..."
            fixer.addMissingAtoms()
            #stage = "Adding missing hydrogens..."
            #fixer.addMissingHydrogens(pH)
            stage = "Writing PDB file..."
            app.PDBFile.writeFile(fixer.topology, fixer.positions, f'struct/{out_pdb}')
            stage = "Done."
            outfile.close()
            build_successful = True
    
        except Watchdog:
            message = "timed out in stage %s" % stage
            print(message)
            failures.append((pdbcode, Exception(message)))
    
        except Exception as e:
            print("EXCEPTION DURING BUILD")
            #import traceback
            #print traceback.print_exc()
            print(str(e))
            failures.append((pdbcode, e))
        
        watchdog.stop()
        del watchdog
        return
    
    def step2(self, inp_pdb2, out_pdb2):
        print(self.ambertools_env)
        os.system(f'conda run -n {self.ambertools_env} pdb4amber -i struct/{inp_pdb2} -o struct/{out_pdb2}')
        
    def step3(self):
        protonate_lig('struct/ligid')
        os.system('grep -v -e "CONECT" struct/ligid.H.pdb > struct/ligid.H.parsed.pdb')
        #try:
        #os.system(f'conda run -n {self.ambertools_env} antechamber -i -at sybyl struct/ligid.H.pdb -fi pdb -o struct/ligid.mol2 -fo mol2 -c bcc -nc 0 -rn GWS -at gaff2')-at gaff2
        if True:
            os.system(f'conda run -n {self.ambertools_env} antechamber -j 5 -at sybyl -dr no -i struct/ligid.H.parsed.pdb -fi pdb -o struct/ligid.mol2 -fo mol2 -c bcc -nc 0 -rn LIG -at gaff2')

        os.system(f'conda run -n {self.ambertools_env} parmchk2 -i struct/ligid.mol2 -f mol2 -o struct/ligid.frcmod -s gaff2')
    
    def step4(self):
        
        with open(self.leapin, 'w') as leapfile:
            leapfile.write('')
        with open(self.leapin, 'a') as leapfile:
            leapfile.write('# Load force field parameters\n \
                                source leaprc.protein.ff19SB\n \
                                source leaprc.water.opc\n\
                                source leaprc.gaff2\n')
    
            leapfile.write('# Load protein, ligand and water molecules\n \
                                loadamberparams struct/ligid.frcmod\n \
                                ligand = loadmol2 struct/ligid.mol2\n\
                                protein = loadPDB struct/protein4amber.pdb\n \
                                waters = loadPDB struct/waters.pdb\n')
             
            leapfile.write('# Build system\n \
                                system = combine {ligand protein waters}\n \
                                savepdb system struct/system.dry.pdb\n\
                                check system\n')
             
            leapfile.write('# Solvate\n \
                                solvateBox system OPCBOX 10\n \
                                # Neutralise\n \
                                #addions2 system Cl- 0 \n \
                                #addions2 system Na+ 0 \n\
                                quit')
        os.system(f'conda run -n ambertools tleap -f {self.leapin}')
        os.system(f'rm {self.leapin}')

    def step5(self):
        n_ions = (obtain_nions('leap.log', 0.15))
        print(n_ions)
        os.system('grep -vwE "CONECT" struct/protein4amber.pdb > struct/protein4amber1.pdb')
        os.system('mv struct/protein4amber1.pdb struct/protein4amber.pdb')
        with open(self.leapin_ion, 'w') as leapfile:
            leapfile.write('')
        with open(self.leapin_ion, 'a') as leapfile:
            leapfile.write('# Load protein, ligand and water molecules\n \
                    source leaprc.water.opc\n\
                    source leaprc.protein.ff19SB\n\
                    source leaprc.gaff2\n\
                    loadamberparams struct/ligid.frcmod\n \
                    ligand1 = loadmol2 struct/ligid.mol2\n\
                    protein1 = loadPDB struct/protein4amber.pdb\n \
                    waters1 = loadPDB struct/waters.pdb\n \
                    system1 = combine {ligand1 protein1 waters1} \n \
                    check system1\n')
            leapfile.write(f'# Solvate\n \
                    solvateBox system1 OPCBOX 10\n \
                    addions2 system1 Cl- 0 \n \
                    addions2 system1 Na+ 0 \n \
                    addions2 system1 Cl- {n_ions} \n \
                    addions2 system1 Na+ {n_ions} \n \
                    savePDB system1 struct/system1.pdb\n\
                    saveAmberParm system1 struct/prmtop1 struct/inpcrd1\n\
                    quit')
        os.system(f'conda run -n ambertools tleap -f {self.leapin_ion}')
        #os.system(f'rm {self.leapin_ion}')
        
        #system = loadPDB struct/system.dry.pdb\n \

def run_steps():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p',
                        '--pdbid',
                        type=str,
                        help='PDB ID')

    parser.add_argument('-l',
                        '--ligid',
                        type=str,
                        help='Ligand ID')
    
    parser.add_argument('-a',
                        '--ambertools',
                        type=str,
                        help='AmberTools Conda')

    parser.add_argument('-i',
                        '--leapin',
                        type=str,
                        help='leapin file')

    parser.add_argument('-I',
                        '--leapion',
                        type=str,
                        help='leapin_ion file')
    
    args = parser.parse_args()

    Creation = CreationSteps(args.pdbid, args.ligid, args.ambertools, args.leapin, args.leapion)
    '''
    Obaining PDB
    '''
    print("Obaining PDB")
    Creation.step0()

    #'''
    #Fixing loops, protonating, adding hydrogens, etc.
    #'''
    #print("Amberizing")
    #Creation.step2('protein.pdb', 'protein_amber1.pdb')

    '''
    Making protein pdb into amber format
    '''
    print("Amberizing protein pdb")
    Creation.step1('protein.pdb', 'protein_fixed.pdb')

    '''
    Fixing loops, protonating, adding hydrogens, etc.
    '''
    print("Fixing PDB")
    Creation.step2('protein_fixed.pdb', 'protein4amber.pdb')

    '''
    Parametrizing ligand
    '''
    print("Parameterizing ligand")
    Creation.step3()

    '''
    Combining protein, ligand, crystal waters
    '''
    print("Combining everything")
    #Creation.step4()

    '''
    Solvating and ionizing
    '''
    print("solvating/ionizing")
    Creation.step5()

    try:
        os.system('rm ANTECHAMBER_AC.AC  ANTECHAMBER_AC.AC0  ANTECHAMBER_BOND_TYPE.AC  ANTECHAMBER_BOND_TYPE.AC0  ATOMTYPE.INF sqm.in sqm.out')
    except:
        pass
    return


if __name__ == "__main__":
    run_steps()










if False:
    pdb_id = "6yhr"
    ligid = "ADP"
    ambertools_env = 'ambertools'
    
    step0_pdb(pdb_id, ligid)
    step0p5_fix_pdb()
    step1_amber_create_protfile()

    subprocess.run(["propka3", "-d", "--protonate-all", f"struct/protein.pdb"])
    
        
    #os.system('conda run -n ambertools source leaprc.protein.ff14SB')
    
    def step1_param():
        return
    os.system('rm leap.log')
    
    with open('leapin', 'w') as leapfile:
        leapfile.write('')
    with open('leapin', 'a') as leapfile:
        leapfile.write('# Load force field parameters\n \
                            source leaprc.protein.ff19SB\n \
                            source leaprc.water.opc\n')
    
        leapfile.write('# Load protein, ligand and water molecules\n \
                            protein = loadPDB struct/protein4amber.pdb\n \
                            waters = loadPDB struct/waters.pdb\n')
    
        leapfile.write('# Build system\n \
                            system = combine {protein waters}\n \
                            savepdb system struct/system.dry.pdb\n\
                            check system\n')
        
        leapfile.write('# Solvate\n \
                            solvateBox system OPCBOX 10\n \
                            # Neutralise\n \
                            addions2 system Cl- 0 \n \
                            addions2 system Na+ 0 \n')
        
        leapfile.write('# Save AMBER input files\n\
                            #savePDB system struct/system_noion.pdb\n\
                            #saveAmberParm system prmtop inpcrd\n\
                            quit')
    
    os.system('conda run -n ambertools tleap -f leapin')
    
    n_ions = (obtain_nions('leap.log', 0.15))
    
    with open('leapin_ion', 'w') as leapfile:
        leapfile.write('')
    with open('leapin_ion', 'a') as leapfile:
        leapfile.write(f'# Load protein, ligand and water molecules\n \
                source leaprc.water.opc\n\
                source leaprc.protein.ff19SB\n\
                system = loadPDB struct/system.dry.pdb\n \
                solvateBox system OPCBOX 10\n \
                addions2 system Cl- {n_ions} \n \
                addions2 system Na+ {n_ions} \n \
                savePDB system struct/system.pdb\n\
                saveAmberParm system struct/prmtop struct/inpcrd\n\
                #struct/system.parm7 struct/system.rst7\n\
                quit')
    
    os.system('conda run -n ambertools tleap -f leapin_ion')
    
    
    
    import argparse
    
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')
    
    args = parser.parse_args()
    
