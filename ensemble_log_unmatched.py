from prody import *
from pylab import *

#ion()

# select PDBs to work with
pathPDBFolder("../PDBs")
pdbids = ['1yax', '3bqa', '3bq8', '4uey']
pdbfiles = fetchPDB(*pdbids, copy=False)

# print len(pdbfiles)

# identify reference chain
ref_structure = parsePDB('4uey')
ref_selection = ref_structure.select('calpha')
ref_chain = ref_selection.getHierView().getChain('A')

#parameters for comparison
sequence_identity = 70
sequence_coverage = 70

# create ensemble
startLogfile('PhoQ_pca_4structs')
ensemble = PDBEnsemble('PhoQ')

# define atoms and coordinates for the ensemble reference
ensemble.setAtoms(ref_chain)
ensemble.setCoords(ref_chain)

NoStructure = []
unmapped = []

for pdbid in pdbids:
    
    hv = parsePDB(pdbid).getHierView()
    for chain in hv:
        structure = parsePDB(pdbid, subset='calpha')
        # here's the new bit        
        if structure is None:
            NoStructure.append(pdbid)
            continue
        
        mappings = mapOntoChain(chain, ref_chain, seqid=sequence_identity, coverage=sequence_coverage)
        if len(mappings) == 0:
            unmapped.append((pdbid,chain))
            continue
        
        atommap = mappings[0][0]
        ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))
 
print NoStructure       
print unmapped
       
        
        
repr(ensemble)

## in the tutorial it says len(ensemble) == len(pdbfiles), but in this case, len(ensemble should be equal to the sum of the number of chains?

ensemble.iterpose()

closeLogfile('PhoQ_pca_4structs')


writePDB('PhoQ_ensemble.pdb', ensemble)



