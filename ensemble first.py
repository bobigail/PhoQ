from prody import *
from pylab import *
#ion()

# select PDBs to work with
pdbids = ['1yax', '3bqa', '3bq8', '3cgz', '3cgy', '4uey']
pdbfiles = fetchPDB(*pdbids, compressed=False)

# print len(pdbfiles)

# identify reference chain
ref_structure = parsePDB('4uey')
ref_selection = ref_structure.select('calpha')
ref_chain = ref_selection.getHierView().getChain('A')

# create ensemble
startLogfile('PhoQ_pca')
ensemble = PDBEnsemble('PhoQ')

# define atoms and coordinates for the ensemble reference
ensemble.setAtoms(ref_chain)
ensemble.setCoords(ref_chain)

for pdbid in pdbids:
    
    hv = parsePDB(pdbid).getHierView()
    for chain in hv:
        #structure = parsePDB(pdbid, subset='calpha')
        mappings = mapOntoChain(chain, ref_chain)
        if mappings != []:
            atommap = mappings[0][0]
            ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))
        
repr(ensemble)

## in the tutorial it says len(ensemble) == len(pdbfiles), but in this case, len(ensemble should be equal to the sum of the number of chains?

ensemble.iterpose()

closeLogfile('PhoQ_pca')


writePDB('PhoQ_ensemble.pdb', ensemble)



