import numpy as np
import deepdish as dd
from pyteomics import pepxml
import sys


def evals(pepxmlfile):
    scores2 = {}
    scores3 = {}
    core = pepxmlfile.split('.')[0].split('/')[-1]
    pep_dict = pepxml.read(pepxmlfile)
    print(f'processing {core}...')
    for spectrum in pep_dict:
        
        if 'DECOY' not in el['search_hit'][0]['proteins'][0]['protein']:
            continue
        if 'search_hit' in spectrum.keys():
            e_val = spectrum['search_hit'][0]['search_score']['expect']
            spec = str(spectrum['spectrum'])
            #spec = ''.join(spec.split('.'))
            pep = spectrum['search_hit'][0]['peptide']
            charge = int(spectrum['assumed_charge'])
            if len(pep) >= 7:

                logged = -0.02 * np.log(float(e_val) / 500.)

                if charge == 2:
                    scores2[spec] = logged
                if charge == 3:
                    scores3[spec] = logged

    dd.io.save(f'{core}_ch2_evals', scores2, compression=('blosc',9))
    dd.io.save(f'{core}_ch3_evals', scores3, compression=('blosc', 9))

if __name__ == '__main__':
    if len(sys.argv) != 2:
        usage()
    if len(sys.argv) == 2:
        evals(sys.argv[1])

