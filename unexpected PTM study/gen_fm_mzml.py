import numpy as np
from pyopenms import *
import loadSptxtFile
from random import shuffle
import sys

def task(mzmlinput, splibtxtinput, fakemod, mod_frac):
    print('Reading mzML file...')
    mzml_org = MSExperiment()
    MzMLFile().load(mzmlinput, mzml_org)
    print('Reading spectral library...')
    speclib = loadSptxtFile.SpTXT()
    speclib.parse(splibtxtinput)

    # create index dictionary where key is mzML index, value is speclib index
    # keep in mind it will work only on original mzML file, if noise added then shift in index due to omission of some spectra
    # needs to be considered
    splib_index = {}
    fm_index = []
    fake_mod = fakemod
    for no in range(0, speclib.counts()):
        mzml_id = speclib.getLibItemOnIDX(no).ScanID
        splib_index[mzml_id] = no
        if fake_mod in speclib.getLibItemOnIDX(no).seq:
            fm_index.append(no)

    # create index dictionary relating position of mzML spectra to its scan number
    mzml_dict_index = {}
    for i in range(0, len(mzml_org.getSpectra())):
        current_index = int(str(mzml_org[i].getNativeID()).split('=')[-1][:-1])
        mzml_dict_index[i] = current_index

    # randomly select proportion of the spectra to be added as a replacement to original mzML
    mod_fraction = float(mod_frac)
    np.random.seed()
    shuffle(fm_index)
    sel_mod = fm_index[:int(mod_fraction * len(fm_index))]
    new_mzml = MSExperiment()

    for i in range(len(mzml_org.getSpectra())):
        mzml_index = mzml_dict_index[i]
        mzml_spec = mzml_org[i]
        if mzml_index in splib_index.keys() and splib_index[mzml_index] in sel_mod:
            splib_ind = splib_index[mzml_index]
            print(f'splib index is: {splib_ind}')
            splib_peaks = speclib.getLibItemOnIDX(splib_ind)
            splib_mz = splib_peaks.mz
            splib_int = splib_peaks.intensity
            mzml_spec.set_peaks([splib_mz, splib_int])
            new_mzml.addSpectrum(mzml_spec)
        else:
            new_mzml.addSpectrum(mzml_spec)

    inputname = mzmlinput.split("/")[-1]
    MzMLFile().store('semi_' + inputname, new_mzml)
    print('Done!')

def usage():
    print('Usage: \n --> to generate mzML with fraction of spectra with fake modification: \n \
    first param: input mzML file \n \
    second param: spectral library (sptxt) \n \
    third param: fake modification considered (e.g. "V[123]") \n \
    fourth param: percentage of all fakely modified spectra to be used as replacement in the original mzML (0-1)')

if __name__ == '__main__':
    if len(sys.argv) != 5:
        usage()

    if len(sys.argv) == 5:
        task(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
