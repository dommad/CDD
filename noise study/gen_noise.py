import pandas as pd
import numpy as np
from pyopenms import *
import random
import sys

def noise(original_mzml, perc_int, perc_peaks):

    mzml_file = original_mzml
    mzml_dict = MSExperiment()
    MzMLFile().load(mzml_file, mzml_dict)

    new_mzml = MSExperiment()

    for i in range(len(mzml_dict.getSpectra())):
        print(f'iteration: {i}')
        np.random.seed()
        per_int = float(perc_int)
        per_rand_peaks = float(perc_peaks)
        current_spectrum = mzml_dict[i]
        peak_list = current_spectrum.get_peaks()
        mz_ar = peak_list[0]
        int_ar = peak_list[1]
        if len(mz_ar) == 0:
            continue
        no_peaks = len(mz_ar)
        mz_min = min(mz_ar)
        max_sorted = sorted(list(mz_ar), reverse=False)

        if len(max_sorted) >= 2:
            sec_max = max_sorted[-2]
        else:
            sec_max = max_sorted[-1]

        rand_ints = []
        rand_mzs = []
        last_peaks = np.array(max_sorted[-100:])
        sd = last_peaks.std()

        for j in range(int(no_peaks * per_rand_peaks)):
            new_int = abs(random.gauss(per_int, sd))
            new_mz = random.uniform(mz_min, sec_max)
            rand_ints.append(new_int)
            rand_mzs.append(new_mz)

        mixed_mz = list(mz_ar) + rand_mzs
        mixed_int = list(int_ar) + rand_ints

        df = pd.DataFrame()
        df['mz'] = mixed_mz
        df['int'] = mixed_int
        df = df.sort_values(['mz'], ascending=True)
        final_mz = list(df['mz'])
        final_int = list(df['int'])

        current_spectrum.set_peaks([final_mz, final_int])
        new_mzml.addSpectrum(current_spectrum)

    MzMLFile().store('noise_' + original_mzml, new_mzml)

def usage():
    print('Usage: \n --> to generate mzML with noise: \n \
            first param: original mzML \n \
            second param: percentage of intensity of 2nd largest peaks used as mean (0-1) \n \
            third param: number of noise peaks as percentage of all peaks (0-1)')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        usage()
    if len(sys.argv) == 4:
        noise(sys.argv[1], sys.argv[2],sys.argv[3])


