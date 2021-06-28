import deepdish as dd
import scipy.stats as st

paths = glob.glob(os.path.join("/data/dominik/q_decoy/function/","*ch3_evals"))

TEVs = []

for i in paths:
    CurDict = dd.io.load(i)
# remove the artifacts caused by Comet hitting 999 limit of e-value
    CurTEVs = [el for el in CurDict.values() if el > -0.0138]
    TEVs += CurTEVs

params = st.gumbel_r.fit(TEVs)

print(params)

