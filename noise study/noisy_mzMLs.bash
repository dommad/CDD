ilim=0.6 #intensity parameter
plim=0.7 #number of peaks parameter

for i in `ls *.mzML`;
 do python3 gen_noise.py $i $ilim $plim; done