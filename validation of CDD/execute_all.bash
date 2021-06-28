#example of file to conduct validation of CDD method

FDRs=(0.001 0.005 0.01 0.015 0.02 0.025 0.03);
for fdr in ${FDRs[@]};
do echo $fdr;
files=(2 3 4 5 6 8);

for i in ${files[@]};
do echo $i;

out=$(find . -name "H*pH${i}*.mzXML");
core=$(echo $out | cut -d '.' -f 2 | cut -d '/' -f 2);
echo $core;

#prepare files for generation of ground truth
/tools/bin/comet.2018014.linux.exe -Pcomet_td $i;
java -Xmx50g -jar /tools/msfragger2019/MSFragger-20190222.jar fragger.params $i;
/tools/bin/xinteract -dDECOY -p0 -OAPd -PPM -N${core}_comet ${core}.pep.xml;
/tools/bin/xinteract -dDECOY -p0 -OAPd -PPM -N${core}_msfrag ${core}.pepXML;

#if synthetic validation data set used, uncomment the line below and hide the next one
#python3 refine_comet_msfrag_synthetic.py interact-${core}_comet.pep.xml interact-${core}_msfrag.pep.xml ${core} pools_pep.csv $id;
python3 refine_comet_msfrag_natural.py interact-${core}_comet.pep.xml interact-${core}_msfrag.pep.xml $core; 
python3 val_files_pmz.py ${core}.pkl $i; 


/tools/bin/comet.2018014.linux.exe -Pcomet_t_human -N${core}_t refined_${i};
/tools/bin/comet.2018014.linux.exe -Pcomet_td -N${core}_td refined_$i;

option1="MNkd";
option2="MNEA";

/tools/bin/xinteract -dDECOY -p0 -O${option1} -N${core}_td_O${option1} ${core}_td.pep.xml;
/tools/bin/xinteract  -dDECOY -p0 -O${option2} -N${core}_td_O${option2} ${core}_td.pep.xml;

#PeptideProphet with CDD as negative model
/tools/tpp51_eval_fixedGumbel/bin/xinteract -p0 -O${option1} -N${core}_t_O${option1} ${core}_t.pep.xml; 
/tools/tpp51_eval_fixedGumbel/bin/xinteract -p0 -O${option2} -N${core}_t_O${option2} ${core}_t.pep.xml;

python3 decoy_original.py ${core}.pkl interact-${core}_t.pep.xml interact-${core}_t_OMNEk.pep.xml $fdr output.csv 1 cdd_params_Elite.yaml; done
done;

