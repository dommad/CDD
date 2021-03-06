fakemod="V[123]"
fakefrac=0.7

for i in `ls *mzML`;
 do corename=$(echo $i | cut -d'.' -f 1);
 python3 gen_fm_mzml.py $i fake_mod_${corename}.sptxt $fakemod $fakefrac;
 /tools/bin/comet semi_${i};
 /tools/tpp51/bin/xinteract -dDECOY -ip -p0 -OPd -N${corename} semi_${corename}.pep.xml;
 /tools/tpp51/cgi-bin/PepXMLViewer.cgi -B exportSpreadsheet 1 -C "Iiprobability,Pprobability,Gspectrum,Gstart_scan,Sexpect,Gions2,Gpeptide,Gprotein,Gcalc_neutral_pep_mass,Sxcorr,Sspscore,Ssprank,Sdeltacn,Sdeltacnstar,Pfval,Gmassdiff" -I "interact-${cpname}.pep.xml"; done
