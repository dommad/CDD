for i in `ls interact*pep.xml | grep -v ipro`;
 do outputname=$(echo $i | cut -d'.' -f "1" | cut -d'-' -f "2,3");
 /tools/tpp51/bin/spectrast -cNpre_${outputname} -cP0 -c_NAA0 $i;
 /tools/tpp51/bin/spectrast -cNbr_${outputname} -cAB -c_XPK1000000 -c_NAA0 pre_${outputname}.splib
 /tools/tpp51/bin/spectrast -Mspectrast.usermods -cNfake_mod_${outputname} -c_XPK1000000 -c_NAA0 -cAM -cx'V[fake]' br_${outputname}.splib; done
