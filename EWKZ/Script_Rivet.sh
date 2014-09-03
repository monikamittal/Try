#!/bin/sh
#root -l
foreach whichfile (`cat AllLHEFiles.txt`)
echo $whichfile
set base = `basename $whichfile .lhe`
echo $base
set outputfilename=OUTAIDA_1Sept/Out_${base}.aida
sed -e "s#INPUTLHEFILE#${whichfile}#g" \
       -e "s#OUTPUTAIDAFILE#${outputfilename}#g" <VBFZ_8TeV.py> Python_RIVET.py
cmsRun Python_RIVET.py
aida2root $outputfilename 
set base1 = `basename $outputfilename .aida`
set outputrootfile=$base1.root
mv $outputrootfile OUTROOT_1Sept/
end
