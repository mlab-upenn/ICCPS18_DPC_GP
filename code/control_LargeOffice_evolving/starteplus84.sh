#!/bin/bash
# Arguments: <directory of IDF file> <IDF file name without idf> <weather file name without extension>

cd $1

# Remove old files
rm -f  eplusout.end
rm -f  eplusout.eso
rm -f  eplusout.rdd
rm -f  eplusout.edd
rm -f  eplusout.mdd
rm -f  eplusout.dbg
rm -f  eplusout.eio
rm -f  eplusout.err
rm -f  eplusout.dxf
rm -f  eplusout.csv
rm -f  eplusout.tab
rm -f  eplusout.txt
rm -f  eplusmtr.csv
rm -f  eplusmtr.tab
rm -f  eplusmtr.txt
rm -f  eplusout.sln
rm -f  epluszsz.csv
rm -f  epluszsz.tab
rm -f  epluszsz.txt
rm -f  eplusssz.csv
rm -f  eplusssz.tab
rm -f  eplusssz.txt
rm -f  eplusout.mtr
rm -f  eplusout.mtd
rm -f  eplusout.cif
rm -f  eplusout.bnd
rm -f  eplusout.dbg
rm -f  eplusout.sci
rm -f  eplusout.cfp
rm -f  eplusmap.csv
rm -f  eplusmap.txt
rm -f  eplusmap.tab
rm -f  eplustbl.csv
rm -f  eplustbl.txt
rm -f  eplustbl.tab
rm -f  eplustbl.htm
rm -f  eplusout.log
rm -f  eplusout.svg
rm -f  eplusout.shd
rm -f  eplusout.wrl
rm -f  eplusoutscreen.csv
rm -f  eplusout.delightin
rm -f  eplusout.delightout
rm -f  eplusout.delighteldmp
rm -f  eplusout.delightdfdmp
rm -f  eplusout.sparklog
rm -f  in.imf
rm -f  in.idf
rm -f  out.idf
rm -f  eplusout.inp
rm -f  in.epw
rm -f  eplusout.audit
rm -f  eplusmtr.inp
rm -f  expanded.idf
rm -f  expandedidf.err
rm -f  readvars.audit

rm -f  eplusout.sql
rm -f  sqlite.err

# Run Simulation
/Applications/EnergyPlus-8-4-0/energyplus -w /Applications/EnergyPlus-8-4-0/WeatherData/$3.epw $2.idf

# Generate the output results
echo eplusout.eso >eplusout.inp
echo eplusout.csv >>eplusout.inp

if [ -f Output/$2.csv ]; then
    rm Output/$2.csv
fi

/Applications/EnergyPlus-8-4-0/PostProcess/ReadVarsESO eplusout.inp unlimited
if [ -f eplusout.csv ]; then
    mkdir -p Output
    mv eplusout.csv Output/$2.csv
fi
