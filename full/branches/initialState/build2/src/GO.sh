#!/bin/bash


cp inputfile_generic.in temp.in 

toReplace='DEFAULT'
replaceBy=$1
sed -i -e 's/'$toReplace'/'$replaceBy'/g' temp.in

toReplace='SAVERUN'
replaceBy=$2
sed -i -e 's/'$toReplace'/'$replaceBy'/g' temp.in

cp temp.in /home/greif/HERWIGNEW/

cp $3 /home/greif/HERWIGNEW/
cd /home/greif/HERWIGNEW/

source bin/activate

#rm myEPLUSEMINUS*
Herwig read temp.in
Herwig run $2.run -N 10 
#-s $2

echo "$1 HERWIG run done."

cp $replaceBy.log /home/greif/Full-PQCD-BAMPS/full/branches/initialState/build2/src/
