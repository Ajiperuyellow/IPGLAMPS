for ID in `seq 1 10`
do
cp inputDefault inputEventID${ID}
toReplace='pythia_file'
replaceBy="pythia_file\t=\tinitialData\/EventID"$ID"\/rejectionSampling\/SmearedDistributionT0.2005fm"
# = initialData/EventID"${ID}
#'/rejectionSampling/SmearedDistributionT0.2005fm'
sed -i -e 's/'$toReplace'/'$replaceBy'/g' inputEventID${ID}
toReplace='jobname'
replaceBy="jobname\t=\tconst10mbEventID"$ID
# = initialData/EventID"${ID}
#'/rejectionSampling/SmearedDistributionT0.2005fm'
sed -i -e 's/'$toReplace'/'$replaceBy'/g' inputEventID${ID}

cp jobscriptDefault jobscriptEventID${ID}
sed -i -e 's/defaultinput/'inputEventID${ID}'/g' jobscriptEventID${ID} 

sbatch jobscriptEventID${ID}
#./cascade inputEventID${ID}
done

