

for ID in `seq 1 1`
do
cp inputDefault inputEventID${ID}
toReplace='pythia_file'
#replaceBy="pythia_file\t=\tinitialData\/EventID"$ID"\/rejectionSampling\/SmearedDistributionT0.2005fm"
replaceBy="pythia_file\t=\tLowInput\/SmearedDistributionT0.2005fm"
# = initialData/EventID"${ID}
#'/rejectionSampling/SmearedDistributionT0.2005fm'
sed -i -e 's/'$toReplace'/'$replaceBy'/g' inputEventID${ID}
toReplace='jobname'
replaceBy="jobname\t=\tClustertestDCA_F"
# = initialData/EventID"${ID}
#'/rejectionSampling/SmearedDistributionT0.2005fm'
sed -i -e 's/'$toReplace'/'$replaceBy'/g' inputEventID${ID}
#./cascade --output.directory=/home/greif/HERWIGNEW/PLOTS/VIDEOS/RESULT/ inputEventID${ID}
./cascade --output.directory=output inputEventID${ID}
done

