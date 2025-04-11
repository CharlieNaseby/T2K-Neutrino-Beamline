#source /cvmfs/beam-physics.cern.ch/bdsim/x86_64-centos7-gcc9-opt/bdsim-env-develop-g4v10.7.2.3-boost.sh


echo "input number of particles to simulate"
read NPROTONS

rm output.root
bdsim --file=survey/optimised.gmad --outfile=output --batch --ngenerate=$NPROTONS --seed=1989

rebdsimOptics output.root ./optics_dump/optics_output.root --emittanceOnTheFly

cd optics_dump/

./Optics optics_output.root optics_output.csv

cd ../
