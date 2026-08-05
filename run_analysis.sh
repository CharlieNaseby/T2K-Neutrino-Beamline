#source /cvmfs/beam-physics.cern.ch/bdsim/x86_64-centos7-gcc9-opt/bdsim-env-develop-g4v10.7.2.3-boost.sh

echo "Make sure make_beamline.py is set to the plot configuration!"
echo "And that the fit result file has the runnumber included e.g. ./bdsim_optimiser/fit_results_910216.root"

echo "input run number"
read RUNNUMBER

echo "input number of particles to simulate MUST use same as bdsim optimiser if you want exact matching (usually 100)"
read NPROTONS


cd survey
python make_beamline.py ../bdsim_optimiser/fit_results_${RUNNUMBER}.root
cp optimised.gmad optimised_${RUNNUMBER}.gmad
cd ..
rm output.root
bdsim --file=survey/optimised_${RUNNUMBER}.gmad --outfile=output --batch --ngenerate=$NPROTONS --seed=1989

rebdsimOptics output.root ./optics_dump/optics_output.root --emittanceOnTheFly


cd bdsim_optimiser
root -l -b -q "plot_beam_result.C(\"fit_results_${RUNNUMBER}\", \"../optics_dump/optics_output.root\", \"\")"
root -l -b -q "plot_fit_result.C(\"fit_results_${RUNNUMBER}\")"
python make_pretty_beam_plots.py ${RUNNUMBER}
