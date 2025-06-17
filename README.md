# T2K-Neutrino-Beamline

This is a respository intended to house all code necessary to simulate the T2K neutrino beamline (preparation section) in the accelerator tracking software BDSIM (L.J. Nevay et al., BDSIM: An Accelerator Tracking Code with Particle-Matter Interactions, Computer Physics Communications 252 107200 (2020) https://doi.org/10.1016/j.cpc.2020.107200). This intends to include simulation of beam loss and comparison to existing beam loss measurements at the beamline with the goal to determine improvements in optics or collimator design to reduce beam loss particularly in the preparation section of the neutrino beamline.

To get started first install BDSIM, alternatively and the recommended option; use the .def file to create a singularity container which will automatically download and install all dependencies. Most scripts are in Python, for this you require a python 3+ install (tested with 3.9.21 and 3.12.3) (also confirmed to work within the container). BDSIM requires GEANT4 and CLHEP, if you are doing the install yourself check the .def file for build information or the respective repositories.

# **Building the container**

You will need sudo permissions on your host to build this container however once built it can be uploaded to another host and used without root permission.

    git clone git@github.com:CharlieNaseby/T2K-Neutrino-Beamline.git

    sudo apptainer build --sandbox*containername* T2K-Neutrino-Beamline/rocky_9.def

This will take approximately 1.5hrs to build. To open the container use the open_container.sh script changing the container name to match your chosen name. You may also wish to bind a directory to be available within the container, first

    mkdir containername/bind_path/

Then add this bind argument to the open_container.sh script with the argument **-B/bind_path/**

if you wish to use this on another machine 

    chmod a+r -R containername

    tar -czf tarball.tar.gz containername

then copy this to the other host, unpack and use as normal but without sudo in the apptainer --shell command.

# Fitting Data

#### Generating Geometry Description

To run a simulation, first create a beamline, this is described in a csv. A csv for the preparation section is provided in fujii-san.csv, for the final focus in fujii-san and is based on the nominal component positions provided by a spreadsheed provided by Fujii-san. Component misalignments are set based on a .csv in survey/ to create the preparation section run the beamline creator with
    cd survey
    python extract_positions.py
this will create a file in the gmad format in the survey directory describing the beamline and beam properties as well as tunnel geometry, physics list etc. 

There are several options available for different descriptions that are useful in different cases within extract_positions.py these are grouped together inside an if with an accompanying explanation for what that combination of settings is used for. Each option is relatively self-explanatory, however it is possible to make combinations of settings that will conflict, so it is advisable to check the usage of the option in the code before deviating from the predefined sets. There are a several options that make use of additional inputs such as the option to use a magnetic fieldmap and external gdml geometry, these inputs are contained in the magnet_responses and CAD directories respectively.

The file unoptimised.gmad can then be run through BDSIM with

    source /opt/setup.sh

    bdsim --file=survey/unoptimised.gmad
to see the interactive Geant4 viewer of the geometry.

#### Optimising the Beamline

To run fits to SSEM data, make use of the ./bdsim_optimiser directory. Running of the fitter requires changes to BDSIM source, for convenience this is contained in the t2k_bdsim_modifications.diff file, this can be applied within the container and recompiled by using apply_patch.sh. The purposes of these patches is to allow for BDSIM to be rerun with different magnetic field strength settings without recreating the geometry every time and being able to directly evaluate the beam position and width at each SSEM without writing and reading from an output root file, for reference, building a geometry, running 100 protons through and running reBdsimOptics on the output takes ~2s, with this approach a simulation can be done in 0.1s allowing for reasonable fitting times.

To build the optimiser, first source the BDSIM setup script:

    source /opt/setup.sh

Then make the necessary directories in bdsim_optimiser and source the setup script, finally run make

    cd bdsim_optimiser

    mkdir bin

    mkdir libs

    mkdir output

    source setup.sh

    make

To run the executable 

    ./bin/optimise_line [optional: filename_out.root]

This will grab the geometry defined in unoptimised.gmad and attempt to optimise the beam properties and magnetic field strengths to match the simulated beam position and width at the 9 SSEM planes to the data contained in the file ssem_data/run0910216_gen.root.

The resulting parameters will be written to a file fit_results.root and if given, the first argument to optimise_line. In addition the best-fit parameters will be output in gmad format to survey/optimised.gmad. This optimised version can then be used as a basis for studies into beam Twiss parameters, beam loss etc. Further the fit_results.root file can be read by the extract_positions.py script to allow for using these optimised settings for other configurations.

#### Experimental Bayesian Optimisation

There is a highly experimental Bayesian Optimisaion framework included in bdsim_optimiser, this makes use of python bindings which can be generated by running

    bash make_python_binding.sh

and the optimisation can be run with 

    python bayesian_optimise.py

At present this optimises a mockup of a simple quad and dipole arrangement and is not applied to the full beamline.

# Custom Geometries and Fieldmaps

Magnet responses uses the field maps made of the magnets at J-PARC in xls format and converts them to a 3d point cloud in BDSIM format and is run with:
    python create_fieldmap.py
Edits to the file for each magnet are at present required but may be made automatic in future.

The CAD directory contains magnet and beamline geometries, the gdml format is used by BDSIM to insert geometries into the beamline, the current production flow for these is: use freeCAD to create the beamline components, z=0 is the center of the gap in the beamline created in the gmad file, hence it is best to make any extrusions reflected about the z=0 plane when designing the components. Each component made of a different material should be designed separately and exported as an STL, e.g. the beampipe, the beampipe vacuum, magnet yolk, support structure etc. Each of these STLs is then combined and a material assigned by the convert_stl_to_gdml.py script, making use of pyg4ometry package. The position of each component is assigned here, for complex geometries I can imagine this will get tricky, but alignment along z is the main difficult thing, best to examine the current values and work backwards, checking with the Geant4 viewer as you go.

To include such a component in the beamline, simply replace the 'type' field in fujii-san.csv with fieldmap1dgeom or fieldmap3dgeom for the component you want to model, create_beamline.py will use the fieldmap at location magnet_response/element.dat where element is the value in the 'element' field of fujii-san.csv the geometry used will similarly be element.gdml. If only the fieldmap is desired not the custom geometry use fieldmap1dsbend, fieldmap1drbend, fieldmap3dquadrupole instead. Note this fieldmap option is currently not supported by the optimiser, however can be used to generate beam orbits in the untuned case.



Author:
Charlie Naseby, c.naseby18@imperial.ac.uk
