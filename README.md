# T2K-Neutrino-Beamline
This is a respository intended to house all code necessary to simulate the T2K neutrino beamline (preparation section) in the accelerator tracking software BDSIM (L.J. Nevay et al., BDSIM: An Accelerator Tracking Code with Particle-Matter Interactions, Computer Physics Communications 252 107200 (2020) https://doi.org/10.1016/j.cpc.2020.107200). This intends to include simulation of beam loss and comparison to existing beam loss measurements at the beamline with the goal to determine improvements in optics or collimator design to reduce beam loss particularly in the preparation section of the neutrino beamline.

To get started first install BDSIM or if available use the version provided by cvmfs. Most scripts and analysis are in Python, for this you require a pyton 3+ install (tested with 3.12.3).
To run a simulation first create a beamline, this is described in a csv, a csv for the preparation section is provided in fujii-san.csv and is based on the nominal component positions provided by a spreadsheed provided by Fujii-san. run the beamline creator with
    python create_beamline.py
this will create a file in the gmad format in the gmad directory describing the beamline and beam properties as well as tunnel geometry, physics list etc. This file can then be run through BDSIM with 
    bdsim --file=gmad/file.gmad
To see the interactive Geant4 viewer of the geometry. To run the optics analysis you can use the run_analysis.sh script, this uses the aforementioned gmad file and a user provided number of protons to simulate, the seed is set for reproducibility purposes. rebdsimOptics is then run over the resulting file to obtain the beam position, width, twiss parameters at each interface in beamline components, this is then converted using the Optics executable to a pandas dataframe csv that can easily be plotted in python.

There are a several options that make use of additional inputs such as the option to use a fieldmap and external gdml geometry, these are contained in the magnet_responses and CAD directories respectively. Magnet responses uses the field maps made of the magnets at J-PARC in xls format and converts them to a 3d point cloud in BDSIM format and is run with:
    python create_fieldmap.py
Edits to the file for each magnet are at present required but may be made automatic in future. 

The CAD directory contains magnet and beamline geometries, the gdml format is used by BDSIM to insert geometries into the beamline, the current production flow for these is: use freeCAD to create the beamline components, z=0 is the center of the gap in the beamline created in the gmad file, hence it is best to make any extrusions reflected about the z=0 plane when designing the components. Each component made of a different material should be designed separately and exported as an STL, e.g. the beampipe, the beampipe vacuum, magnet yolk, support structure etc. Each of these STLs is then combined and a material assigned by the convert_stl_to_gdml.py script, making use of pyg4ometry package. The position of each component is assigned here, for complex geometries I can imagine this will get tricky, but alignment along z is the main difficult thing, best to examine the current values and work backwards, checking with the Geant4 viewer as you go.

To include such a component in the beamline, simply replace the 'type' field in fujii-san.csv with fieldmapgeom for the component you want to model, create_beamline.py will use the fieldmap at location magnet_response/element.dat where element is the value in the 'element' field of fujii-san.csv the geometry used will similarly be element.gdml.


Author:
Charlie Naseby, c.naseby18@imperial.ac.uk
