import pyg4ometry

r = pyg4ometry.gdml.Reader("PQ4.gdml")

l = r.getRegistry().getWorldVolume()

v = pyg4ometry.visualisation.VtkViewerColouredMaterial()

v.addLogicalVolume(l)

v.view()
