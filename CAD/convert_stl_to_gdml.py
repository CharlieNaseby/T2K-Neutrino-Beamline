import pyg4ometry

class magnet:
    def __init__(self, size):
        self.reg = pyg4ometry.geant4.Registry()
        self.world_solid = pyg4ometry.geant4.solid.Box("worldsolid",size[0],size[1],size[2],self.reg)
        self.world_logical = pyg4ometry.geant4.LogicalVolume(self.world_solid,"G4_AIR","worldlogical",self.reg)
        self.reg.setWorld(self.world_logical.name)

    def addStl(self, file, name, material, offset=[0, 0, 0]):
        r = pyg4ometry.stl.Reader(file, name, registry=self.reg)
        s = r.getSolid()
        mat = pyg4ometry.geant4.MaterialPredefined(material, self.reg)
        l = pyg4ometry.geant4.LogicalVolume(s, mat, name+"_lv", registry=self.reg)
        p = pyg4ometry.geant4.PhysicalVolume([0, 0, 0], offset, l, name+"_pv", self.world_logical, registry=self.reg)
    
    def save(self, filename):
        w = pyg4ometry.gdml.Writer()
        w.addDetector(self.reg)
        w.write(filename)


PQ4 = magnet([550, 550, 3622])
mark = 1782
offsetz = mark - 3622/2.
PQ4.addStl("stls/PQ4-Yolk.stl", "yolk", "G4_Fe",offset=[0, 0, offsetz])
PQ4.addStl("stls/PQ4-BeamPipe.stl","beampipe", "G4_Fe")
PQ4.addStl("stls/PQ4-Vacuum.stl","vacuum", "G4_Galactic")
PQ4.save("gdml/QPQ4.gdml")

