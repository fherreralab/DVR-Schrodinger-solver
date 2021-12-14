import numpy as np
import numpy.f2py as f2py
import rovibpy


class Rovib():
    """
    
    Attributes
    ----------

    """

    def __init__(self, rmin, rmax, step, mass, NMax, vmax):
        self.r = []
        self.Vr = []
        self.rmin = rmin
        self.rmax = rmax
        self.step = step
        self.mass = mass
        self.NMax = NMax
        self.vmax = vmax
        self.energies = []
        self.wavefunctions = []
        self.molecule_name = ''

    def set_params(self, rmin, rmax, step, mass, NMax, vmax):
        self.rmin = rmin
        self.rmax = rmax
        self.step = step
        self.mass = mass
        self.NMax = NMax
        self.vmax = vmax
        print("Parameters set")

    def load_PEC(self, input_file, molecule_name):
        input_data = np.loadtxt(input_file)
        self.r = input_data[:,0]
        self.Vr = input_data[:,1]
        self.molecule_name = molecule_name
        print(molecule_name + " potential energy curve loaded")

    def use_rovib(self):
        rovibpy.rovib.use_rovib(self.rmin, self.rmax, self.step,
                                self.mass, self.NMax, self.vmax)
        print("Rovibrational structure calculated")
        
    def get_energies(self):
        self.energies = rovibpy.rovib.rovib_energies[0]
        return self.energies

    def get_wvefunctions(self):
        self.wavefunctions = rovibpy.rovib.rovib_wavefunctions[0]
        return self.wavefunctions