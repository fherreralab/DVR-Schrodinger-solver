import numpy as np
import rovibpy


class Rovib():
    """
    This class is used to simulate the vibrational structure of a molecule.
    Given a potential energy curve, it will calculate the vibrational structure
    of the molecule and return the molecule's wavefunctions and energies for
    different states.

    Attributes
    ----------
    rmin: float
        The minimum value of the radial coordinate.
    rmax: float
        The maximum value of the radial coordinate.
    step: float
        The step size of the radial coordinate.
    mass: float
        The mass of the molecule.
    NMax: int
        The maximum number of states to be considered.
    vmax: float
        The maximum value of the potential energy curve.
    """

    def __init__(self, rmin, rmax, step, mass, NMax, vmax):
        """Initializes the Rovib class.

        Parameters
        ----------
        rmin: float
            The minimum value of the radial coordinate.
        rmax: float
            The maximum value of the radial coordinate.
        step: float
            The step size of the radial coordinate.
        mass: float
            The mass of the molecule.
        NMax: int
            The maximum number of states to be considered.
        vmax: float
            The maximum value of the potential energy curve.
        """
        self.r = []
        self.Vr = []
        self.PEC = []
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
        """Sets the parameters of the simulation.

        Parameters
        ----------
        rmin: float
            The minimum value of the radial coordinate.
        rmax: float
            The maximum value of the radial coordinate.
        step: float
            The step size of the radial coordinate.
        mass: float
            The mass of the molecule.
        NMax: int
            The maximum number of states to be considered.
        vmax: float
            The maximum value of the potential energy curve.
        """
        self.rmin = rmin
        self.rmax = rmax
        self.step = step
        self.mass = mass
        self.NMax = NMax
        self.vmax = vmax
        print("Parameters set")

    def load_PEC(self, input_file, molecule_name):
        """Loads the potential energy curve from a file.

        Parameters
        ----------
        input_file: str
            The full path of the file containing the potential energy curve.
        molecule_name: str
            The name of the molecule.
        """
        input_data = np.loadtxt(input_file)
        self.r = input_data[:,0]
        self.Vr = input_data[:,1]
        self.knots = [[self.r],[self.Vr]]
        self.molecule_name = molecule_name
        print(molecule_name + " potential energy curve loaded")

    def use_rovib(self, folder):
        """Uses the Rovib package to calculate the vibrational structure of the molecule.
        
        Parameters
        ----------
        folder: str
            The full path of the folder where the input PEC file is located.
        """
        rovibpy.rovib.use_rovib(self.rmin, self.rmax, self.step,
                                self.mass, self.NMax, self.vmax, folder)
                                #self.mass, self.NMax, self.vmax, self.knots)
        print("Rovibrational structure calculated")
        
    def get_energies(self):
        """Returns the energies of the states.

        Returns
        -------
        energies: list
            The energies of the states.
        """
        self.energies = rovibpy.rovib.rovib_energies[0]
        return self.energies

    def get_wavefunctions(self):
        """Returns the wavefunctions of the states.

        Returns
        -------
        wavefunctions: list
            The wavefunctions of the states.
        """
        self.wavefunctions = rovibpy.rovib.rovib_wavefunctions[0]
        return self.wavefunctions