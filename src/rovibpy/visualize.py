import numpy as np
import pandas as pd
import plotly.graph_objs as go

class Visualize:
    """
    A class used to visualize the results of the simulation.

    Attributes
    ----------
    wavefunctions: pd.DataFrame
        The wavefunctions obtained after running the simulation.
    energies: numpy.ndarray
        The energies obtained after running the simulation.

    Methods
    -------
    plot_data(rmin, rmax, step, Vr, r)
        Plots the data given the same parameters that were used in the simulation.
    """

    def __init__(self, wavefunctions, energies):
        """Initialize the Visualize class.

        Parameters
        ----------
        wavefunctions: pd.DataFrame
            The wavefunctions obtained after running the simulation.
        energies: numpy.ndarray
            The energies obtained after running the simulation.
        """
        self.wavefunctions = pd.DataFrame(data=wavefunctions)
        self.energies = energies

    def plot_data(self, rmin, rmax, step, Vr, r):
        """Plots the data using the default renderer engine for the current environment.
        
        Based on Plotly's choice of renderer, this function will plot the data as a cell output
        in Jupyter notebooks or open a static HTML file in Python scripts, but won't save this HTML file.

        Parameters
        ----------
        rmin: float
            The minimum value of the radial coordinate.
        rmax: float
            The maximum value of the radial coordinate.
        step: float
            The step size of the radial coordinate.
        Vr: numpy.ndarray
            The potential energy curve.
        r: numpy.ndarray
            The radial coordinate.
        """
        X = np.linspace(rmin,rmax, int((rmax-rmin)/step+1))
        fig = go.Figure()
        #w = np.linspace((Vr/1000).min(), 0, 16)
        for i in range(self.wavefunctions.shape[1]-1):
            fig.add_trace(go.Scatter(y=self.wavefunctions.iloc[:,i]/step+self.energies[i]*219474.63067, x=X))
        fig.add_trace(go.Scatter(x=[x for x in r if x < 50], y=Vr))
        fig.update_layout(
            autosize=False,
            width=1800,
            height=800,
            yaxis_title="[cm⁻¹]",
            xaxis_title="[Bohr]",
        )
        fig.show()

    def plot_in_browser(self, rmin, rmax, step, Vr, r):
        """Plots the data as an interactive JavaScript object in a static HTML file that opens in the default browser.

        Parameters
        ----------
        rmin: float
            The minimum value of the radial coordinate.
        rmax: float
            The maximum value of the radial coordinate.
        step: float
            The step size of the radial coordinate.
        Vr: numpy.ndarray
            The potential energy curve.
        r: numpy.ndarray
            The radial coordinate.
        """
        X = np.linspace(rmin,rmax, int((rmax-rmin)/step+1))
        fig = go.Figure()
        #w = np.linspace((Vr/1000).min(), 0, 16)
        for i in range(self.wavefunctions.shape[1]-1):
            fig.add_trace(go.Scatter(y=self.wavefunctions.iloc[:,i]/step+self.energies[i]*219474.63067, x=X))
        fig.add_trace(go.Scatter(x=[x for x in r if x < 50], y=Vr))
        fig.update_layout(
            autosize=False,
            width=1800,
            height=800,
            yaxis_title="[cm⁻¹]",
            xaxis_title="[Bohr]",
        )
        fig.show(renderer='browser')
