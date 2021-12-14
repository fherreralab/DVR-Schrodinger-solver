import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go

class Visualize:
    """
    Visualize the results of the simulation.

    Attributes
    ----------
    data : numpy.ndarray
        The data to be visualized.
    """

    def __init__(self, wavefunctions, energies):
        """
        Initialize the Visualize class.

        Parameters
        ----------
        data : numpy.ndarray
            The data to be visualized.
        """
        self.wavefunctions = pd.DataFrame(data=wavefunctions)
        self.energies = energies

    def plot_data(self, rmin, rmax, step, Vr):
        """
        Plot the data.
        """
        X = np.linspace(rmin,rmax, int((rmax-rmin)/step+1))
        fig = go.Figure()
        #w = np.linspace((Vr/1000).min(), 0, 16)
        for i in range(self.df.shape[1]-1):
            fig.add_trace(go.Scatter(y=self.df.iloc[:,i]/step+self.energies[i]*219474.63067, x=X))
        fig.add_trace(go.Scatter(x=[x for x in r if x < 50], y=Vr))
        fig.update_layout(
            autosize=False,
            width=1800,
            height=800)
        fig.show()
        plotly.offline.plot({}, auto_open=True)
