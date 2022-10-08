from re import T
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import readurl

class MechVib1d():
    def __init__(self):
        self.mass = 10.
        self.damping = 3.
        self.stiffn = 100.
        
        self.t_start = 0.
        self.t_end = 20.
        self.dt = 1e-3
        
        self.init_disp = 0
        self.init_vel = 0
        
        self.forcedata = None  # Includes time
        self.spectrum_force = None
        self.spectrum_disp = None

        self.solution = None

    def __motion_equation_derivatives(self, t, y, mass, damping, stiffn, forcefun):
        disp, velo = y[0], y[1]
        acc = (forcefun(t) - damping * velo - stiffn * disp) / float(mass)
        return np.array([velo, acc])
    
    def solve(self):
        self.solution = solve_ivp(self.__motion_equation_derivatives, [self.t_start, self.t_end], [self.init_disp, self.init_vel],
            t_eval=self.get_time_pnts(), args=[self.mass, self.damping, self.stiffn, self.forcefun])
        self.calculate_spectrum()

    def get_time_pnts(self):
        n_pnts = int((self.t_end - self.t_start) / float(self.dt))
        t, dt = np.linspace(self.t_start, self.t_end, n_pnts, retstep=True)
        self.dt = dt
        return t
    
    def forcefun(self, t):
        return np.interp(t, self.forcedata[0], self.forcedata[1])
    
    def set_test_force(self):
        t = self.get_time_pnts()
        amp = 20.
        freq = 8.
        force = amp * np.sin(2 * np.pi * freq * t)
        self.forcedata = np.array([t, force])
    
    def plot_results(self, maxfreq=None):
        fig, (ax0, ax1, ax2) = plt.subplots(3)
        
        ax0.plot(self.solution.t, self.solution.y[0])
        ax0.set(xlabel='Time [s]', ylabel='Displacement [m]')

        ax1.plot(self.forcedata[0], self.forcedata[1])
        ax1.set(xlabel='Time [s]', ylabel='Force [N]')
        
        ax2.plot(self.get_spectrum_freq('load'), abs(self.spectrum_force))
        ax2.set(xlabel='Frequency [Hz]', ylabel='Force [N]')
        if maxfreq:
            ax2.set(xlim=[0, maxfreq])

        fig.tight_layout()
        plt.show()

    def read_loading_from_url(self, url, type='seismic', separator=None):
        lines = readurl.get_url_text(url).split('\n')
    
        if type == 'seismic':
            x, y = [], []
            baseacc_to_force = -self.mass * 9.81
            for line in lines:
                if line:  # False if string is empty
                    vals = line.split(separator)
                    x.append(float(vals[0]))
                    y.append(baseacc_to_force * float(vals[1]))    
            self.forcedata = np.array([x, y])
            self.t_end = x[-1]
        else:
            raise NotImplementedError('Currently only seismic data supported.')

    # TODO: To be implemented, plot this too
    def calculate_spectrum(self):
        disp = self.solution.y[0]
        self.spectrum_disp = (2. / disp.size) * np.fft.rfft(disp)

        t = self.forcedata[0]
        if not np.all(t[1:] - t[:-1] == t[1] - t[0]):
            print('Warning: Force data possibly not equally spaced in time. Spectrum for force possibly unreliable.')
        self.spectrum_force = (2. / disp.size) * np.fft.rfft(self.forcedata[1])

    def get_spectrum_freq(self, type):
        if type == 'kinematic':
            dt = self.dt
            n = self.get_time_pnts().size
        elif type == 'load':
            dt = self.forcedata[0][1] - self.forcedata[0][0]
            n = self.forcedata[0].size
        return np.fft.rfftfreq(n, dt)

if __name__ == '__main__':
    vib = MechVib1d()
    # vib.set_test_force()
    url = 'http://www.vibrationdata.com/elcentro_NS.dat'
    vib.read_loading_from_url(url)
    vib.solve()
    vib.plot_results()
    