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
        
        self.forcedata = None

        self.solution = None

    def __motion_equation_derivative(self, t, y, mass, damping, stiffn, forcefun):
        disp, velo = y[0], y[1]
        acc = (forcefun(t) - damping * velo - stiffn * disp) / float(mass)
        return np.array([velo, acc])
    
    def solve(self):
        t_vals = self.get_time_vals()
        self.solution = solve_ivp(self.__motion_equation_derivative, [self.t_start, self.t_end], [self.init_disp, self.init_vel],
            t_eval=t_vals, args=[self.mass, self.damping, self.stiffn, self.forcefun])

    def get_time_vals(self):
        n_pnts = int((self.t_end - self.t_start) / float(self.dt))
        return np.linspace(self.t_start, self.t_end, n_pnts)
    
    def forcefun(self, t):
        return np.interp(t, self.forcedata[0], self.forcedata[1])

    def set_test_force(self):
        t = self.get_time_vals()
        amp = 20.
        freq = 8.
        force = amp * np.sin(2 * np.pi * freq * t)
        self.forcedata = np.array([t, force])
    
    def plot_solution(self):
        fig, (ax0, ax1) = plt.subplots(2)
        
        ax0.plot(self.solution.t, self.solution.y[0])
        ax0.set(xlabel='Time [s]', ylabel='Displacement [m]')

        ax1.plot(self.forcedata[0], self.forcedata[1])
        ax1.set(xlabel='Time [s]', ylabel='Force [N]')
        
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
    def get_spectrum():
        pass

if __name__ == '__main__':
    vib = MechVib1d()
    # vib.set_test_force()
    url = 'http://www.vibrationdata.com/elcentro_NS.dat'
    vib.read_loading_from_url(url)
    vib.solve()
    vib.plot_solution()
    