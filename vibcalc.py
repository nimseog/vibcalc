import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

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
    
    def plot_disp(self):
        plt.plot(self.solution.t, self.solution.y[0])
        plt.xlabel('Time [s]')
        plt.ylabel('Displacement [m]')
        plt.show()

if __name__ == '__main__':
    vib = MechVib1d()
    vib.set_test_force()
    vib.solve()
    vib.plot_disp()
    