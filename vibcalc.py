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
        self.spectrum_force = {}
        self.spectrum_disp = {}

        self.solution = None
        self.acc = None

    def __motion_equation_derivatives(self, t, y, mass, damping, stiffn, forcefun):
        disp, velo = y[0], y[1]
        acc = (forcefun(t) - damping * velo - stiffn * disp) / float(mass)
        return np.array([velo, acc])
    
    def calculate_acc(self):
        disp = self.solution.y[0]
        velo = self.solution.y[1]
        t = self.solution.t
        self.acc = (self.forcefun(t) - self.damping * velo - self.stiffn * disp) / float(self.mass)
    
    def solve(self):
        self.solution = solve_ivp(self.__motion_equation_derivatives, [self.t_start, self.t_end], [self.init_disp, self.init_vel],
            t_eval=self.get_time_pnts(), args=[self.mass, self.damping, self.stiffn, self.forcefun])
        self.calculate_acc()
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
        fig, axs = plt.subplots(3, 2)
        
        axs[0, 0].plot(self.solution.t, self.solution.y[0])
        axs[0, 0].set(xlabel='Time [s]', ylabel='Displacement [m]', title='Displacement response')

        axs[1, 0].plot(self.solution.t, self.solution.y[1])
        axs[1, 0].set(xlabel='Time [s]', ylabel='Velocity [m/s]', title='Velocity response')

        axs[2, 0].plot(self.solution.t, self.acc)
        axs[2, 0].set(xlabel='Time [s]', ylabel='Acceleration [mm/s^2]', title='Acceleration response')

        axs[0, 1].plot(self.forcedata[0], self.forcedata[1])
        axs[0, 1].set(xlabel='Time [s]', ylabel='Force [N]', title='Input force')
        
        axs[1, 1].plot(self.spectrum_force['x'], abs(self.spectrum_force['y']))
        axs[1, 1].set(xlabel='Frequency [Hz]', ylabel='Force [N]', title='Spectrum of input force')
        if maxfreq:
            axs[1, 1].set(xlim=[0, maxfreq])
        
        axs[2, 1].axis('off')

        peakdisp, t_peakdisp = self.peakdisp() 
        fig.suptitle('MECHANICAL OSCILLATOR\nMass = %s kg, stiffness = %s N/m, damping = %s Ns/m\n'
            'Peak displacement: %.3f m @ %.3f s' % (self.mass, self.stiffn, self.damping, peakdisp, t_peakdisp))
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

    def calculate_spectrum(self):
        # TODO: Check correct scaling of FFT
        
        disp = self.solution.y[0]
        n_pnts = disp.size
        self.spectrum_disp['y'] = (2. / disp.size) * np.fft.rfft(disp)
        self.spectrum_disp['x'] = np.fft.rfftfreq(n_pnts, self.dt)

        force = self.forcedata[1]
        t = self.forcedata[0]
        dt = t[1] - t[0]
        if not np.all(t[1:] - t[:-1] == dt):
            t_new = np.linspace(t[0], t[-1], t.size)
            force = np.interp(t, t_new, force)
            t = t_new
        n_pnts = force.size
        self.spectrum_force['y'] = (2. / n_pnts) * np.fft.rfft(force)
        self.spectrum_force['x'] = np.fft.rfftfreq(n_pnts, dt)

    def peakdisp(self):
        i_peakdisp = np.argmax(abs(self.solution.y[0]))
        peakdisp = self.solution.y[0][i_peakdisp]
        t_peakdisp = self.solution.t[i_peakdisp]
        return peakdisp, t_peakdisp

if __name__ == '__main__':
    vib = MechVib1d()
    # vib.set_test_force()
    url = 'http://www.vibrationdata.com/elcentro_NS.dat'
    vib.read_loading_from_url(url)
    vib.solve()
    vib.plot_results()
    