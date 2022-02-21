from astropy import units as u
from astropy import constants as c
import numpy as np
import rebound
import time

import globals
import heartbeat

class MetaSim:
    
    def __init__(self,filestem='test/test0000.',sim=None):
        
        print('Simulation '+filestem)
        print()
        
        if sim is None:
            self.filestem = filestem
            self.archive = filestem + 'bin'
            self.log = filestem + 'log'
            globals.glob_archive = self.archive
            globals.glob_log = self.log

            # observational data for EPIC 220208795 from van der Kamp et al (2021)
            Mstar = 0.85 * u.Msun
            Rstar = 0.830 * u.Rsun
            self.name_star = 'EPIC220208795'
            # present-day disc object
            Pdisc = 290 * u.d
            edisc = 0.72
            adisc = ((c.G*Mstar*Pdisc**2/(4*np.pi**2))**(1/3)).to(u.au)

            rng = np.random.default_rng(12676)

            # initial planet architecture
            self.Npl = 2
            Rpl = 1 * u.Rjup * np.ones((self.Npl))
            Mpl = 2 * u.Mjup * np.ones((self.Npl))
            a1 = adisc*2
            rH = (a1 * (Mpl[0]/(3*Mstar))**(1/3)).to(u.au)
            delta = rng.uniform(1,5)
            apl = np.zeros((self.Npl)) * u.au
            apl[0] = a1
            apl[1] = a1 + rH*delta
            epl = rng.uniform(0,0.1,(self.Npl))
            ipl = np.radians(rng.uniform(0,3,(self.Npl)))
            wpl = np.radians(rng.uniform(0,360,(self.Npl))) #arg of peri
            Opl = np.radians(rng.uniform(0,360,(self.Npl)))
            mapl = np.radians(rng.uniform(0,360,(self.Npl)))
            self.name_pl = ['Planet1','Planet2']

    #        globals.glob_names = [name_star] + name_pl

            # times to run for
            self.tmax = 1e6
            self.TINY = 1e-3
            self.tmoons = 1e3
            
            # create rebound simulation
            print('Creating new simulation...')
            with open(self.log,'a') as f:
                print('Creating new simulation...',file=f)
            sim = rebound.Simulation()
            sim.integrator = "ias15"
            sim.units = ('yr', 'AU', 'Msun')

            sim.exit_max_distance = 50000
            sim.track_energy_offset = 0

            sim.heartbeat = heartbeat.heartbeat


            sim.automateSimulationArchive(globals.glob_archive, interval=1000.,deletefile=True)

            sim.add(m=Mstar.to_value(u.Msun),r=Rstar.to_value(u.au),hash=self.name_star)
            for j in range(self.Npl):
                sim.add(a=apl[j].to_value(u.au),e=epl[j],inc=ipl[j],omega=wpl[j],Omega=Opl[j],M=mapl[j],
                        m=Mpl[j].to_value(u.Msun),r=Rpl[j].to_value(u.au),primary=sim.particles[self.name_star],
                        hash=self.name_pl[j])

            sim.move_to_com()
            self.sim = sim
            globals.glob_planets = [n for n in self.name_pl]
            globals.glob_npl_start = len(globals.glob_planets)
            globals.glob_npl = globals.glob_npl_start
            globals.glob_darr = [[[9999.9,9999.9,9999.9],[9999.9,9999.9,9999.9]],
                                 [[9999.9,9999.9,9999.9],[9999.9,9999.9,9999.9]]]
            globals.glob_is_eject = False
            globals.glob_is_stop = False


            
        else:
            # need to know what point we're at...
            pass
            
        return
    
    def run_planets(self):
        
        sim_t0 = self.sim.t
        clock_t0 = time.time()
        
        dt = 1
    
        print('Running planets-only simulation...')
        with open(self.log,'a') as f:
            print('Running planets-only simulation...',file=f)
    
        while ((self.sim.t < self.tmax) and (globals.glob_is_close == False) and (globals.glob_is_eject == False)):
            try:
# this is pretty ugly but seems to be the only way I can abort on a close encounter.
# I guess that because the heartbeat function is called from the C code it can't pass exceptions
# back to the main Python program?
                self.sim.integrate(self.sim.t+dt)
            except rebound.Escape as error:
                print(error)
                self.sim.remove # XXX not sure what's going on here...haven't finished the ejection bit....
                break
                
        sim_t1 = self.sim.t
        clock_t1 = time.time()
    
        globals.glob_is_stop = True

        print('Planets-only simulation complete')
        print(f'{sim_t1-sim_t0} years took {clock_t1-clock_t0} seconds')
        with open(self.log,'a') as f:
            print('Planets-only simulation complete',file=f)
            print(f'{sim_t1-sim_t0} years took {clock_t1-clock_t0} seconds',file=f)
        
        return
    
        
    def rewind(self,t=1):
        
        
        
        return
    
    
    def add_moons(self):
        
        
        return