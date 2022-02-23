from astropy import units as u
from astropy import constants as c
from matplotlib import pyplot as plt
from matplotlib import patches
from matplotlib import lines
import numpy as np
import rebound
import time

import globs
import heartbeat
import logged_merge
import SimEvent
import find_primary

class MetaSim:
    
    def __init__(self,filestem='test/test0000.',sim=None):
        
        print('Simulation '+filestem)
        print()
        
        if sim is None:
            self.filestem = filestem
            self.archive = filestem + 'bin'
            self.log = filestem + 'log'
            globs.glob_archive = self.archive
            globs.glob_log = self.log

            # observational data for EPIC 220208795 from van der Kamp et al (2021)
            Mstar = 0.85 * u.Msun
            Rstar = 0.830 * u.Rsun
            self.name_star = 'EPIC220208795'
            # present-day disc object
            Pdisc = 290 * u.d
            edisc = 0.72
            adisc = ((c.G*Mstar*Pdisc**2/(4*np.pi**2))**(1/3)).to(u.au)

            # check if saved random seed exists
            seed_file_pl = self.filestem+'pl.seed'
            try:
                with open(seed_file_pl,'r') as f:
                    seed = int(f.read())
                rng = np.random.default_rng(seed)
            except:
                seed = int(time.time()*1000)
                rng = np.random.default_rng(seed)
                with open(seed_file_pl,'w') as f:
                    print(str(seed),file=f)


            # initial planet architecture
            self.Npl = 2
            Mpl = 2 * u.Mjup * np.ones((self.Npl))
            Rpl = (((9/4*np.pi)*Mpl/(1.8*u.g/u.cm**3))**(1/3)).to(u.au)

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


            sim.automateSimulationArchive(globs.glob_archive, interval=1000.,deletefile=True)

            sim.add(m=Mstar.to_value(u.Msun),r=Rstar.to_value(u.au),hash=self.name_star)
            for j in range(self.Npl):
                sim.add(a=apl[j].to_value(u.au),e=epl[j],inc=ipl[j],omega=wpl[j],Omega=Opl[j],M=mapl[j],
                        m=Mpl[j].to_value(u.Msun),r=Rpl[j].to_value(u.au),primary=sim.particles[self.name_star],
                        hash=self.name_pl[j])

            sim.move_to_com()
            self.sim = sim
            globs.glob_planets = [n for n in self.name_pl]
            globs.glob_npl_start = len(globs.glob_planets)
            globs.glob_npl = globs.glob_npl_start
            globs.glob_darr = [[[9999.9,9999.9,9999.9],[9999.9,9999.9,9999.9]],
                               [[9999.9,9999.9,9999.9],[9999.9,9999.9,9999.9]]]
            globs.glob_is_eject = False
            globs.glob_is_stop = False
            globs.glob_names = [self.name_star] + self.name_pl


            
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
    
        while ((self.sim.t < self.tmax) and (globs.glob_is_close == False) and (globs.glob_is_eject == False)):
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
    
        globs.glob_is_stop = True

        print('Planets-only simulation complete')
        print(f'{sim_t1-sim_t0} years took {clock_t1-clock_t0} seconds')
        with open(self.log,'a') as f:
            print('Planets-only simulation complete',file=f)
            print(f'{sim_t1-sim_t0} years took {clock_t1-clock_t0} seconds',file=f)
        
        return
    
        
    def rewind(self,t=1):
        
        print('Rewinding simulation...')
        with open(self.log,'a') as f:
            print('Rewinding simulation...',file=f)
        self.sim.dt = -self.sim.dt
        tstart = self.sim.t - t
        self.sim.integrate(tstart)
        self.sim.dt = -self.sim.dt
        print(f'Simulation rewound {t} years')
        with open(self.log,'a') as f:
            print(f'Simulation rewound {t} years',file=f)

        return
    
    
    def add_moons(self):
        
        print('Adding moons...')
        with open(self.log,'a') as f:
            print('Adding moons..',file=f)
            
        # check if saved random seed exists
        seed_file_moon = self.filestem+'moon.seed'
        try:
            with open(seed_file_moon,'r') as f:
                seed = int(f.read())
            rng = np.random.default_rng(seed)
        except:
            seed = int(time.time()*1000)
            rng = np.random.default_rng(seed)
            with open(seed_file_moon,'w') as f:
                print(str(seed),file=f)
        
        # moons. Data from JPL horizons 2022-01-26
        Nmoonsppl = 4 #per planet
        amo = (np.array([4.220279893619238E+05,6.712942692154385E+05,1.070941243876099E+06,
                         1.883987339695996E+06]) * u.km).to(u.au)
        emo = np.array([3.779947596328282E-03,9.515599858678114E-03,1.437934948643515E-03,
                        7.429672504931133E-03])
        imo = np.radians([2.233261880507145E+00,2.499911845043659E+00,2.320207222240029E+00,
                          1.964066225350317E+00])
        wmo = np.radians([2.230114025429823E+02,5.871258018281303E+01,3.520404542215327E+02,
                          3.411426411826523E+01]) # arg of peri
        Omo = np.radians([3.367638701787085E+02,3.289828238020597E+02,3.401167587299205E+02,
                          3.370179889677478E+02])
        mamo = np.radians([2.411569358232406E+02,1.160555854430442E+02,1.129206855751383E+02,
                           6.172934365391181E+01])
        self.name_mo_stem = ['Io','Europa','Ganymede','Callisto']
        Rmo = (np.array([1821.49,1560.8,2631.2,2410.3])*u.km).to(u.au)
        Mmo = (np.array([5959.9155,3202.7121,9887.8328,
                         7179.2834])*u.km**3/u.s**2 / c.G).to(u.Msun) #JPL gives product GM

        Omo_offset = np.radians(rng.uniform(0,360,2)) #we rotate Omega for all moons in the same moon system by the same amount
        self.name_moons = [[m+str(i+1) for m in self.name_mo_stem] for i in range(self.Npl)]
        self.name_moons_flat = list(np.array(self.name_moons).flat)
        globs.glob_names = globs.glob_names + self.name_moons_flat

        for i in range(self.Npl):
            for j in range(Nmoonsppl):
                self.sim.add(a=amo[j].to_value(u.au),e=emo[j],inc=imo[j],omega=wmo[j],Omega=Omo[j]+Omo_offset[i],
                             M=mamo[j],m=Mmo[j].to_value(u.Msun),r=Rmo[j].to_value(u.au),
                             primary=self.sim.particles[self.name_pl[i]],hash=self.name_moons[i][j])

        self.sim.move_to_com()

        # save a snapshot
        self.sim.simulationarchive_snapshot(self.archive)
        self.tend = self.sim.t + self.tmoons

        print('Moons added')
        with open(self.log,'a') as f:
            print('Moons added',file=f)
        
        return
    
    def run_moons(self):
        
        print('Running moons simulation...')
        with open(self.log,'a') as f:
            print('Running moons simulation...',file=f)
            
        sim_t0 = self.sim.t
        clock_t0 = time.time()
                
        self.sim.collision = "line"
        self.sim.collision_resolve = logged_merge.logged_merge
        self.sim.track_energy_offset = 0

        self.sim.automateSimulationArchive(self.archive,interval=1.,deletefile=False)

        while self.sim.t < self.tend:
            try:
                self.sim.integrate(self.tend)
            except rebound.Escape as error:
# save at this point
                self.sim.simulationarchive_snapshot(self.archive)

                print(error)
                hashes = set()
# Rebound example just allows one body to be removed; might have more    
                for h in globs.glob_names:
                    try:
                        p = self.sim.particles[h]
                        d2 = p.x**2 + p.y**2 + p.z**2
                        if d2 > self.sim.exit_max_distance**2:
                            hashes.add(h)
                    except:
                        if verbose:
                            print(f'{h} has already been removed')

# Here, we also want to remove any moons bound to a removed planet
                for h in globs.glob_names:
                    try:
                        if rebound.hash(find_primary.find_primary(h,self.name_pl,
                                                                  globs.glob_names,self.sim)) in [ha.value for ha in hashes]:
                            hashes.add(h)
                    except:
                        if verbose:
                            print(f'{h} has already been removed')


                Ein = self.sim.calculate_energy()

                for h in hashes:
                    print(f'{h} ejected at {sim.t} years')
                    print(self.sim.particles[h])
                    with open(globs.glob_log,'a') as f:
                        print(f'{h} ejected at {sim.t} years',file=f)
                        print(self.sim.particles[h],file=f)
                    if h in globs.glob_planets:
                        globs.glob_planets.remove(h)
                        globs.glob_npl = globs.glob_npl-1
                    self.sim.remove(hash=h)
                self.sim.move_to_com()    

                Eout = self.sim.calculate_energy()

#        if sim.track_energy_offset:
                self.sim.energy_offset += (Ein-Eout)
        
        sim_t1 = self.sim.t
        clock_t1 = time.time()
    
        globs.glob_is_stop = True

        print('Moons simulation complete')
        print(f'{sim_t1-sim_t0} years took {clock_t1-clock_t0} seconds')
        with open(self.log,'a') as f:
            print('Moons simulation complete',file=f)
            print(f'{sim_t1-sim_t0} years took {clock_t1-clock_t0} seconds',file=f)

        
        return
    
    def read_particle(self,string):
        
        chunks = string.split()

        for c in chunks:
            if "m=" in c:
                m = float(c[2:])
            if "x=" in c:
                x = float(c[2:])
            if "y=" in c:
                y = float(c[2:])
            if "z=" in c:
                z = float(c[2:])
            if "vx=" in c:
                vx = float(c[3:])
            if "vy=" in c:
                vy = float(c[3:])
            if "vz=" in c:
                vz = float(c[3:])
            
        return (m,x,y,z,vx,vy,vz)

    def parse_log(self):
        
        CEs = []
        colls = []
        ejecs = []

        try:
            with open(self.log,'r') as f:
                lines = f.readlines()
        except:
            print('Error reading '+self.log)
            return
        else:
            for i,l in enumerate(lines):
                if "CE between" in l:
                    chunks = l.split()
                    CEs.append(SimEvent.SimEvent('CE',(chunks[2],chunks[4]),float(chunks[-2])))

                if "collision" in l:
                    chunks = l.split()
                    colls.append(SimEvent.SimEvent('Coll',(chunks[0],chunks[5]),float(chunks[-2])))

            return (CEs, colls, ejecs)

    def analyse(self):
        
        (CEs, colls, ejecs) = self.parse_log()
        self.CEs = CEs
        self.colls = colls
        self.ejecs = ejecs

        self.sa = rebound.SimulationArchive(self.archive)
        
        self.t0 = min([s.t for s in self.sa if len(s.particles) > 3])
        self.tend = self.sa[-1].t
        self.post_moons = np.where([s.t >= self.t0 for s in self.sa])[0]
        self.post_moons = [int(i) for i in self.post_moons] # have to do int casts explicitly
        self.t = self.sa[self.post_moons[0]].t
        self.t = [self.sa[i].t for i in self.post_moons]
        self.p1moons = []
        self.p2moons = []
        self.stmoons = []
        self.ejmoons = []
        self.p1mset = [] # these are not actually sets as I want to preserve order
        self.p2mset = []
        self.stmset = []
        self.ejmset = []
        self.mhost = [[] for m in self.name_moons_flat]

        for i in self.post_moons:
            p1m = []
            p2m = []
            stm = []
            ejm = []
            for j,m in enumerate(self.name_moons_flat):
                try:
                    pl = find_primary.find_primary(m,self.name_pl,globs.glob_names,self.sa[i])
                except:
                    pl = None
                if pl == self.name_pl[0]:
                    p1m.append(m)
                    if not m in self.p1mset:
                        self.p1mset.append(m)
                if pl == self.name_pl[1]:
                    p2m.append(m)
                    if not m in self.p2mset:
                        self.p2mset.append(m)
                if pl == self.name_star:
                    stm.append(m)
                    if not m in self.stmset:
                        self.stmset.append(m)
                if pl == 'Galaxy':
                    ejm.append(m)
                    if not m in self.ejmset:
                        self.ejmset.append(m)
                        
                self.mhost[j].append(pl)


            self.p1moons.append(p1m)
            self.p2moons.append(p2m)
            self.stmoons.append(stm)
            self.ejmoons.append(ejm)
        
        self.p1nm = [len(n) for n in self.p1moons]
        self.p2nm = [len(n) for n in self.p2moons]
        self.Nmoonmax = [len(self.p1mset),len(self.p2mset)]
        self.stnm = [len(n) for n in self.stmoons]
        self.Nmoonstar = len(self.stmset)
        self.ejnm = [len(n) for n in self.ejmoons]
        self.Nmoonej = len(self.ejmset)
        self.mhost = np.array(self.mhost)
        for i in range(len(self.mhost[0])-1):
            if (self.mhost[:,i+1] != self.mhost[:,i]).any():
                print(self.t[i],self.mhost[:,i+1])
        
        return

    def make_timeline(self):
        
        CE_col = 'red'
        col = {self.name_moons_flat[0]:'sienna',
               self.name_moons_flat[1]:'olivedrab',
               self.name_moons_flat[2]:'blue',
               self.name_moons_flat[3]:'deeppink',
               self.name_moons_flat[4]:'darksalmon',
               self.name_moons_flat[5]:'lightgreen',
               self.name_moons_flat[6]:'dodgerblue',
               self.name_moons_flat[7]:'pink',
               'Planet1':'black',
               'Planet2':'gainsboro',
               self.name_star:'yellow'}
        
        margin = 0.025 #figure coordinates
        axscale = 2
        moonwidth = (1 - (8*margin)) / (axscale+self.Nmoonmax[0]+self.Nmoonmax[1]+self.Nmoonstar+self.Nmoonej)
        axwidth = axscale*moonwidth
        
        xsize = 7
        ysize = 4
        fig = plt.figure(figsize=(xsize,ysize))
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        
        ytracker = margin
        galbox = patches.Rectangle((margin,ytracker),1-2*margin,1-2*margin,edgecolor='k',fill=None)
        ax.add_patch(galbox)
        
        ytracker += (margin + self.Nmoonej*moonwidth + axwidth)
        starbox = patches.Rectangle((2*margin,ytracker),1-4*margin,3*margin+
                                    (self.Nmoonmax[0]+self.Nmoonmax[1]+self.Nmoonstar)*moonwidth,
                                    edgecolor='k',fill=None)
        ytracker += (margin + self.Nmoonstar*moonwidth)
        ax.add_patch(starbox)

        plbox = []
        Npl = 2
        for i in range(Npl):
            plbox.append(patches.Rectangle((3*margin,ytracker),1-6*margin,self.Nmoonmax[i-1]*moonwidth,
                                           edgecolor='k',fill=None))
            ytracker += (margin + self.Nmoonmax[i-1]*moonwidth)
            ax.add_patch(plbox[i])
        
        circsize = 0.025
        starcirc = patches.Ellipse((starbox.get_x(),starbox.get_y()+0.5*starbox.get_height()),circsize,
                                   circsize*xsize/ysize,
                                   facecolor=col[self.name_star],edgecolor='k')
        ax.add_patch(starcirc)
        p1circ = patches.Ellipse((plbox[1].get_x(),plbox[1].get_y()+0.5*plbox[1].get_height()),circsize,
                                  circsize*xsize/ysize,
                                  facecolor=col['Planet1'],edgecolor='k')
        ax.add_patch(p1circ)
        p2circ = patches.Ellipse((plbox[0].get_x(),plbox[0].get_y()+0.5*plbox[0].get_height()),circsize,
                                  circsize*xsize/ysize,
                                  facecolor=col['Planet2'],edgecolor='k')
        ax.add_patch(p2circ)
            
        xstart = 4*margin
        xend = 1-4*margin
        xaxis = lines.Line2D([xstart,xend],[margin+0.5*axwidth,margin+0.5*axwidth],c='k')
        ax.add_line(xaxis)
        
        # add CE events
        CElines = []
        for CE in self.CEs:
            trel = (CE.t-self.t0)/(self.tend-self.t0)
            CElines.append(lines.Line2D([xstart+(xend-xstart)*trel,xstart+(xend-xstart)*trel],
                                        [plbox[0].get_y()+plbox[0].get_height(),plbox[1].get_y()],c=CE_col,
                                        zorder=99))
            ax.add_line(CElines[-1])
        
        #set up slots for moon timelines
        pl1slots = {}
        count = 0
        for p in self.p1mset:
            pl1slots[p] = plbox[1].get_y() + plbox[1].get_height() - count*moonwidth - margin
            count += 1
        count = 0
        pl2slots = {} 
        for p in self.p2mset:
            pl2slots[p] = plbox[0].get_y() + plbox[0].get_height() - count*moonwidth - margin
            count += 1
        stslots = {}
        count = 0
        for p in self.stmset:
            stslots[p] = plbox[0].get_y()-(count+0.5)*moonwidth
            count += 1
        ejslots = {}
        count = 0
        for p in self.ejmset:
            ejslots[p] = starbox.get_y()-(count+0.5)*moonwidth
            count+=1

        moon_y = [[] for i in self.name_moons_flat]
        
        # add the moon timelines
        trel = (np.array(self.t)-self.t0)/(self.tend-self.t0)
        print(self.mhost)

        for j,m in enumerate(self.name_moons_flat):
            for i in range(len(self.t)):
                if self.mhost[j][i] == 'Planet1':
                    moon_y[j].append(pl1slots[m])
                if self.mhost[j][i] == 'Planet2':
                    moon_y[j].append(pl2slots[m])
                if self.mhost[j][i] == self.name_star:
                    moon_y[j].append(stslots[m])
                if self.mhost[j][i] is None:
                    moon_y[j].append(np.nan)
                if self.mhost[j][i] == 'Galaxy':
                    moon_y[j].append(ejslots[m])
            line = lines.Line2D(xstart+(xend-xstart)*trel,moon_y[j],c=col[m])
            ax.add_line(line)
            
        # add collision events
        for c in self.colls:
            if 'Planet' in c.names[0] and 'Planet' in c.names[1]:
                w = xstart - plbox[0].get_x() + (c.t-self.t0)/(self.tend-self.t0)*(xend-xstart)
                if '1' in c.names[0]:
                    plbox[1].set_width(w)
                if '2' in c.names[0]:
                    plbox[0].set_width(w)
                continue
            
            try:
                yind = self.name_moons_flat.index(c.names[0])
                n1 = c.names[0]
                n2 = c.names[1]
            except:
                print('planet lost?',c.names)
            else:
                try:
                    xind = min(np.where([i is None for i in self.mhost[yind]])[0][1:])
                except:
                    yind = self.name_moons_flat.index(c.names[1]) #just in case you removed the wrong body
                    xind = min(np.where([i is None for i in self.mhost[yind]])[0][1:])
                    n1 = c.names[1]
                    n2 = c.names[0]
                xy = ((c.t-self.t0)/(self.tend-self.t0) * (xend-xstart) + xstart, moon_y[yind][xind-1])
                print(xy)
                ax.add_patch(patches.Ellipse(xy,circsize,circsize*xsize/ysize,lw=3,
                                             fc=col[n2],ec=col[n1]))
                
        plt.savefig(self.filestem+'_timeline.pdf')
                
        return