import unhash

def heartbeat(sim_pointer):
    global glob_dclo
    global glob_archive
    global glob_is_close
    
    sim = sim_pointer.contents
        
    names = [unhash(p.hash,glob_names) for p in sim.particles[1:]]
    
    for p in names:
            pp = sim.particles[p]
            rh2 = pp.d**2 * (pp.m/(3*sim.particles[0].m))**(2/3)
            for q in names:
                qq = sim.particles[q]
                if p == q:
                    continue
                dx = pp.x-qq.x
                dy = pp.y-qq.y
                dz = pp.z-qq.z
                d2 = dx*dx + dy*dy + dz*dz
                if d2 <= rh2*glob_dclo:
                    if not glob_is_close:
                        print(f'CE between {p} and {q} with dist '
                              f'{np.sqrt(d2)} au at {sim.t} years')
                        print(sim.particles[p])
                        print(sim.particles[q])
                        sim.simulationarchive_snapshot(glob_archive)
                        glob_is_close = True