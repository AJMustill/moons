import unhash

def heartbeat(sim_pointer):
    global glob_dclo
    global glob_planets
    global glob_npl
    global glob_darr
    global glob_log
    
    sim = sim_pointer.contents
    
    pls = glob_planets[:glob_npl]
    if glob_npl >= 2:
        for i,p in enumerate(pls):
            pp = sim.particles[p]
            rh2 = pp.d**2 * (pp.m/(3*sim.particles[0].m))**(2/3)
            for j,q in enumerate(pls):
                qq = sim.particles[q]
                if p == q:
                    continue
#                d2 = (p.x-q.x)**2 + (p.y-q.y)**2 + (p.z-q.z)**2
#                d2 = (p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y) + (p.z-q.z)*(p.z-q.z)
                dx = pp.x-qq.x
                dy = pp.y-qq.y
                dz = pp.z-qq.z
                d2 = dx*dx + dy*dy + dz*dz
                if d2 <= rh2*glob_dclo:
                    glob_darr[i][j][2] = glob_darr[i][j][1]
                    glob_darr[i][j][1] = glob_darr[i][j][0]
                    glob_darr[i][j][0] = d2
#                    print(f't = {sim.t} CE in progress {glob_darr[i][j]}')
                    if d2 > glob_darr[i][j][1] and glob_darr[i][j][1] < glob_darr[i][j][2]:
                        print(f'CE between {p} and {q} with dist '
                              f'{np.sqrt(d2)} au at {sim.t} years')
                        print(pp)
                        print(qq)
                        with open(glob_log,'a') as f:
                            print(f'CE between {p} and {q} with dist '
                                  f'{np.sqrt(d2)} au at {sim.t} years',file=f)
                            print(pp,file=f)
                            print(qq,file=f)
                else:
                    glob_darr[i][j] = [9999.9,9999.9,9999.9]