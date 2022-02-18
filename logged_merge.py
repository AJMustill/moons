def logged_merge(sim_pointer,collided_particles_index):
    
    global glob_planets
    global glob_npl
    global archive
    global glob_log
    
    sim = sim_pointer.contents # retrieve the standard simulation object
    ps = sim.particles # easy access to list of particles

# save at this point
    sim.simulationarchive_snapshot(archive)
    
    i = collided_particles_index.p1   # Note that p1 < p2 is not guaranteed.    
    j = collided_particles_index.p2 

    # Merging Logic 
    total_mass = ps[i].m + ps[j].m
    merged_planet = (ps[i] * ps[i].m + ps[j] * ps[j].m)/total_mass # conservation of momentum

    # merged radius assuming keep the same mean density
    merged_radius = (total_mass**2/(ps[i].m**2/ps[i].r**3 + ps[j].m**2/ps[j].r**3))**(1/3)

    # keep the more massive body and remove the lesser
    if ps[i].m >= ps[j].m:
        keep = i
        kill = j
        ret = 2
    else:
        keep = j
        kill = i
        ret = 1
    print(f'{unhash(ps[kill].hash,name_all)} removed by collision with '
          f'{unhash(ps[keep].hash,name_all)} at {sim.t} years')
    print(ps[i])
    print(ps[j])
    orbit = ps[kill].calculate_orbit(primary=ps[keep])
    print(f'Orbelts: a = {orbit.a} au; e = {orbit.e} ; I = {orbit.inc} ; Omega = {orbit.Omega}; '
          f'omega = {orbit.omega} ; MA = {orbit.M}')
    with open(glob_log,'a') as f:
        print(f'{unhash(ps[kill].hash,name_all)} removed by collision with '
              f'{unhash(ps[keep].hash,name_all)} at {sim.t} years',file=f)
        print(ps[i],file=f)
        print(ps[j],file=f)
        print(f'Orbelts: a = {orbit.a} au; e = {orbit.e} ; I = {orbit.inc} ; Omega = {orbit.Omega}; '
              f'omega = {orbit.omega} ; MA = {orbit.M}',file=f)
    
    Ein = sim.calculate_energy()
    
    ps[keep] = merged_planet   # update p1's state vector (mass and radius will need corrections)
    ps[keep].m = total_mass    # update to total mass
    ps[keep].r = merged_radius # update to joined radius
# We have to remove the body here to get the energy logging correct, and then tell
# REBOUND NOT to remove something as an effect of the function return value
    if unhash(ps[kill].hash,name_all) in glob_planets:
        glob_planets.remove(unhash(ps[kill].hash,name_all))
        glob_npl = glob_npl-1
    sim.remove(kill)

    Eout = sim.calculate_energy()
    
#    if sim.track_energy_offset:
    sim.energy_offset += (Ein-Eout)
        
    return 0