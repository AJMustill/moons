# moons
Planet–planet scattering with moons

Setup and analysis scripts for REBOUND simulations of planet--planet scattering with moons orbiting the planets

Systems are integrated without moons until the first close encounter. The simulation is then backwards-integrated 
for a brief period, then moons inserted so that forward integration with moons can resume

Generates standard REBOUND simulation archive binary files, and a .log text file containing status updates for the 
simulation (add/remove particle, close encounters, etc)

Includes the following files
 - moon_sims.ipynb           Sets up, runs and analyses the simulations
 - find_primary.py           Checks each object to find which planet it is bound to, or star, or unbound ('Galaxy')
 - globals.py                Stores some global variables needed by heartbeat.py and logged_merge.py
 - heartbeat.py              REBOUND heartbeat function to check for and log close encounters
 - logged_merge.py           Replacement for REBOUND merger routine: will log details when mergers take place
 - metasim.py                "Metasimulation" object containing the rebound simulation as well as metadata
 - parse_log.py              Reads the *.log log files
 - read_particle.py          Reads the particles details in the .log files
 - SimEvent.py               Class file for merger/close encounter events read from the .log files
 - unhash.py                 "Unhashes" the Rebound particle hash identifiers to give the human-readable name. Needs list of names to function
