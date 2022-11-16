import sys
from random import random
import math
import argparse

def create_grid_setup(N_x,N_y,N_z, min_mass, max_mass, scenario, snapshots,final_time, dt):
  particle_string  = str(snapshots)
  particle_string += " " + str(final_time)
  particle_string += " " + str(dt)
  particle_count = 0
  for x in range(0,N_x):
   for y in range(0,N_y):
    for z in range(0,N_z):
      xPos = (random() - 0.5)*0.9*(1.0/N_x) + x * (1.0/N_x)
      yPos = (random() - 0.5)*0.9*(1.0/N_y) + y * (1.0/N_y)
      zPos = (random() - 0.5)*0.9*(1.0/N_z) + z * (1.0/N_z)
      
      mass = random()*(max_mass-min_mass) + min_mass
      
      xVel = 0
      yVel = 0
      zVel = 0
      if scenario=="random-grid":      
        pass        
      if scenario=="no-noise":
        xPos = x * (1.0/N_x)
        yPos = y * (1.0/N_y)
        zPos = z * (1.0/N_z)
      if scenario=="shock":      
        dist = math.sqrt( (xPos-0.5)*(xPos-0.5) + (yPos-0.5)*(yPos-0.5) + (zPos-0.5)*(zPos-0.5) )
        if dist<0.1:
          xVel = xPos-0.5 / (dist+0.00001)
          yVel = yPos-0.5 / (dist+0.00001)
          zVel = zPos-0.5 / (dist+0.00001)
        
        
      particle_string += " " + str(xPos) + " " + str(yPos) + " " + str(zPos) 
      particle_string += " " + str(xVel) + " " + str(yVel) + " " + str(zVel)
      particle_string += " " + str(mass)
      
      particle_count += 1
  return particle_count,particle_string


if __name__ =="__main__":
  parser = argparse.ArgumentParser(description='Particles - initial condition generator')
  parser.add_argument("--final-time",        dest="final_time",   type=float, help="The simulation runs from 0 through final time.", required=True )
  parser.add_argument("--snapshots",         dest="snapshots",    type=float, help="The simulation writes a snapshot every snapshot time units, i.e. it does not write a snapshot after each time step. Set this value to zero to switch off any I/O.", required=True )
  parser.add_argument("--dt",                dest="dt",           type=float, help="Time step size.", required=True )
  parser.add_argument("--executable-name",   dest="executable",               help="Name of you executable, i.e. something alike ./a.out or ./assignment-code. Ensure you add the relative path, i.e. the ./ prefix, on Linux systems.", required=True )
  parser.add_argument("--min-mass",          dest="min_mass",     type=float, help="Minimal mass of particles.", required=True )  
  parser.add_argument("--max-mass",          dest="max_mass",     type=float, help="Maximal mass of particles.", required=True )  
  parser.add_argument("--dim",               dest="dim",          type=int,   help="You can create 1d, 2d and 3d setups. 1d and 2d are primarily there for debugging", default=3 )  
  parser.add_argument("--N",                 dest="N",            type=int,   help="The script generates an initial setup with NxNxN particles, or NxN particles (dim=2) or N particles (dim=1).", required=True )  
  parser.add_argument("--scenario",          dest="scenario",                 help="There are different scenarios that you can play with. Supported values are random-grid, shock and no-noise", default="random-grid" )  
  args = parser.parse_args()

  N_x = args.N
  N_y = args.N
  N_z = args.N
  
  if args.dim<3:
    N_z = 1
  if args.dim<2:
    N_y = 1

  if args.scenario!="random-grid" and args.scenario!="shock" and args.scenario!="no-noise":
    raise Exception( "ERROR: This scenario is not supported." )
  
  particle_count,particle_string = create_grid_setup(N_x,N_y,N_z, args.min_mass, args.max_mass, args.scenario, args.snapshots, args.final_time,args.dt)
  
  print( "Created setup with " + str(particle_count) + " particles"  )

  script_name = args.executable + ".sh"
  dumpFile = open( script_name, "w" )
  dumpFile.write( args.executable + " " + particle_string )

  print( "Wrote a script file to current directly which you can launch directly via " + script_name + " on a Linux terminal" )
  print( "You might have to give it executable rights first via chmod u+x " + script_name )
  print( "If you don't want to use the script directly (or work in Windows, e.g., open " + script_name + " in a text editor. You'll see the exact program invocation call" )

  
