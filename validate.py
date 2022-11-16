#
# This script works only with Python 3
#
import os
import re

import create_initial_conditions
import Test

MaxParticlesInSequentialUpscalingStudies = 200


def step1():       
  test = Test.Test( "step-1.cpp" )
  
  arguments = "g++ --std=c++0x"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  _,arguments = create_initial_conditions.create_grid_setup( 2,1,1, 1, 1, "no-noise" , 0.1, 0.1, 1 )
  print( "Run code with " + arguments )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()

  search_pattern = "([-+]?\d*\.\d+|\d+), *([-+]?\d*\.\d+|\d+), *([-+]?\d*\.\d+|\d+)"
  result = test.search_for_pattern_in_output(search_pattern)
  if result!="":
    print( "Got " + result + " as last output which I interprete to be a particle position ... ok (though that does not mean that the data is correct; that's something I don't validate here)" )
  else:
    print( "Last line of output should still be the one I used in the template ... failed" )
    exit()

def step2():
  test = Test.Test( "step-2.cpp" )
  
  arguments = "g++ -O3 -fopenmp --std=c++0x"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  _,arguments = create_initial_conditions.create_grid_setup( 2,1,1, 1, 1, "no_noise" , 0.1, 0.1, 1 )
  print( "Run code with " + arguments )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok (but this does not mean that the outcome is correct)" )
  else:
    print( "Run failed: " + test.last_output )
    exit()

def step3():     
  test = Test.Test( "step-3.cpp" )
  
  arguments = "g++ --std=c++0x -fno-tree-vectorize"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  particles_counts = [11,11,11]
  _,arguments = create_initial_conditions.create_grid_setup( particles_counts[0],particles_counts[1],particles_counts[2], 1, 1, "random" , 0.1, 0.0001, 10 )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()
  no_vec_time = test.runtime
  
  arguments = "g++ -O3 -fopenmp --std=c++0x"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  particles_counts = [11,11,11]
  _,arguments = create_initial_conditions.create_grid_setup( particles_counts[0],particles_counts[1],particles_counts[2], 1, 1, "random" , 0.1, 0.0001, 10 )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()
  with_vec_time = test.runtime
  
  if no_vec_time<=with_vec_time:
    print( "Code is slower with vectorisation, so you might want to tune it ... check" )
    exit()
  else:
    print( "Code is already faster by a factor of " + str(no_vec_time/with_vec_time) + " through vectorisation but you might want to tune it further ... ok" )


def step4():     
  test = Test.Test( "step-4.cpp" )
  
  arguments = "g++ -O3  --std=c++0x -fopenmp"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  particles_counts = [11,11,11] 
  _,arguments = create_initial_conditions.create_grid_setup( particles_counts[0],particles_counts[1],particles_counts[2], 1, 1, "random" , 0.1, 0.0001, 10 )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )

  print( "Test for one core" )
  result = test.run( arguments, {"OMP_NUM_THREADS": "1"} )
  if result==1:
    print( "Run terminated successfully after " + str(test.runtime) + "s ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()
  serial_runtime = test.runtime
  
  for p in range(2,28,2):
    print( "Test for " + str(p) + " cores" )
    result = test.run( arguments, {"OMP_NUM_THREADS": str(p)} )
    speedup = serial_runtime/test.runtime
    if result==1 and speedup>p*0.8:
      print( "Run terminated successfully after " + str(test.runtime) + "s, i.e. with speedup of " + str(speedup) + " ... ok  (but you might want to tune it further)" )
    elif result==1:
      print( "Run terminated successfully after " + str(test.runtime) + "s, i.e. without expected speedup ... check runtimes" )
    else:
      print( "Run failed: " + test.last_output )
      exit()




if __name__=="__main__":
	print("""
 This is a small Python script to validate that the format of a submission
 is valid. It does not check for any correctness (though it gives some clues
 on the performance). If you use the script on Hamilton - which is the type of
 machine
 I plan to use to assess the submissions - you have to load appropriate modules.
 I use python/3.6.8 and then need an Intel compiler. The environment has recently
 been updated, so use the latest Intel compiler to get best performance:

 module load python/3.6.8 intel/2020.4;

 As always, I recommend to use a compute node rather than the login nodes. To do 
 this, you have to call something like

 salloc -N 1 -p test.q python3 validate.py test.zip;

 Disclaimer: This is only a sanity check, i.e. I'll run further tests on correctness and scalability 
 when I mark the coursework. But the script ensures that your submission is formally correct plus it
 also does some very basic checks.

""")

	try:
		print("Checking step 1")
		step1()   
		print("Checking step 2")     
		step2()        
		print("Checking step 3")
		step3()
		print("Checking step 4")
		step4()
     
	except BaseException as err:
   		print(f"Unexpected {err=}, {type(err)=}")
   		raise
