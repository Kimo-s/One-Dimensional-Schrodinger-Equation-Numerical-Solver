# One Dimensional Schrodinger Equation Numerical Solver
![ExampleVideo](https://user-images.githubusercontent.com/35817819/194800308-5077eef4-3daf-4ffd-84db-89439ec337c0.gif)

<p align="center">
    <src="https://user-images.githubusercontent.com/35817819/194800308-5077eef4-3daf-4ffd-84db-89439ec337c0.gif" alt="Example Video">
</p>

A python program that solves Schrödinger Equation numerically given an arbitrary defined potential function with animation using matplotlib.

Inside the quantumParticle.py file there is few lines and the class that is called quatumParticaleTridig which involves all you need to animate the time evolution. When creating a new solver class you have many parameters to choose and adjust to your liking. 

Those paremeters are as follows:
  - dt: the time step. Big values will animate it faster. A good value is 0.004.
  - potential: a lambda function that will be used to make the potential of the space.
  - startX: starting X position on the right.
  - endX: ending X position on the left.
  - m: mass of the particle (currently unused)
  - h: planck constant (currently unused)
  - numOfPoints: How many points to simulate between startX and endX. The higher the better but will take longer to run.

  ### Gaussian Wave Packet parameters
  - a: Decides the initial width of the gaussian wave packet.
  - center: decides the center of the initial wave packet.
  - k: the initial speed of the particle. Positive value means going in the positive direction of the X-axis.
  


After deciding all those parameters, you must call the initilizePsi() function on the solver so it actually computes the guassian wave packet. Now you can choose to direclty manually call computeNextTimeStep() to get the wave function after one time step which is returned by the function. On other hand, you can call the draw command that will generate a video using matplotlib.

The paremeters of the draw() function are:
  - frames: how many time steps to compute.
  - videoName: the name of the output video file. **Note: you must provide an extension at the end of the name so something like solution.mp4**
  - upperYLimit and lowerYLimit: decide the range of the graph on the Y axis.
  
  
  Feel free to email me on ksalem@clemson.edu for any suggestions or feedback. Also free free to edit this code and add to it.
