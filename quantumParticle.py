from solverClass import quatumParticaleTridig


if __name__ == "__main__":

    # choose all the arugments needed through the solver
    solver = quatumParticaleTridig(h = 1, dt = 0.004, m = 1, center = 0, a=0.5, k=12, startX =-10, endX= 10, potential=lambda x: 0, numOfPoints=1024)
     # Initilizes the wave function, must be done
    solver.initilizePsi()


    solver.draw(frames = 500)
