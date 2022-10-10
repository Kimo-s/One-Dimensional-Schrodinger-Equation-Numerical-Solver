import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as an
from abc import abstractmethod
import scipy



class quatumParticaleSolver:

    def __init__(self, dx = 0.1, dt = 1e-4, center = 2, potential = lambda x: 0, numOfPoints = 1024, a = 1, startX=-10 , endX = 30, m = 1.0, h = 1, k = 8) -> None:
        self.dt = dt
        self.potential = potential
        self.k = k
        self.startX = startX
        self.endX = endX
        self.m = m
        self.h = h
        self.currentTimeStep = 1


        # gaussian wave packet parameters
        self.a = a
        self.center = center

        # The x interval array
        self.inputX = np.linspace(startX, endX, numOfPoints, dtype=complex)

        # Size and spatial step size
        self.dx = np.real(self.inputX[1]-self.inputX[0])
        self.N = self.inputX.shape[0]


        # The wave function variables
        self.psiX = np.copy(self.inputX)
        self.psiP = np.copy(self.inputX)
        self.V = np.array(list(map(self.potential, np.real(self.inputX))))    # The potential energy function array

        # Animation matplotlib varaiables
        self.fig, self.ax = None, None
        self.upperLimit = 1.5
        self.lowerLimit = 1.5



    def initilizePsi(self):
        gaussianPacket = lambda x: ((2/(self.a*self.a*np.pi))**(1/4)*np.exp(1j*self.k*(x-self.center))*np.exp(-(x-self.center)**2/(self.a*self.a)))
        gaussianPacketMomentum = lambda x: ((self.a/(2*np.pi))**(1/4))*np.exp((-self.a/4)*(x-self.k)**2)

        for i in range(self.N):
            self.psiX[i] = gaussianPacket(self.psiX[i])
            self.psiP[i] = gaussianPacketMomentum(self.psiP[i])

    @abstractmethod
    def computeNextTimeStep(self):
        pass


    def animate(self, i):

        self.ax.clear()
        self.ax.set_xlim(self.inputX[0],self.inputX[-1])
        self.ax.set_ylim(self.lowerLimit, self.upperLimit)

        self.computeNextTimeStep()    # compute the next time step of the wavefunction and plot it
        Psi = self.psiX

        line = self.ax.plot(np.real(self.inputX), Psi.real, label = "real(ψ)")
        self.ax.plot(np.real(self.inputX), Psi.imag, label = "im(ψ)")
        self.ax.plot(np.real(self.inputX), np.absolute(Psi), label = "|ψ|^2")

        V = np.array(list(map(self.potential, np.real(self.inputX))))
        self.ax.plot(np.real(self.inputX), V, label = "potential function V(x)")

        plt.legend()
        return line


    def draw(self, frames: int = 300, videoName: str = "schrodingerEquationSolution.mp4", upperYLimit: float = 1.5, lowerYLimit: float = -1.5):
        '''
        Args:
            frames: int = 300 -> How many frames to animate 
            videoName: str = "schrodingerEquationSolution.mp4" -> The name of the output video file
            upperYLimit and lowerYLimit -> The lower and upper limit on the plot

        Return:
            None
        '''
        plt.rc('lines')
        self.fig, self.ax = plt.subplots()

        print(f"X step = {self.dx}, time step = {self.dt}")
        self.upperLimit = upperYLimit
        self.lowerLimit = lowerYLimit

        animator = an.FuncAnimation(self.fig, func=self.animate,
                            frames = frames,
                            interval = 100,
                            blit = False)

        animator.save(videoName, writer = 'ffmpeg', fps = 30)








class quatumParticaleTridig(quatumParticaleSolver):

    def __init__(self, dx = 0.1, dt = 1e-4, center = 2, potential = lambda x: 0, numOfPoints = 1024, a = 1, startX=-10 , endX = 30, m = 1.0, h = 1, k = 8) -> None:
        super(quatumParticaleTridig, self).__init__(dx, dt, center, potential, numOfPoints, a, startX, endX, m, h, k)

        V = np.array(list(map(self.potential, np.real(self.inputX))))
        
        inital = np.ones((self.N),complex)
        alp = (1j)*self.dt/(2*self.dx**2)*inital 
        xi = inital + 1j*self.dt/2*(2/(self.dx**2)*inital + V) 
        diags = np.array([-1,0,+1]) # positions of the vectors in the matrix
        vecs1 = np.array([-alp,xi,-alp])
        self.U1 = scipy.sparse.spdiags(vecs1,diags,self.N,self.N) 
        self.U1 = self.U1.tocsc()
        self.LU = scipy.sparse.linalg.splu(self.U1) 


        self.U2 = np.zeros((self.N,self.N), complex)
        # self.U1 = np.zeros((self.N,self.N), complex)

        alpha = 1j*self.dt/(2*self.dx**2)
        

        # U1mid = np.copy(1+(1j*self.dt/2)*((2/self.dx**2)+V))
        # U1Up = np.full(self.N, -alpha)
        # U1Down = np.full(self.N, -alpha)


        for i in range(self.N):
            # U1[i,i] = 1+(1j*self.dt/2)*((2/self.dx**2)+V[i])
            self.U2[i,i] = 1-(1j*self.dt/2)*((2/self.dx**2)+V[i])
            if i != 0 and i != self.N-1:
                self.U2[i+1,i] = alpha
                self.U2[i-1,i] = alpha
                # U1[i+1,i] = -alpha
                # U1[i-1,i] = -alpha
        self.U2[1,0] = alpha
        self.U2[self.N-2, self.N-1] = alpha
        # U1[1,0] = -alpha
        # U1[self.N-2, self.N-1] = -alpha





    def computeNextTimeStep(self):
        # V = np.array(list(map(self.potential, np.real(self.inputX))))

        # U2 = np.zeros((self.N,self.N), complex)
        # U1 = np.zeros((self.N,self.N), complex)

        # alpha = 1j*self.dt/(2*self.dx**2)
        

        # U1mid = np.copy(1+(1j*self.dt/2)*((2/self.dx**2)+V))
        # U1Up = np.full(self.N, -alpha)
        # U1Down = np.full(self.N, -alpha)


        # for i in range(self.N):
        #     U1[i,i] = 1+(1j*self.dt/2)*((2/self.dx**2)+V[i])
        #     U2[i,i] = 1-(1j*self.dt/2)*((2/self.dx**2)+V[i])
        #     if i != 0 and i != self.N-1:
        #         U2[i+1,i] = alpha
        #         U2[i-1,i] = alpha
        #         U1[i+1,i] = -alpha
        #         U1[i-1,i] = -alpha
        # U2[1,0] = alpha
        # U2[self.N-2, self.N-1] = alpha
        # U1[1,0] = -alpha
        # U1[self.N-2, self.N-1] = -alpha


        psiTemp = self.U2@self.psiX

        # self.psiX = self.TDMAsolver(U1Down, U1mid, U1Up, psiTemp)
        self.psiX = self.LU.solve(psiTemp)


        print(f"Time step {self.currentTimeStep} finished, Time = {self.currentTimeStep*self.dt}")
        self.currentTimeStep += 1
        return self.psiX


    def TDMAsolver(self, a, b, c, d):
        '''
        TDMA solver, a b c d can be NumPy array type or Python list type.
        refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
        and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
        '''
        nf = len(d) # number of equations
        ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
        for it in range(1, nf):
            mc = ac[it-1]/bc[it-1]
            bc[it] = bc[it] - mc*cc[it-1] 
            dc[it] = dc[it] - mc*dc[it-1]
                    
        xc = bc
        xc[-1] = dc[-1]/bc[-1]

        for il in range(nf-2, -1, -1):
            xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

        return xc









# None working examples of the solver. Dont use.

# class quatumParticaleTridiagonal(quatumParticaleSolver):

#     def computeNextTimeStep(self):
#         H = self.computeTridiagonalMatrix()
#         # print(np.diagonal(H))
#         X = np.linalg.solve(H, self.psiX)
#         # print(X)
#         temp = np.copy(self.psiX)
#         temp = 2*X-self.psiX
#         temp = temp/np.linalg.norm(temp)
#         self.psiX = temp
#         # print(self.psiX)

#         print(f"Time step {self.currentTimeStep} finished, Time = {self.currentTimeStep*self.dt}")
#         self.currentTimeStep += 1
#         return temp

#     def computeTridiagonalMatrix(self):

#         V = 1j*(self.dt/(2*self.h))*np.diag(np.array(list(map(self.potential, np.real(self.inputX)))))

#         K = np.zeros((self.N,self.N))
#         for i in range(self.N):
#             K[i,i] = 2
#             if i != 0 and i != self.N-1:
#                 K[i+1,i] = -1
#                 K[i-1,i] = -1
#         K[1,0] = -1
#         K[self.N-2, self.N-1] = -1

#         K = 1j*((self.h*self.dt)/(4*self.m*self.dx**2))*K
#         # print(np.diagonal(K))
#         H = 1-K+V

#         return H
        

# class quatumParticaleVideo(quatumParticaleSolver):

#     def computeNextTimeStep(self):

#         Tconj, t= self.computeTridiagonalMatrix()
#         Tconj = np.conjugate(np.transpose(Tconj))

#         Tpsi = np.matmul(Tconj, self.psiX)


#         a = np.zeros(self.N, dtype=complex)
#         for i in range(1,self.N):
#             a[i] = -(1/(t[i]+a[i-1]))

#         for i in range(1, self.N):
#             self.prevb[i] = -Tpsi[i] + a[i-1]*self.prevb[i]

#         temp = np.copy(self.psiX)
#         for i in range(self.N-2, 0, -1):
#             temp[i] = a[i]*(temp[i+1]-self.prevb[i])

#         self.psiX = np.copy(temp)/np.linalg.norm(temp)
#         print(f"Time step {self.currentTimeStep} finished, Time = {self.currentTimeStep*self.dt}")
#         self.currentTimeStep += 1
#         return self.psiX


#     def computeTridiagonalMatrix(self):
#         V = np.array(list(map(self.potential, np.real(self.inputX))))
#         V[0] = 100000
#         V[self.N-1] = 100000


#         K = np.zeros((self.N,self.N), dtype=complex)
#         for i in range(self.N):
#             K[i,i] = 1j*(4*self.m*self.dx**2/(self.h*self.dt))-((2*self.m*self.dx**2)/(self.h**2))*V[i]-2
#             if i != 0 and i != self.N-1:
#                 K[i+1,i] = 1
#                 K[i-1,i] = 1
#         K[1,0] = 1
#         K[self.N-2, self.N-1] = 1


#         return K, np.diagonal(K)


# class quatumParticaleF(quatumParticaleSolver):

#     def computeNextTimeStep(self):
#         V = np.array(list(map(self.potential, np.real(self.inputX))))


#         temp = np.copy(self.psiX)
#         for i in range(self.N):
#             psi1 = self.psiX[i-1] if (i-1 > 0) else 0
#             psi2 = self.psiX[i]
#             psi3 = self.psiX[i+1] if (i+1 < self.N-1) else 0  
#             # print(f"{psi1}   {psi2}    {psi3}")

#             self.psiX[i] = self.psiXold[i] - 2*((1j*self.dt)/(self.h))*(V[i]*psi2- ((self.h**2)/(2*self.m*self.dx**2))*(psi3+psi1-2*psi2))


#         self.psiX = self.psiX/np.linalg.norm(self.psiX)
#         self.psiXold = np.copy(temp)
#         print(f"Time step {self.currentTimeStep} finished, Time = {self.currentTimeStep*self.dt}")
#         self.currentTimeStep += 1
#         return self.psiX
