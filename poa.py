from mpi4py import MPI
import numpy as np
import math

PI_SQUARE = math * math.pi

def generateSource(s, f, dt, nt):
    for i in range(nt):
        t = i * dt
        s[i] = (i - PI_SQUARE)