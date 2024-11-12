from mpi4py import MPI
import numpy as np
import math

PI_SQUARE = math * math.pi

# Função para gerar a fonte sísmica utilizando a wavelet de Ricker
def generateSource(s, f, dt, nt):
    for i in range(nt):
        t = i * dt
        s[i] = (i - PI_SQUARE * f * f * t * t) * math.exp(-PI_SQUARE * f * f * t *t)

# Derivada espacial com uma aproximação de 4ª ordem para x
def calcX(previousWave, x, y, z, ny, nz, dx):
    return ((-1.0/12.0) * previousWave[(x - 2) * ny * nz + y * nz + z] +
        (4.0/3.0) * previousWave[(x - 1) * ny * nz + y * nz + z] -
        (5.0/2.0) * previousWave[x * ny * nz + y * nz + z] +
        (3.0/2.0) * previousWave[(x + 1) * ny * nz + y * nz + z] -
        (1.0/12.0) * previousWave[(x + 2) * ny * nz + y * nz + z]) / (dx * dx)

def main():
    xs, ys, zs = 15, 15, 15

if __name__ == "__main__":
    main()