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

# Derivada espacial com uma aproximação de 4ª ordem para y
def calcY(previousWave, x, y, z, ny, nz, dy):
    return ((-1.0/12.0) * previousWave[x * ny * nz + (y - 2) * nz + z] +
            (4.0/3.0) * previousWave[x * ny * nz + (y - 1) * nz + z] -
            (5.0/2.0) * previousWave[x * ny * nz + y * nz + z] +
            (4.0/3.0) * previousWave[x * ny * nz + (y + 1) * nz + z] -
            (1.0/12.0) * previousWave[x * ny * nz + (y + 2) * nz + z]) / (dy * dy)

# Derivada espacial com uma aproximação de 4ª ordem para z
def calcZ(previousWave, x, y, z, ny, nz, dz):
    return ((-1.0/12.0) * previousWave[x * ny * nz + y * nz + (z - 2)] +
            (4.0/3.0) * previousWave[x * ny * nz + y * nz + (z - 1)] -
            (5.0/2.0) * previousWave[x * ny * nz + y * nz + z] +
            (4.0/3.0) * previousWave[x * ny * nz + y * nz + (z + 1)] -
            (1.0/12.0) * previousWave[x * ny * nz + y * nz + (z + 2)]) / (dz * dz)



def main():
    xs, ys, zs = 15, 15, 15
    dx, dy, dz = 10, 10, 10

if __name__ == "__main__":
    main()