from mpi4py import MPI
import numpy as np
import math
import time

PI_SQUARE = math.pi * math.pi

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

def wavePropagation(s, c, dx, dy, dz, dt, nx, ny, nz, nt, xs, ys, zs, rank, size):
    local_nx = nx // size
    start_x = rank * local_nx
    end_x = start_x + local_nx if rank != size - 1 else nx

    previousWave = np.zeros(nx * ny * nz)
    nextWave = np.zeros(nx * ny * nz)
    u = np.zeros(nx * ny * nz)
    
    for t in range(nt):
        for x in range(2, nx - 2):
            for y in range(2, ny - 2):
                for z in range(2, nz - 2):
                    dEx = calcX(previousWave, x, y, z, ny, nz, dx)
                    dEy = calcY(previousWave, x, y, z, ny, nz, dy)
                    dEz = calcZ(previousWave, x, y, z, ny, nz, dz)

                    nextWave[x * ny * nz + y * nz + z] = c * c * dt * dt * (dEx + dEy + dEz) - \
                         previousWave[x * ny * nz + y * nz + z] + 2 * u[x * ny * nz + y * nz + z]

        # Atualizar a posição da fonte
        nextWave[xs * ny * nz + ys * nz + zs] -= c * c * dt * dt * s[t]

        # Troca de wavefields
        temp = u
        u = nextWave
        nextWave = previousWave
        previousWave = temp
        
        # if t % 50 == 0:
        #     filename = f"samples/sample_t{t}.bin"
        #     nextWave.astype(np.float32).tofile(filename)


def main():
    start_time = time.time()
    xs, ys, zs = 5, 5, 5
    dx, dy, dz = 10, 10, 10
    dt = 0.001
    nx, ny, nz = 20, 20, 20
    nt = 501
    f = 10
    c = 1500

    s = np.zeros(nt, dtype=np.float32)
    generateSource(s, f, dt, nt)
    wavePropagation(s, c, dx, dy, dz, dt, nx, ny, nz, nt, xs, ys, zs)

    end_time = time.time()
    execution_time = (end_time - start_time)
    print(f"O tempo de execução é: {execution_time:.2f} segundos")


if __name__ == "__main__":
    main()