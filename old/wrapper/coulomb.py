from numpy import *
import ctypes
solver=ctypes.CDLL("libcoulomb.so")

n = 100
x = random.rand(3*n) # coordinates
q = ones(n) / n      # charges
p = zeros(n)         # potential
f = zeros(3*n)       # force
solver.FMMcalccoulomb(n, x.ctypes.data, q.ctypes.data, p.ctypes.data, f.ctypes.data, 0) # FMM coulomb solver

# compare with direct summation in python
diffp = normp = difff = normf = 0
for i in range(n):
    pd = fx = fy = fz = 0
    for j in range(n):
        dx = x[3*i+0] - x[3*j+0]
        dy = x[3*i+1] - x[3*j+1]
        dz = x[3*i+2] - x[3*j+2]
        R2 = dx * dx + dy * dy + dz * dz
        if R2 == 0:
            invR = 0
        else:
            invR = 1 / sqrt(R2)
        invR3 = q[j] * invR * invR * invR
        pd += q[j] * invR
        fx += dx * invR3
        fy += dy * invR3
        fz += dz * invR3
    diffp = (p[i] - pd) * (p[i] - pd)
    normp = pd * pd
    difff = (f[3*i+0] - fx) * (f[3*i+0] - fx) \
          + (f[3*i+1] - fy) * (f[3*i+1] - fy) \
          + (f[3*i+2] - fz) * (f[3*i+2] - fz)
    normf = fx * fx + fy * fy + fz * fz
print ['potential error : ',sqrt(diffp/normp)]
print ['force     error : ',sqrt(difff/normf)]
