import matplotlib.pyplot as pl
import numpy as np
from scipy.integrate import odeint

np.set_printoptions(threshold=np.inf)

def abs(la, lb):
    return np.sqrt(la**2 + lb**2)

def expRot2(lz, ly, lx):
    length = abs(abs(lz, ly), lx)
    l2 = length**2

    if (length == 0):
        return np.array([[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]])
    else:
        if (length < 0.000001):
            print('problem')
        m = (1 - np.cos(length)) / length**2
        n = (length + np.sin(length)) / length**3
        return np.array([[1 - m * (lx**2 + ly**2), lx + m * lz * ly - lx * n * l2, -ly + m * lx * lz + ly * n * l2],
                        [-lx + m * lz * ly + lx * n * l2, 1 - m * (lx**2 + lz**2), lz + m * lx * ly - lz * n * l2],
                        [ly + m * lx * lz - ly * n * l2, -lz + m * lx * ly + lz * n * l2, 1 - m * (ly**2 + lz**2)]])

def logRot2(a):
    length = np.acos((a[0][0] + a[1][1] + a[2][2] - 1) / 2.0)

    if (length == 0.0):
        return np.array([0, 0, 0])
    else:
        if (length < 0.0001):
            print('problem')
        m = (1 - np.cos(length)) / length**2
        n = (length + np.sin(length)) / length**3
        lz = (a[1][2] - a[2][1]) / (2 - 2 * length**2 * n)
        ly = (a[2][0] - a[0][2]) / (2 - 2 * length**2 * n)
        lx = (a[0][1] - a[1][0]) / (2 - 2 * length**2 * n)
        return np.array([lz, ly, lx])

def kreisel(z, dt):
    c1 = 0.1
    #c1 = 0
    c3 = 0.005
    #c3 = 0
    lz, ly, lx, dlz, dly, dlx = z

    rot = np.matmul(expRot2(dt * dlz, dt * dly, dt * dlx), expRot2(lz, ly, lx))
    nz, ny, nx = logRot2(rot)
    #rotOld = expRot2(lz, ly, lx)

    # drehreibung
    fz = 0 - c3 * np.sign(dlz) * 10
    #fz = 0 - c3 * (dlz)

    rotZ = np.matmul(rot, np.array([1, 0, 0]))
    # drehmoment durch gewichtsverteilung
    fy = c1 * rotZ[2]
    fx = c1 * -rotZ[1]

    c4 = 0.04
    c5 = 0.2
    #c5 = 0
    #c4 = 0
    '''
    fz += -c4 * dlz * c5 * (rotZ[1]**2 + rotZ[2]**2)
    fy += -c4 * dlz * rotZ[1]
    fx += -c4 * dlz * rotZ[2]
    '''
    # nicht triviale reibung da center of mass nicht in der mitte ist
    fz += -c4 * 10 * np.sign(dlz) * c5 * (rotZ[1]**2 + rotZ[2]**2)
    fy += -c4 * 10 * np.sign(dlz) * rotZ[1]
    fx += -c4 * 10 * np.sign(dlz) * rotZ[2]

    # the bigger c6 the smaller the cutoff has to be
    if (dly**2 + dlx**2 > 0.00000001):
        c6 = 0.002
        #c6 = 0
        # rollreibung
        fy += -c6 * dly / np.sqrt(dly**2 + dlx**2)
        fx += -c6 * dlx / np.sqrt(dly**2 + dlx**2)

    return np.array([(nz - lz) / dt, (ny - ly) / dt, (nx - lx) / dt, fz, fy, fx])

def solve(res, dt, z0, func):
    t = np.linspace(0, 1, res + 1)
    z = np.eye(res + 1, len(z0))
    #dz = np.eye(res + 1, z0.size())
    z[0] = z0
    #dz[0] = kreisel(z0, 0)
    for n in range(0, res):
        z[n + 1] = z[n] + dt * func(z[n], dt)
    return t, z


# rotZ rotY rotX, angular velocity(z, y, x)
i0 = [0, 0.1, 0, 16, 0, 0]
#i0 = [0, np.pi / 2, 0, 0, 0, 0]

t2, z2 = solve(40000, 0.02, i0, kreisel)

def toDirection(z, v):
    ret = np.eye(len(z), 3)
    for n in range(0, len(z)):
        ret[n] = np.matmul(expRot2(z[n][0], z[n][1], z[n][2]), v)
        #print(abs(abs(ret[n][0], ret[n][1]), ret[n][2]))
        #print(ret[n])
    return ret

def toForce(z):
    ret = np.eye(len(z), 3)
    for n in range(0, len(z)):
        ret[n] = kreisel(z[n], 0.01)[3:6]
        #print(abs(abs(ret[n][0], ret[n][1]), ret[n][2]))
        #print(ret[n])
    return ret


pl.figure(1)
#pl.plot(t2, z2[:,3], t2, z2[:,4], t2, z2[:,5])
pl.plot(t2[0:10000], toDirection(z2[0:10000], np.array([1, 0, 0])))
#pl.plot(t2, toForce(z2))
pl.legend()
pl.savefig('z1.png')

pl.figure(3)
pl.plot(t2[10000:20000], toDirection(z2[10000:20000], np.array([1, 0, 0])))
pl.legend()
pl.savefig('z2.png')

pl.figure(4)
pl.plot(t2[20000:30000], toDirection(z2[20000:30000], np.array([1, 0, 0])))
pl.legend()
pl.savefig('z3.png')

pl.figure(5)
pl.plot(t2[30000:40000], toDirection(z2[30000:40000], np.array([1, 0, 0])))
pl.legend()
pl.savefig('z4.png')

pl.figure(6)
pl.plot(t2[0:10000], toDirection(z2[0:10000], np.array([0, 0, 1])))
pl.legend()
pl.savefig('x1.png')

pl.figure(7)
pl.plot(t2[10000:20000], toDirection(z2[10000:20000], np.array([0, 0, 1])))
pl.legend()
pl.savefig('x2.png')

pl.figure(8)
pl.plot(t2[20000:30000], toDirection(z2[20000:30000], np.array([0, 0, 1])))
pl.legend()
pl.savefig('x3.png')

pl.figure(9)
pl.plot(t2[30000:40000], toDirection(z2[30000:40000], np.array([0, 0, 1])))
pl.legend()
pl.savefig('x4.png')

pl.figure(10)
pl.plot(t2[0:100], toDirection(z2[0:100], np.array([1, 0, 0])))
pl.legend()
pl.savefig('zAnfang.png')

pl.figure(11)
pl.plot(t2[0:100], toDirection(z2[0:100], np.array([0, 0, 1])))
pl.legend()
pl.savefig('xAnfang.png')
