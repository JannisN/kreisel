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
        if (length < 0.0001):
            print('problem')
        m = (1 - np.cos(length)) / length**2
        n = (length + np.sin(length)) / length**3
        return np.array([[1 - m * (lx**2 + ly**2), lx + m * lz * ly - lx * n * l2, -ly + m * lx * lz + ly * n * l2],
                        [-lx + m * lz * ly + lx * n * l2, 1 - m * (lx**2 + lz**2), lz + m * lx * ly - lz * n * l2],
                        [ly + m * lx * lz - ly * n * l2, -lz + m * lx * ly + lz * n * l2, 1 - m * (ly**2 + lz**2)]])

def logRot2(a):
    length = np.acos((a[0][0] + a[1][1] + a[2][2] - 1) / 2.0)

    #length = length % (2 * np.pi)
    #print('length', length)
    if (length == 0.0):
        return np.array([0, 0, 0])
    else:
        if (length < 0.0001):
            print('problem')
        m = (1 - np.cos(length)) / length**2
        n = (length + np.sin(length)) / length**3
        # abs weil es sonst problee mit 0 ~ -6.6e-16 gibt
        '''lz2 = np.abs(length**2 + (a[0][0] - 1) / m)
        ly2 = np.abs(length**2 + (a[1][1] - 1) / m)
        lx2 = np.abs(length**2 + (a[2][2] - 1) / m)'''
        #print(lz2, ly2, lx2)
        #return np.array([np.sqrt(lz2), np.sqrt(ly2), np.sqrt(lx2)])
        lz = (a[1][2] - a[2][1]) / (2 - 2 * length**2 * n)# * np.sqrt(lz2)
        ly = (a[2][0] - a[0][2]) / (2 - 2 * length**2 * n)# * np.sqrt(ly2)
        lx = (a[0][1] - a[1][0]) / (2 -  2 * length**2 * n)# * np.sqrt(lx2)'''
        '''
        lz = np.sign(a[1][2] - a[2][1]) * np.sign(2 - length**2 * n - lz2 * n) * np.sqrt(lz2)
        ly = np.sign(a[2][0] - a[0][2]) * np.sign(2 - length**2 * n - ly2 * n) * np.sqrt(ly2)
        lx = np.sign(a[0][1] - a[1][0]) * np.sign(2 - length**2 * n - lx2 * n) * np.sqrt(lx2)
        '''
        '''
        print('norm', np.linalg.norm(a - expRot2(lz, ly, lx)))
        if (np.linalg.norm(a - expRot2(lz, ly, lx)) > 0.01):
            length = 2 * np.pi - length
            m = (1 - np.cos(length)) / length**2
            n = (length + np.sin(length)) / length**3
            lz = (a[1][2] - a[2][1]) / (2 - 2 * length**2 * n)# * np.sqrt(lz2)
            ly = (a[2][0] - a[0][2]) / (2 - 2 * length**2 * n)# * np.sqrt(ly2)
            lx = (a[0][1] - a[1][0]) / (2 -  2 * length**2 * n)# * np.sqrt(lx2)'''
        return np.array([lz, ly, lx])

def expRot(lz, ly, lx):
    ret = np.array([
        [np.cos(abs(ly, lx)), lx * np.sin(abs(lz, lx)) / abs(lz, lx), -ly * np.sin(abs(lz, ly)) / abs(lz, ly)],
        [-lx * np.sin(abs(ly, lx)) / abs(ly, lx), np.cos(abs(lz, lx)), lz * np.sin(abs(lz, ly)) / abs(lz, ly)],
        [ly * np.sin(abs(ly, lx)) / abs(ly, lx), -lz * np.sin(abs(lz, lx)) / abs(lz, lx), np.cos(abs(lz, ly))]])
    if (abs(lz, lx) < 0.01):
        ret[0][1] = lx
        ret[2][1] = -lz
    if (abs(ly, lz) < 0.01):
        ret[0][2] = -ly
        ret[1][2] = lz
    if (abs(lx, ly) < 0.01):
        ret[2][0] = ly
        ret[1][0] = -lx
    '''if (np.abs(lz) < 0.01):
        ret[1][2] = lz
        ret[2][1] = -lz
    if (np.abs(ly) < 0.01):
        ret[0][2] = -ly
        ret[2][0] = ly
    if (np.abs(lx) < 0.01):
        ret[0][1] = lx
        ret[1][0] = -lx'''
    return ret

def logRot(a):
    ret = [a[1][2] * np.acos(a[2][2]) / np.sqrt(1 - a[2][2]**2),
        a[2][0] * np.acos(a[0][0]) / np.sqrt(1 - a[0][0]**2),
        a[0][1] * np.acos(a[1][1]) / np.sqrt(1 - a[1][1]**2)]
    #print('bla ', a[1][1])
    if (a[2][2] >= 0.99):
        ret[0] = a[1][2]
    if (a[2][2] <= -0.99):
        if (a[1][1] ** 2 <= 0.99):
            ret[0] = -a[2][1] * np.acos(a[1][1]) / np.sqrt(1 - a[1][1]**2)
        else:
            ret[0] = np.sqrt(0.5 * (np.acos(a[1][1]) ** 2 + np.acos(a[2][2]) ** 2 - np.acos(a[0][0]) ** 2))

    if (a[0][0] >= 0.99):
        ret[1] = a[2][0]
    if (a[0][0] <= -0.99):
        if (a[2][2] ** 2 <= 0.99):
            ret[1] = -a[0][2] * np.acos(a[2][2]) / np.sqrt(1 - a[2][2]**2)
        else:
            ret[1] = np.sqrt(0.5 * (np.acos(a[2][2]) ** 2 + np.acos(a[0][0]) ** 2 - np.acos(a[1][1]) ** 2))

    if (a[1][1] >= 0.99):
        ret[2] = a[0][1]
    if (a[1][1] <= -0.99):
        if (a[0][0] ** 2 <= 0.99):
            ret[2] = -a[1][0] * np.acos(a[1][1]) / np.sqrt(1 - a[1][1]**2)
        else:
            ret[2] = np.sqrt(0.5 * (np.acos(a[0][0]) ** 2 + np.acos(a[1][1]) ** 2 - np.acos(a[2][2]) ** 2))

    return ret

def func(z,t):
    x, y=z
    print(t,z)
    return [6*y, (2*t-3*x)*y/(4*y**2+1e-12)]    

def kreisel(z, dt):
    c1 = 1
    #c1 = 0
    c3 = 0.05
    #c3 = 0
    lz, ly, lx, dlz, dly, dlx = z
    '''
    rotLengthOld = abs(abs(lz, ly), lz)
    rotLengthNew = rotLengthOld % np.pi
    if (rotLengthNew >= np.pi / 2):
        rotLengthNew -= np.pi
    if (rotLengthOld > 1):
        lz = lz / rotLengthOld * rotLengthNew
        ly = ly / rotLengthOld * rotLengthNew
        lx = lx / rotLengthOld * rotLengthNew
    '''
    #lz = lz % np.pi
    #ly = ly % np.pi
    #lx = lx % np.pi
    #a = np.sin(np.sqrt(ly**2 + lx**2)) / np.sqrt(ly**2 + lx**2)
    rot = np.matmul(expRot2(dt * dlz, dt * dly, dt * dlx), expRot2(lz, ly, lx))
    nz, ny, nx = logRot2(rot)
    #print('problem ', nx)
    rotOld = expRot2(lz, ly, lx)
    fz = 0 - c3 * np.sign(dlz)
    fz = 0 - c3 * (dlz)
    rotZ = np.matmul(rot, np.array([1, 0, 0]))
    # warum geht das nicht?
    #fz = 0 - c3 * ((nz - lz) / dt)
    fy = c1 * rotZ[2]
    fx = c1 * -rotZ[1]
    '''
    # zum testen ob stabilität ein problem ist
    fz = 0
    fy = c1 * (rot[0][2] + rotOld[0][2]) / 2
    fx = c1 * -(rot[0][1] + rotOld[0][1]) / 2
    '''
    return np.array([(nz - lz) / dt, (ny - ly) / dt, (nx - lx) / dt, fz, fy, fx])
    #return [dlz, dly, dlx, 0, -ly * a * c1, -lx * a * c1]

def solve(res, dt, z0, func):
    t = np.linspace(0, 1, res + 1)
    z = np.eye(res + 1, len(z0))
    #dz = np.eye(res + 1, z0.size())
    z[0] = z0
    #dz[0] = kreisel(z0, 0)
    for n in range(0, res):
        z[n + 1] = z[n] + dt * func(z[n], dt)

        '''
        rotLengthOld = abs(abs(z[n + 1][0], z[n + 1][1]), z[n + 1][2])
        rotLengthNew = rotLengthOld % (2 * np.pi)
        z[n + 1][0] = z[n + 1][0] / rotLengthOld * rotLengthNew
        z[n + 1][1] = z[n + 1][1] / rotLengthOld * rotLengthNew
        z[n + 1][2] = z[n + 1][2] / rotLengthOld * rotLengthNew
        '''

        #print(z[n], func(z[n], dt))
    return t, z


z0=[1,2]
i0 = [0, np.pi / 2, 0, 4, 0, 0]
#i0 = [0, 0.5, 0, 4, 0, 0]

t = np.linspace(0,1,501)
xx=odeint(kreisel, i0, t)

t2, z2 = solve(2000, 0.01, i0, kreisel)

'''
print(xx)

print(logRot(expRot(0.04, 0, 0)))
print(5.0 % np.pi)

print(kreisel([0.041, 0, 0, 4, 0, 0], 0))
'''
t0 = np.linspace(0,1,101)
lz = np.linspace(0,1,101)
ly = np.linspace(0,1,101)
lx = np.linspace(0,1,101)
for n in range(101):
    dt = n / 100.0
    rot = np.matmul(expRot2(dt * i0[3], dt * i0[4], dt * i0[5]), expRot2(i0[0], i0[1], i0[2]))
    #rot = np.matmul(expRot(1.5, 0, 0), expRot(i0[0], i0[1], i0[2]))
    #print(np.matmul(rot, np.array([1.0, 0, 0])))
    #print(rot[0][0], rot[1][0], rot[2][0])
    nz, ny, nx = logRot2(rot)
    lz[n] = nz
    ly[n] = ny
    lx[n] = nx
    nz2, ny2, nx2 = logRot2(expRot2(nz, ny, nx))
    #print('-----------------')
    #print(nz, ny, nx)
    #print(nz2, ny2, nx2)

    #print(rot - expRot2(nz, ny, nx))
    '''
    lz[n] = rot[0, 2]
    ly[n] = rot[1, 2]
    lx[n] = rot[2, 2]
    '''
    '''lz[n] = rot[0][0]
    ly[n] = rot[1][0]
    lx[n] = rot[2][0]'''

#print(logRot2(expRot2(i0[0], i0[1], i0[2])))
'''
print(expRot2(*logRot2(expRot2(1, 0, 0))))
print(expRot2(1, 0, 0))
print(expRot2(*logRot2(expRot2(2, 0, 0))))
print(expRot2(2, 0, 0))
print(expRot2(*logRot2(expRot2(3, 0, 0))))
print(expRot2(3, 0, 0))
print(expRot2(*logRot2(expRot2(4, 0, 0))))
print(expRot2(4, 0, 0))
print(expRot2(*logRot2(expRot2(5, 0, 0))))
print(expRot2(5, 0, 0))
print(expRot2(*logRot2(expRot2(6, 0, 0))))
print(expRot2(6, 0, 0))
print(expRot2(*logRot2(expRot2(7, 0, 0))))
print(expRot2(7, 0, 0))
'''
'''
print(expRot2(-0.6561691023012153, 1.4693449200145123, -0.6561691023012153))
print(expRot2(-0.6870432252704635, 1.4593705250372655, -0.6870432252704635))

print(logRot2(expRot2(-1.79, -0.72, -1.76)))
print(abs(abs(-1.79, -0.72), -1.76))
print(logRot2(expRot2(1.79, -0.72, 1.76)))
print(abs(abs(1.79, -0.72), 1.76))
'''
'''
print(logRot2(expRot2(0.04, 0.2, 0.41)))
print(logRot2(expRot2(1, 1, 0)))
print(logRot2(expRot2(-1, 1, 0)))
print(logRot2(expRot2(-1, 0, -1)))

print((expRot2(0.04, 0.2, 0.41)))
print((expRot2(1, 1, 0)))
print((expRot2(-1, 1, 0)))
print((expRot2(-1, 0, -1)))
'''
'''
dt = 0.0001
rot = np.matmul(expRot2(dt * i0[3], dt * i0[4], dt * i0[5]), expRot2(i0[0], i0[1], i0[2]))
print((logRot2(rot) - np.array([i0[0], i0[1], i0[2]])) / dt)
'''
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
pl.plot(t2, toDirection(z2, np.array([1, 0, 0])))
#pl.plot(t2, toForce(z2))

#pl.plot(t0, lz, t0, ly, t0, lx)
#pl.plot(t, np.cos(np.sqrt(xx[:,1]**2 + xx[:,2]**2)))
#pl.plot(t, xx[:,0], t, xx[:,1], t, xx[:,2])
#pl.plot(t, xx[:,0] % np.pi, t, xx[:,1], t, xx[:,2])
pl.legend()
pl.show()
