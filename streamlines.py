""" Streamline calculations """
import numpy as np
import matplotlib.pyplot as plt

class VelocityField:
    def __init__(self, x, y, u, v):
        xa = np.asarray(x)
        ya = np.asarray(y)
        self.x = xa if xa.ndim == 1 else xa[0]
        self.y = ya if ya.ndim == 1 else ya[:,0]
        self.u = u
        self.v = v
        self.dx = (self.x[-1]-self.x[0])/(self.x.size-1)
        self.dy = (self.y[-1]-self.y[0])/(self.y.size-1)

        # marker for which regions have contours
        self.used = np.zeros(u.shape, dtype=bool)
        self.used[0] = True
        self.used[-1] = True
        self.used[:,0] = True
        self.used[:,-1] = True

        self.markEmpty()

    def extents(self):
        return self.x[0], self.x[-1], self.y[0], self.y[-1]

    def markEmpty(self):
        for i in range(self.x.size):
            for j in range(self.y.size):
                if self.u[j,i] == 0.0 and self.v[j,i] == 0.0:
                    self.used[j,i] = True


def makeStreamline(field, x0, y0, res, spacing, maxLen, detectLoops=True):

    xmin, xmax, ymin, ymax = field.extents()
    d = res*np.sqrt(field.dx*field.dy)

    # forwards:
    sx = [x0]
    sy = [y0]

    x = x0
    y = y0
    i = 0
    while True:
        if x >= xmax or x <= xmin or y >= ymax or y <= ymin:
            break

        u, v = fixedGridInterp(field, x, y, spacing)
        theta = np.arctan2(v,u)

        x += d*np.cos(theta)
        y += d*np.sin(theta)
        sx.append(x)
        sy.append(y)

        i += 1

        # Detect closed loops and nodes
        if detectLoops and i%10 == 0:
            D = np.ndarray(len(sx)-1)
            for j in range(D.size):
                D[j] = np.hypot(x-sx[j], y-sy[j])

            if (D < 0.9*d).any():
                break

        if i > maxLen/2:
            break

    # backwards
    rx = []
    ry = []

    x = x0
    y = y0
    i = 0
    while True:
        if x >= xmax or x <= xmin or y >= ymax or y <= ymin:
            break

        u, v = fixedGridInterp(field, x, y, spacing)
        theta = np.arctan2(v,u)

        x -= d*np.cos(theta)
        y -= d*np.sin(theta)
        rx.append(x)
        ry.append(y)

        i += 1

        # Detect closed loops and nodes
        if detectLoops and  i % 10 == 0:
            D = np.ndarray(len(rx)-1)
            for j in range(D.size):
                D[j] = np.hypot(x-rx[j], y-ry[j])

            if (D < 0.9*d).any():
                break

        if i > maxLen/2:
            break

    rx.reverse()
    ry.reverse()

    return rx+sx, ry+sy


def fixedGridInterp(field, x, y, spacing=None):
    xmin, xmax, ymin, ymax = field.extents()

    i = (x-xmin)/field.dx
    ai = i % 1

    j = (y-ymin)/field.dy
    aj = j % 1

    try:
        # Bilinear interpolation
        u = (field.u[j,i]*(1-ai)*(1-aj) +
             field.u[j,i+1]*ai*(1-aj) +
             field.u[j+1,i]*(1-ai)*aj +
             field.u[j+1,i+1]*ai*aj)

        v = (field.v[j,i]*(1-ai)*(1-aj) +
             field.v[j,i+1]*ai*(1-aj) +
             field.v[j+1,i]*(1-ai)*aj +
             field.v[j+1,i+1]*ai*aj)
    except IndexError:
        print j, y, field.dy
        print ymin, ymax
        raise

    if spacing:
        field.used[j:j+spacing,i:i+spacing] = True

    return u,v


def allStreamlines(X, Y, U, V, res=0.1, spacing=2, maxLen=100, detectLoops=True):
    F = VelocityField(X, Y, U, V)

    S = []
    while not F.used.all():
        nz = np.transpose(np.logical_not(F.used).nonzero())
        S.append(makeStreamline(F, F.x[nz[0][1]], F.y[nz[0][0]], res, spacing, maxLen, detectLoops))

    return S


def drawStreamlines(S):
    plt.figure()
    plt.ioff()
    color = 'k'
    for s in S:
        foo = plt.plot(s[0],s[1],color)

    plt.axis('tight')
    plt.ion()
    plt.draw()


def animateStreamlines(SS):
    plt.figure()
    color = 'k'
    for S in SS:
        plt.ioff()
        plt.clf()
        for s in S:
            foo = plt.plot(s[0],s[1],color)
        plt.axis('tight')
        plt.ion()
        plt.draw()
        raw_input('...')
