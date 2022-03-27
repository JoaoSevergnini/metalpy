from math import sqrt, pi

class SecaoGenerica:

    def __init__(self, a, ix, iy, j, material, wx=None, wy=None, xo=None, yo=None, cw=None, simetria=None):
        self.a = a
        self.ix = ix
        self.iy = iy
        self.j = j
        self.wx = wx
        self.wy = wy
        self.xo = xo
        self.yo = yo
        self.cw = cw
        self.material = material
        self.simetria = simetria

    @property
    def rx(self):
        return sqrt(self.ix / self.a)

    @property
    def ry(self):
        return sqrt(self.iy / self.a)

    @property
    def ro(self):
        return sqrt(self.ry ** 2 + self.rx ** 2 + self.xo ** 2 + self.yo ** 2)

    def calcularNex(self, klx):
        nex = pi ** 2 * self.material.e * self.ix / (klx ** 2)
        return nex

    def calcularNey(self, kly):
        ney = pi ** 2 * self.material.e * self.iy / (kly ** 2)
        return ney

    def calcularNez(self, klz):
        nez = (pi ** 2 * self.material.e * self.cw / (klz ** 2) + self.material.g * self.j) / self.ro
        return nez

    def calcularNe(self, klx, kly, klz):
        return min(self.calcularNex(klx), self.calcularNey(kly), self.calcularNez(klz))