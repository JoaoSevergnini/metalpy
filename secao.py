from math import sqrt, pi
from material import Material


class SecaoGenerica:
    """
    Esta classe define uma seção tranversal de barra de formato genérico
    de acordo com suas propriedades geométricas e seu material.

    Parameter
    ----------

    A: 'float'
        área total da seção transversal

    Ix: 'float'
        momento de inércia a flexão da seção em relação ao eixo X (horizontal)
        que passa pelo centroide da seção.

    Iy: 'float'
        momento de inércia a flexão da seção em relação ao eixo Y (Vertical)
        que passa pelo centroide da seção.

    J: 'float'
        constante de torção da seção em relação ao centróide da seção

    material: 'material'
        material que compõe a seção em relação ao eixo X (horizontal)
        que passa pelo centroide da seção.

    Wx: 'float' (default: None)
        módulo elástico da seção em relação ao eixo X (horizontal)

    Wy: 'float' (default: None)
        módulo elástico da seção em relação ao eixo Y (Vertical)

    xo: 'float' (default: None)
        coordenada x do centro de corte da seção trasnversal em relação ao
        centróide da seção

    yo: 'float' (default: None)
        coordenada y do centro de corte da seção trasnversal em relação ao
        centróide da seção

    Cw: 'float' (default: None)
        constante de empenamento da seção tranversal

    simetria:
        indica se a seção apresenta eixos de simetria
    """

    def __init__(self, A, Ix, Iy, J, material, Wx=None, Wy=None, xo=None, yo=None, Cw=None, simetria=None):
        self.A = A
        self.Ix = Ix
        self.Iy = Iy
        self.J = J
        self.Wx = Wx
        self.Wy = Wy
        self.xo = xo
        self.yo = yo
        self.Cw = Cw
        self.material = material
        self.simetria = simetria

    @property
    def rx(self):
        """
        Método que determina o raio de giração em relação ao X da seção transversal
        """
        return sqrt(self.Ix / self.A)

    @property
    def ry(self):
        """
        Método que determina o raio de giração em relação ao Y da seção transversal
        """
        return sqrt(self.Iy / self.A)

    @property
    def ro(self):
        """
        Método que determina o raio de giração polar de inércia da seção em relação
        ao centro de cisalhamento
        """
        return sqrt((self.Iy + self.Ix) / self.A + self.xo ** 2 + self.yo ** 2)

    # Métodos para definição da carga critica de flambagem de um elemento de barra de comprimentos
    # de flambagem klx, kly e klz

    def Nex(self, klx):
        """
        Método que determina a carga critica de flambagem a flexão em torno do eixo X

        Parameter
        ---------
        klx:'float'
            comprimento de flambagem por flexão em relação ao eixo X

        Return
        -----
        Nex: 'float'
            carga crítica de flambagem por flexão em relação ao eixo X
        """
        Nex = pi ** 2 * self.material.E * self.Ix / (klx ** 2)
        return Nex

    def Ney(self, kly):
        """
        Método que determina a carga critica de flambagem a flexão em torno do eixo X

        Parameter
        ---------
        kly:'float'
            comprimento de flambagem por flexão em relação ao eixo Y

        Return
        -----
        Ney: 'float'
            carga crítica de flambagem por flexão em relação ao eixo Y
        """
        Ney = pi ** 2 * self.material.E * self.Iy / (kly ** 2)
        return Ney

    def Nez(self, klz):
        """
        Método que determina a carga critica de flambagem por torção em torno do eixo longitudinal Z

        Parameter
        ---------
        klz:'float'
            comprimento de flambagem por torção em relação ao eixo longitudinal Z

        Return
        -----
        Nez: 'float'
            carga crítica de flambagem por torção em relação ao eixo longitudinal Z
        """
        Nez = (pi ** 2 * self.material.E * self.Cw / (klz ** 2) + self.material.G * self.J) / self.ro ** 2
        return Nez

    def Ne(self, klx, kly, klz):
        """
        Método que determina a carga critica de flambagem por torção em torno do eixo longitudinal Z

        Parameter
        ---------
        klx:'float'
            comprimento de flambagem por flexão em relação ao eixo x

        kly:'float'
            comprimento de flambagem por flexão em relação ao eixo Y

        klz:'float'
            comprimento de flambagem por torção em relação ao eixo longitudinal Z

        Return
        -----
        Ne: 'float'
            carga crítica de flambagem
        """
        return min(self.Nex(klx), self.Ney(kly), self.Nez(klz))
