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
    """

    def __init__(self, A, Ix, Iy, J, material):

        self.A = A
        self.Ix = Ix
        self.Iy = Iy
        self.J = J
        self.material = material

        # Rigidezes
        self.EA = self.material.E * self.A
        self.EIx = self.material.E * self.Ix
        self.EIy = self.material.E * self.Iy
        self.GJ = self.material.G * self.J



