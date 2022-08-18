from metalpy.material import Material
from metalpy.util.util import PropGeo


class SecaoGenerica:
    """
    Esta classe define uma seção tranversal de barra de formato genérico
    de acordo com suas propriedades geométricas e seu mat.

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

    mat: 'Material'
        material que compõe a seção em relação ao eixo X (horizontal)
        que passa pelo centroide da seção.
    """

    A  = PropGeo('A', 2)
    Ix = PropGeo('Ix', 4)
    Iy = PropGeo('Iy', 4)
    J  = PropGeo('J', 4)

    def __init__(self, A, Ix, Iy, J, mat: Material):

        self.A = A
        self.Ix = Ix
        self.Iy = Iy
        self.J = J
        self.mat = mat
        
        # Rigidezes
        self.EA = self.mat.E * self.A
        self.EIx = self.mat.E * self.Ix
        self.EIy = self.mat.E * self.Iy
        self.GJ = self.mat.G * self.J



