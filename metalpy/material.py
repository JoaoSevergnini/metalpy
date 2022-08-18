from metalpy.util.util import NumPositivo


class Material:
    """
    Esta classe define um material isotrópico de acordo com as suas propriedades mecânicas.

    Parameter
    ---------
    E: float, opcional
        Módulo de elasticidade (módulo de young) do material

    poisson: float, opcional
        Coeficiente de poisson do material

    fy: float, opcional
        Tensão de escoamento do material

    fu: float, opcional
        Tensão de ruptura do material

    rho: float, opcional
        Massa especifica do material

    alfa = float, opcional
        Coeficiente de dilatação térmica
    """
    
    E = NumPositivo('E')
    G = NumPositivo('G')
    poisson = NumPositivo('poisson')
    rho = NumPositivo('rho')
    alfa = NumPositivo('alfa')

    def __init__(self, E, poisson, p=None, a=None):
        self.E = E
        self.G = E / (2 * (1 + poisson))
        self.poisson = poisson
        self.p = p
        self.a = a


class Aco(Material):

    def __init__(self, E, poisson, fy, fu, cv, p=None, a=None):

        super().__init__(E, poisson, p, a)

        self.fy = fy
        self.fu = fu
        self.cv = cv


class AcoNBR8800(Aco):

    def __init__(self, fy, fu):
        super().__init__(200000, 0.3, fy, fu, 0.6, 1.2E-5, 7850)
        self.G = 77000
