class Material:
    """
    Esta classe define um material isotrópico de acordo com as suas propriedades mecânicas.

    Parameter
    ---------
    E: 'float'
        módulo de elasticidade (módulo de young) do material

    poisson: 'float'
        coeficiente de poisson do material

    fy: 'float'
        tensão de escoamento do material

    fu: 'float'
        tensão de ruptura do material
    """

    def __init__(self, E, poisson, fy=None, fu=None):
        self.E = E
        self.G = E / (2 * (1 + poisson))
        self.poisson = poisson
        self.fy = fy
        self.fu = fu
