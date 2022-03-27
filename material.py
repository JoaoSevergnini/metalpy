class Material:
    '''
    Classe que define um material de acordo com as suas propriedades mecânicas


    '''

    def __init__(self, e, g, poisson, fy, fu):
        self.e = e
        self.g = g
        self.poisson = poisson
        self.fy = fy
        self.fu = fu