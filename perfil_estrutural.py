from secao import SecaoGenerica
from math import pi, sqrt
from material import Material, Aco


class PerfilEstrutural(SecaoGenerica):

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

    Wx: 'float'
        módulo elástico da seção em relação ao eixo X (horizontal)

    Wy: 'float'
        módulo elástico da seção em relação ao eixo Y (Vertical)

    Zx: 'float',
        módulo plástico da seção em relação ao eixo X (horizontal)

    Zy: 'float'
        módulo plástico da seção em relação ao eixo Y (vertical)

    xo: 'float', opcional
        coordenada x do centro de corte da seção trasnversal em relação ao
        centróide da seção

    yo: 'float', opcional
        coordenada y do centro de corte da seção trasnversal em relação ao
        centróide da seção

    Cw: 'float' (default: None)
        constante de empenamento da seção tranversal

    simetria:
        indica se a seção apresenta eixos de simetria
    """

    def __init__(self, A, Ix, Iy, J, Wx, Wy, Zx, Zy, Awx, Awy, xo, yo, Cw, material, simetria):

        if isinstance(material, list):
            material = Aco(*material)
        if isinstance(material, dict):
            material = Aco(**material)

        super().__init__(A, Ix, Iy, J, material)

        self.Wx = Wx
        self.Wy = Wy
        self.Zx = Zx
        self.Zy = Zy
        self.Awx = Awx
        self.Awy = Awy
        self.xo = xo
        self.yo = yo
        self.Cw = Cw
        self.material = material
        self.simetria = simetria

        if self.bi_simetrica:
            self.Ne = self.__Ne_bissimetria
        elif self.simetria_x and not self.simetria_y:
            self.Ne = self.__Ne_monossimetrica_x
        elif self.simetria_y and not self.simetria_x:
            self.Ne = self.__Ne_monossimetrica_y
        else:
            self.Ne = self.__Ne_assimetrica

    @property
    def rx(self):
        """ Método que determina o raio de giração em relação ao X da seção transversal """
        return sqrt(self.Ix / self.A)

    @property
    def ry(self):
        """ Método que determina o raio de giração em relação ao Y da seção transversal """
        return sqrt(self.Iy / self.A)

    @property
    def ro(self):
        """ Método que determina o raio de giração polar de inércia da seção em relação
        ao centro de cisalhamento """
        return sqrt((self.Iy + self.Ix) / self.A + self.xo ** 2 + self.yo ** 2)

    # Capacidade resistente do perfil sem consideração dos efeitos de instabilidade

    @property
    def Afy(self):
        """ Resitência axial ao escoamento"""
        return self.A * self.material.fy

    @property
    def Afu(self):
        """ Resitência axial a ruptura"""
        return self.A * self.material.fu

    @property
    def Vplx(self):
        """ Força cortante resistênte na direção X """
        return self.Awx * self.material.fy * self.material.cv

    @property
    def Vply(self):
        """ Força cortante resistênte na direção Y"""
        return self.Awy * self.material.fy * self.material.cv

    @property
    def Mplx(self):
        """ Momento de plastificação da seção em relação ao eixo X"""
        return self.Zx * self.material.fy

    @property
    def Mply(self):
        """ Momento de plastificação da seção em relação ao eixo Y"""
        return self.Zy * self.material.fy

    @property
    def Mrx(self):
        """ Momento de início de escoamento da seção em relação ao eixo X """
        return self.Wx * self.material.fy

    @property
    def Mry(self):
        """ Momento de início de escoamento da seção em relação ao eixo Y """
        return self.Wy * self.material.fy

    @property
    def simetria_x(self):
        return self.simetria[0]

    @property
    def simetria_y(self):
        return self.simetria[1]

    @property
    def bi_simetrica(self):
        return True if self.simetria_x and self.simetria_y else False

    # Métodos para definição da carga critica de flambagem de um elemento de barra de comprimentos
    # de flambagem klx, kly e klz

    def indice_esbeltez_X(self, Lb):
        """
        Retorna o indice de esbeltez de uma barra de comprimento destravado Lb
        formado pelo perfil em relação ao eixo X

        Parameter
        ---------
        Lb: 'float'
            comprimento destravado da barra

        Return
        ------

        """
        return Lb / self.ry

    def indice_esbeltez_Y(self, Lb):
        """
        Retorna o indice de esbeltez de uma barra de comprimento destravado Lb
        formado pelo perfil em relação ao eixo Y

        Parameter
        ---------
        Lb: 'float'
            comprimento destravado da barra

        Return
        ------

        """
        return Lb / self.rx

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

    def Neyz(self, kly, klz):
        """
        Método que determina a carga critica de flambagem elástica por flexo-torção para perfis monossimetrico
        com simetria em relação ao eixo Y

        Parameter
        ---------
        kly:'float'
            comprimento de flambagem por flexão em relação ao eixo longitudinal Z

        klz:'float'
            comprimento de flambagem por torção em relação ao eixo longitudinal Z

        Return
        -----
        Neyz: 'float'
            carga crítica de flambagem elástica por flexo-torção
        """
        if self.simetria_y:
            Ney = self.Ney(kly)
            Nez = self.Nez(klz)
            Neyz = (Ney + Nez) / (2 * (1 - (self.yo / self.ro) ** 2)) * \
                   (1 - sqrt(1 - 4 * Ney * Nez * (1 - (self.yo / self.ro)**2) / (Ney + Nez) ** 2))
            return Neyz
        else:
            return None

    def Nexz(self, klx, klz):
        """
        Método que determina a carga critica de flambagem elástica por flexo-torção para perfis monossimetrico
        com simetria em relação ao eixo Y

        Parameter
        ---------
        klx:'float'
            comprimento de flambagem por flexão em relação ao eixo longitudinal X

        klz:'float'
            comprimento de flambagem por torção em relação ao eixo longitudinal Z

        Return
        -----
        Nexz: 'float'
            carga crítica de flambagem elástica por flexo-torção
        """
        if self.simetria_x:
            Nex = self.Nex(klx)
            Nez = self.Nez(klz)
            Nexz = (Nex + Nez) / (2 * (1 - (self.xo / self.ro) ** 2)) * \
                   (1 - sqrt(1 - 4 * Nex * Nez * (1 - (self.xo / self.ro)**2) / (Nex + Nez) ** 2))
            return Nexz
        else:
            return None

    def __Ne_bissimetria(self, klx, kly, klz):
        """
        Método que determina a carga critica de flambagem da barra

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

    def __Ne_monossimetrica_x(self, klx, kly, klz):
        """
        Método que determina a carga critica de flambagem da barra

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
        return min(self.Ney(kly), self.Nexz(klx, klz))

    def __Ne_monossimetrica_y(self, klx, kly, klz):
        """
        Método que determina a carga critica de flambagem da barra

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
        return min(self.Nex(klx), self.Neyz(kly, klz))

    def __Ne_assimetrica(self, klx, kly, klz):
        """
        Método que determina a carga critica de flambagem da barra

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
        pass

    def fe(self, klx, kly, klz):
        """
        Método que determina a tensão critica de flambagem
        """
        return self.Ne(klx, kly, klz) / self.A

    def Mex(self, kly, klz):
        """
        Método que determina o momento fletor de flambagem lateral com torção em regime elástico
        """
        return self.ro * sqrt(self.Ney(kly) * self.Nez(klz))

    def Mey(self, klx, klz):
        """
        Método que determina o momento fletor de flambagem lateral com torção em regime elástico
        """
        return self.ro * sqrt(self.Nex(klx) * self.Nez(klz))


#class PerfilI(PerfilEstrutural):

    #def __init__(self, d, dl, tw, bfs, bfi, tfs, tfi):
        #pass


#class PerfilTuboRet(PerfilEstrutural):

    #def __init__(self, ca):
        #pass


#class PerfilTuboCir(PerfilEstrutural):
    #pass