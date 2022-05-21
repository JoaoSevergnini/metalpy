from secao import SecaoGenerica
from math import pi, sqrt
from material import Material, Aco
import pandas as pd

perfis_AISC = pd.read_excel('aisc-shapes-database-v15.0.xlsx', 1).iloc[:, 84:]


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

    def __init__(self, A, Ix, Iy, J, Wxs, Wxi, Wys, Wyi, Zx, Zy, Awx, Awy, xo, yo, Cw, material, simetria,
                 tipo='GENERICO'):

        if isinstance(material, list):
            material = Aco(*material)
        if isinstance(material, dict):
            material = Aco(**material)

        super().__init__(A, Ix, Iy, J, material)

        self.Wxs = Wxs
        self.Wxi = Wxi
        self.Wys = Wys
        self.Wyi = Wyi
        self.Zx = Zx
        self.Zy = Zy
        self.Awx = Awx
        self.Awy = Awy
        self.xo = xo
        self.yo = yo
        self.Cw = Cw
        self.material = material
        self.simetria = simetria
        self.tipo = tipo


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
        return min(self.Wxs, self.Wxi) * self.material.fy

    @property
    def Mry(self):
        """ Momento de início de escoamento da seção em relação ao eixo Y """
        return min(self.Wys, self.Wyi) * self.material.fy

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

    def indice_esbeltez(self, Lbx, Lby):
        """
        Retorna o indice de esbeltez de uma barra de comprimento destravado Lb
        formado pelo perfil em relação aos eixos X e Y

        Parameter
        ---------
        Lb: 'float'
            comprimento destravado da barra

        Return
        ------

        """
        return Lbx / self.ry, Lby / self.rx

    def par_estabilidade(self, klx, kly, klz):
        """
        Método que determina as cargas críticas de flambagem e o momento de fletor de flambagem
        elástica com flexo torção para perfis monossimétricos.

        Parameter
        ---------
        klx:'float'
            comprimento de flambagem por flexão em relação ao eixo x

        kly:'float'
            comprimento de flambagem por flexão em relação ao eixo Y

        klz:'float'
            comprimento de flambagem por torção em relação ao eixo longitudinal Z

        Return
        ------

        """

        Nex = pi ** 2 * self.EIx / (klx ** 2)
        Ney = pi ** 2 * self.EIy / (kly ** 2)
        Nez = (pi ** 2 * self.material.E * self.Cw / (klz ** 2) + self.GJ) / self.ro ** 2
        Nexz = None
        Neyz = None
        Me = None

        if self.bi_simetrica:
            Ne = min(Nex, Ney, Nez)

        elif not self.bi_simetrica and self.simetria_y:

            Neyz = (Ney + Nez) / (2 * (1 - (self.yo / self.ro) ** 2)) * \
                   (1 - sqrt(1 - 4 * Ney * Nez * (1 - (self.yo / self.ro) ** 2) / (Ney + Nez) ** 2))

            Ne = min(Nex, Neyz)
            fe = Ne / self.A
            Me = self.ro * sqrt(self.Nex(klx) * self.Nez(klz))

        elif not self.bi_simetrica and self.simetria_x:

            Nexz = (Nex + Nez) / (2 * (1 - (self.xo / self.ro) ** 2)) * \
                   (1 - sqrt(1 - 4 * Nex * Nez * (1 - (self.xo / self.ro) ** 2) / (Nex + Nez) ** 2))

            Ne = min(Ney, Nexz)
            Me = self.ro * sqrt(self.Nex(klx) * self.Nez(klz))

        else:
            pass

        fe = Ne / self.A

        return {'Ne': Ne, 'fe': fe, 'Nex': Nex, 'Ney': Ney, 'Nez': Nez, 'Nexz': Nexz, 'Neyz': Neyz, 'Me': Me}

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
                   (1 - sqrt(1 - 4 * Ney * Nez * (1 - (self.yo / self.ro) ** 2) / (Ney + Nez) ** 2))
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
                   (1 - sqrt(1 - 4 * Nex * Nez * (1 - (self.xo / self.ro) ** 2) / (Nex + Nez) ** 2))
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


class PerfilI(PerfilEstrutural):
    """
    Está classe define um perfil I soldado.

    Parameters
    ----------
    d: float

    bfs: float

    bfi: float

    tw: float

    tfs: float

    tfi = float

    """

    def __init__(self, d, bfs, bfi, tw, tfs, tfi, material):

        self.d = d
        self.dl = d - tfi - tfs
        self.tw = tw
        self.bfs = bfs
        self.bfi = bfi
        self.tfs = tfs
        self.tfi = tfi

        self.esb_mesa_s = bfs / (2 * tfs)
        self.esb_mesa_i = bfi / (2 * tfi)
        self.esb_alma = self.dl / tw

        simetria = [True, True]
        simetria[1] = False if bfs == bfi or tfi == tfs else True

        super().__init__(**self.prop_geo(), material=material, simetria=simetria, tipo='I SOLDADO')

        self.Wy = self.Wys

    def prop_geo(self):
        """
        Método que determina as propriedades geométricas de um perfil I soldado bissimetrico ou
        monossimétrico em relação ao eixo Y
        """

        # Áreas
        # -----

        Am_sup = self.bfs * self.tfs
        Am_inf = self.tfi * self.bfi
        Aalma = self.dl * self.tw

        A = Am_sup + Am_inf + Aalma

        # Altura do centro geométrico da seção
        # -------------------------------------

        ycg = (Am_inf * self.tfi / 2 + Aalma * (self.tfi + self.dl / 2) + Am_sup * (self.d - self.tfs / 2)) / A

        # Momentos de inércia e constante de torção
        # -------------------------------------------

        Ix = self.bfi * self.tfi ** 3 / 12 + Am_inf * (ycg - self.tfi / 2) ** 2 + \
             self.bfs * self.tfs ** 3 / 12 + Am_sup * (ycg - (self.d - self.tfs / 2)) ** 2 + \
             self.tw * self.dl ** 2 / 12 + Aalma * (ycg - self.d / 2) ** 2

        Iy = self.tfi * self.bfi ** 3 / 12 + \
             self.tfs * self.bfs ** 3 / 12 + \
             self.dl * self.tw ** 3 / 12

        J = 1 / 3 * (self.bfs * self.tfi ** 3 + self.bfi * self.tfi ** 3 +
                     (self.d - self.tfi / 2 - self.tfs / 2) ** 3)

        # Módulos elásticos
        # -----------------

        Wxs = Ix / ycg

        Wxi = Ix / (self.d - ycg)

        Wys = Wyi = Iy / max(self.bfi, self.bfs)

        # Módulos plásticos
        # -----------------

        if Am_sup > Aalma + Am_inf:

            ypl = self.d - (Aalma + Am_inf + self.bfs * self.tfs) / 2 * self.bfs

            ys = (ypl + self.d) / 2

            yi = (Am_inf * self.tfi / 2 + (self.tfi + self.dl / 2) + self.bfs * (self.tfs - (self.d - ypl)) *
                  (ypl + self.tfi + self.dl / 2)) / (Am_inf + Aalma + self.bfs * (self.tfs - (self.d - ypl)))

        elif Am_inf > Aalma + Am_sup:

            ypl = (Aalma + Am_sup + self.bfi * self.tfi) / 2 * self.bfi

            ys = (Am_sup * (self.d - self.tfs / 2) + Aalma * (self.tfi + self.dl / 2) + self.bfi * (
                    ypl + self.tfi) / 2) / \
                 (Am_sup + Aalma + self.bfi * (self.tfi - ypl))

            yi = (self.tfi - ypl) / 2

        else:
            ypl = self.dl + self.tfi - (Am_inf - Am_sup + self.tw * self.dl) / 2 * self.tw

            ys = (Am_sup * (self.d - self.tfs / 2) + (Aalma - self.tw * (ypl - self.tfi)) * (
                    ypl + self.dl + self.tfi) / 2) \
                 / (Am_sup + (Aalma - self.tw * (ypl - self.tfi)))

            yi = (Am_inf * self.tfi / 2 + self.tw * (ypl - self.tfi) * (self.tfi + ypl) / 2) \
                 / (Am_inf + self.tw * (ypl - self.tfi))

        Zx = A * (ys - yi) / 2

        Zy = A * (Am_inf * self.bfi + Aalma * self.tw + Am_sup * self.tfs) / 4

        # Área de cisalhamento
        # -------------------

        Awx = self.dl * self.tw

        Awy = self.tfi * self.bfi + self.tfs * self.bfs

        # Cordenadas do centro de cisalhamento em relação ao centro geométrico
        # ---------------------------------------------------------------------

        xo = 0

        yo = (self.d - self.tfi / 2 + self.tfs / 2) * self.tfi * self.bfi ** 3 / 12 \
             + self.tfs * self.bfs ** 3 / 12

        # Constante de empenamento
        # ------------------------

        Cw = (((self.tfi + self.tfs) / 2) * (self.d - (self.tfs - self.tfs) / 2) / 12) * self.bfi ** 3 * \
             self.bfs ** 3 / (self.bfi ** 3 + self.bfs ** 3)

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wxs, 'Wxi': Wxi, 'Wyi': Wyi, 'Wys': Wys,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}


class PerfilILam(PerfilEstrutural):
    """
    parameters
    ----------
    nome: str

    material: Material, list, dict, str
    """

    def __init__(self, nome, material):
        perfil = perfis_AISC[perfis_AISC['EDI_Std_Nomenclature.1'] == nome]

        self.d = float(perfil['d.1'])
        self.bf = float(perfil['bf.1'])
        self.tf = float(perfil['tf.1'])
        self.tw = float(perfil['tw.1'])
        self.r = float(perfil['kdes.1']) - self.tf

        self.esb_alma = float(perfil['h/tw.1'])
        self.esb_mesa = float(perfil['bf/2tf.1'])

        simetria = [True, True]

        super().__init__(**self.prop_geo(perfil), material=material, simetria=simetria, tipo='I LAMINADO')

        self.Wx = self.Wxs
        self.Wy = self.Wys

    def prop_geo(self, perfil):
        A = float(perfil['A.1'])
        Ix = float(perfil['Ix.1']) * 1E6
        Iy = float(perfil['Iy.1']) * 1E6
        J = float(perfil['J.1']) * 1E3
        Wx = float(perfil['Sx.1']) * 1E3
        Wy = float(perfil['Sy.1']) * 1E3
        Zx = float(perfil['Zx.1']) * 1E3
        Zy = float(perfil['Zy.1']) * 1E3
        Cw = float(perfil['Cw.1']) * 1E9

        Awy = self.d * self.tw
        Awx = 2 * self.bf * self.tf

        xo = 0
        yo = 0

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wx, 'Wxi': Wx, 'Wyi': Wy, 'Wys': Wy,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}

    def _validar_nome(self, nome):
        pass


class Caixao(PerfilEstrutural):
    """
    Está classe define uma perfil caixão.

    Parameters
    ----------
    h: float

    b: float

    tw: float

    tf: float

    material: float

    """

    def __init__(self, h, b, tw, tf, material):
        self.h = h
        self.b = b
        self.tw = tw
        self.tf = tf
        self.hint = self.h - 2 * tf
        self.bint = self.b - 2 * tw

        self.esb_alma = (h - 2 * tf) / tw
        self.esb_mesa = (b - 2 * tw) / tf

        simetria = [True, True]

        super().__init__(**self.prop_geo(), material=material, simetria=simetria, tipo='CAIXAO')

        self.Wx = self.Wxs
        self.Wy = self.Wys

    def prop_geo(self):
        Amesa = self.b * self.tf
        Aalma = self.hint * self.tw

        A = 2 * Amesa + 2 * Aalma

        Ix = 2 * (self.b * self.tf ** 3 / 12 + Amesa * (self.h / 2 - self.tf / 2) ** 2 + self.tw * self.hint ** 3 / 12)
        Iy = 2 * (self.hint * self.tw ** 3 / 12 + Aalma * (self.b / 2 - self.tw / 2) ** 2 + self.tf * self.b ** 3 / 12)

        u = 2 * (self.b + self.h - self.tw - self.tf)
        K = (self.b - self.tw) * (self.h - self.tf) * (self.tf + self.tw) / u
        J = (self.tw ** 3 + self.tf ** 3) * u / 6 + 2 * K * (self.b - self.tw) * (self.h - self.tf)

        Wx = 2 * Ix / self.b
        Wy = 2 * Iy / self.b

        Zx = (self.b * self.h ** 2 - self.bint * self.hint ** 2) / 4
        Zy = (self.h * self.b ** 2 - self.hint * self.bint ** 2) / 4

        Awy = 2 * self.hint * self.tw
        Awx = 2 * self.bint * self.tf

        xo = 0

        yo = 0

        Cw = 0

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wx, 'Wxi': Wx, 'Wyi': Wy, 'Wys': Wy,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}


class TuboRet(PerfilEstrutural):

    def __init__(self, nome, material):

        perfil = perfis_AISC[perfis_AISC['EDI_Std_Nomenclature.1'] == nome]

        try:
            self.h = float(perfil['Ht.1'])
            self.b = float(perfil['B.1'])
            self.esb_alma = float(perfil['h/tdes.1'])
            self.esb_mesa = float(perfil['b/tdes.1'])

        except:
            self.h = float(perfil['OD.1'])
            self.b = float(perfil['OD.1'])
            self.esb_alma = float(perfil['D/t.1'])
            self.esb_mesa = float(perfil['D/t.1'])

        self.t = float(perfil['tdes.1'])

        simetria = [True, True]

        super().__init__(**self.prop_geo(perfil), material=material, simetria=simetria, tipo='TUBO RET')

        self.Wx = self.Wxs
        self.Wy = self.Wys

    def prop_geo(self, perfil):

        A = float(perfil['A.1'])
        Ix = float(perfil['Ix.1']) * 1E6
        Iy = float(perfil['Iy.1']) * 1E6
        J = float(perfil['J.1']) * 1E3
        Wx = float(perfil['Sx.1']) * 1E3
        Wy = float(perfil['Sy.1']) * 1E3
        Zx = float(perfil['Zx.1']) * 1E3
        Zy = float(perfil['Zy.1']) * 1E3

        Awy = 2 * (self.h - 3 * self.t) * self.t
        Awx = 2 * (self.b - 3 * self.t) * self.t

        xo = 0
        yo = 0

        Cw = 0

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wx, 'Wxi': Wx, 'Wyi': Wy, 'Wys': Wy,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}

    def _validar_nome(self, nome):
        pass


class TuboCir(PerfilEstrutural):

    def __init__(self, nome, material):

        self._validar_nome(nome)

        perfil = perfis_AISC[perfis_AISC['EDI_Std_Nomenclature.1'] == nome]

        self.D = float(perfil['OD.1'])
        self.t = float(perfil['tdes.1'])
        self.esb = float(perfil['D/t.1'])

        self.Dint = self.D - self.t

        simetria = [True, True]

        super().__init__(**self.prop_geo(perfil), material=material, simetria=simetria, tipo='TUBO CIR')

    def prop_geo(self, perfil):

        A = float(perfil['A.1'])
        Ix = float(perfil['Ix.1']) * 1E6
        Iy = float(perfil['Iy.1']) * 1E6
        J = float(perfil['J.1']) * 1E3
        Wx = float(perfil['Sx.1']) * 1E3
        Wy = float(perfil['Sy.1']) * 1E3
        Zx = float(perfil['Zx.1']) * 1E3
        Zy = float(perfil['Zy.1']) * 1E3

        Awy = 0.5 * A
        Awx = 0.5 * A

        xo = 0
        yo = 0

        Cw = 0

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wx, 'Wxi': Wx, 'Wyi': Wy, 'Wys': Wy,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}

    def _validar_nome(self, nome):
        pass
