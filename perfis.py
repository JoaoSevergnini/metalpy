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

        self.raiz_E_fy = sqrt(self.material.E / self.material.fy)
        self.raiz_fy_E = sqrt(self.material.fy / self.material.E)

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
            Me = self.ro * sqrt(Nex * Nez)

        elif not self.bi_simetrica and self.simetria_x:

            Nexz = (Nex + Nez) / (2 * (1 - (self.xo / self.ro) ** 2)) * \
                   (1 - sqrt(1 - 4 * Nex * Nez * (1 - (self.xo / self.ro) ** 2) / (Nex + Nez) ** 2))

            Ne = min(Ney, Nexz)
            Me = self.ro * sqrt(Nex * Nez)

        else:
            pass

        fe = Ne / self.A

        return {'Ne': Ne, 'fe': fe, 'Nex': Nex, 'Ney': Ney, 'Nez': Nez, 'Nexz': Nexz, 'Neyz': Neyz, 'Me': Me}


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
        simetria[1] = False if bfs != bfi or tfi != tfs else True

        super().__init__(**self.prop_geo(), material=material, simetria=simetria, tipo='I SOLDADO')

        self.Wy = self.Wys

    def prop_geo(self):
        """
        Método que determina as propriedades geométricas de um perfil I soldado bissimetrico ou
        monossimétrico em relação ao eixo Y
        """

        # Área(A)
        # -------

        Am_sup = self.bfs * self.tfs
        Am_inf = self.tfi * self.bfi
        Aalma = self.dl * self.tw

        A = Am_sup + Am_inf + Aalma

        # Altura do centro geométrico da seção (ycg)
        # ------------------------------------------

        ycg = (Am_inf * self.tfi / 2 + Aalma * (self.tfi + self.dl / 2) + Am_sup * (self.d - self.tfs / 2)) / A

        # Momentos de inércia e constante de torção (Ix, Iy e J)
        # ------------------------------------------------------

        I1x = self.bfs * self.tfs ** 3 / 12
        I2x = self.bfi * self.tfi ** 3 / 12
        Iax = self.tw * self.dl ** 3 / 12

        I1y = self.tfs * self.bfs ** 3 / 12
        I2y = self.tfi * self.bfi ** 3 / 12
        Iay = self.dl * self.tw ** 3 / 12

        d1y = ycg - (self.d - self.tfs / 2)
        d2y = ycg - self.tfi / 2
        day = ycg - self.d / 2

        Ix = I2x + Am_inf * d2y ** 2 + \
             I1x + Am_sup * d1y ** 2 + \
             Iax + Aalma * day ** 2

        Iy = I2y + I1y + Iay
        J = 1 / 3 * (self.bfs * self.tfs ** 3 + self.bfi * self.tfi ** 3 +
                     (self.d - self.tfi / 2 - self.tfs / 2) * self.tw ** 3)

        # Módulos elásticos (Wxs, Wxi, Wys e Wyi)
        # --------------------------------------

        Wxs = Ix / ycg
        Wxi = Ix / (self.d - ycg)
        Wys = Wyi = 2 * Iy / max(self.bfi, self.bfs)

        # Módulos plásticos (Zx e Zy)
        # -------------------------

        if Am_sup > Aalma + Am_inf:

            ypl = self.d - (Aalma + Am_inf + Am_sup) / (2 * self.bfs)

            ys = (ypl + self.d) / 2

            yi = (Am_inf * self.tfi / 2 + (self.tfi + self.dl / 2) + self.bfs * (self.tfs - (self.d - ypl)) *
                  (ypl + self.tfi + self.dl / 2)) / (Am_inf + Aalma + self.bfs * (self.tfs - (self.d - ypl)))

        elif Am_inf > Aalma + Am_sup:

            ypl = (Aalma + Am_sup + self.bfi * self.tfi) / (2 * self.bfi)

            ys = (Am_sup * (self.d - self.tfs / 2) + Aalma * (self.tfi + self.dl / 2) + self.bfi * (
                    ypl + self.tfi) / 2) / \
                 (Am_sup + Aalma + self.bfi * (self.tfi - ypl))

            yi = (self.tfi - ypl) / 2

        else:
            ypl = self.dl + self.tfi - (Am_inf - Am_sup + self.tw * self.dl) / (2 * self.tw)

            ys = (Am_sup * (self.d - self.tfs / 2) + (Aalma - self.tw * (ypl - self.tfi)) * (
                ypl + self.dl + self.tfi) / 2) \
                / (Am_sup + (Aalma - self.tw * (ypl - self.tfi)))

            yi = (Am_inf * self.tfi / 2 + self.tw * (ypl - self.tfi) * (self.tfi + ypl) / 2) \
                / (Am_inf + self.tw * (ypl - self.tfi))

        Zx = A * (ys - yi) / 2
        Zy = (Am_inf * self.bfi + Aalma * self.tw + Am_sup * self.bfs) / 4

        # Áreas de cisalhamento (Awx e Awy)
        # ---------------------------------

        Awx = self.dl * self.tw
        Awy = self.tfi * self.bfi + self.tfs * self.bfs

        # Coordenadas do centro de cisalhamento em relação ao centro geométrico (xo e yo)
        # ---------------------------------------------------------------------------

        xo = 0

        h = self.d - self.tfi / 2 + self.tfs / 2
        ycc = (self.d - self.tfs/2) - h * I2y / (I1y + I2y)
        yo = ycc - ycg

        # Constante de empenamento (Cw)
        # -----------------------------

        C = (self.bfi * self.bfs) ** 3 / (self.bfi ** 3 + self.bfs ** 3)
        Cw = C * (self.tfi + self.tfs) * h ** 2 / 24

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
        self.h = self.d - 2 * self.tf
        self.dl = self.h - 2 * self.r

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

        self.esb_alma = self.hint / tw
        self.esb_mesa = self.bint / tf

        simetria = [True, True]

        super().__init__(**self.prop_geo(), material=material, simetria=simetria, tipo='CAIXAO')

        self.Wx = self.Wxs
        self.Wy = self.Wys

    def prop_geo(self):

        # Área (A)
        # ---------

        Amesa = self.b * self.tf
        Aalma = self.hint * self.tw

        A = 2 * Amesa + 2 * Aalma

        # Momentos de inércia e constante de torção (Ix, Iy e J)
        # ------------------------------------------------------

        Imx = self.b * self.tf ** 3 / 12
        Iax = self.tw * self.hint ** 3 / 12

        Imy = self.tf * self.b ** 3 / 12
        Iay = self.hint * self.tw ** 3 / 12

        dmy = self.h / 2 - self.tf / 2
        dax = self.b / 2 - self.tw / 2

        Ix = 2 * (Imx + Amesa * dmy ** 2 + Iax)
        Iy = 2 * (Iay + Aalma * dax ** 2 + Imy)

        # Módulos elásticos (Wxs, Wxi, Wys e Wyi)
        # --------------------------------------

        Wx = 2 * Ix / self.b
        Wy = 2 * Iy / self.b

        u = 2 * (self.b + self.h - self.tw - self.tf)
        K = (self.b - self.tw) * (self.h - self.tf) * (self.tf + self.tw) / u
        J = (self.tw ** 3 + self.tf ** 3) * u / 6 + 2 * K * (self.b - self.tw) * (self.h - self.tf)

        # Módulos plásticos (Zx e Zy)
        # -------------------------

        Zx = (self.b * self.h ** 2 - self.bint * self.hint ** 2) / 4
        Zy = (self.h * self.b ** 2 - self.hint * self.bint ** 2) / 4

        # Áreas de cisalhamento (Awy e Awx)
        # ---------------------------------

        Awy = 2 * self.hint * self.tw
        Awx = 2 * self.bint * self.tf

        # Centro de corte em relação ao centro geométrico (xo e yo)
        # ---------------------------------------------------------

        xo = 0
        yo = 0

        # Constante de empenamento (Cw)
        # -----------------------------

        Cw = 0

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wx, 'Wxi': Wx, 'Wyi': Wy, 'Wys': Wy,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}


class TuboRet(PerfilEstrutural):

    def __init__(self, nome, material):
        perfil = perfis_AISC[perfis_AISC['EDI_Std_Nomenclature.1'] == nome]

        self.h = float(perfil['Ht.1'])
        self.b = float(perfil['B.1'])
        self.esb_alma = float(perfil['h/tdes.1'])
        self.esb_mesa = float(perfil['b/tdes.1'])
        self.t = float(perfil['tdes.1'])
        self.tw = self.t
        self.tf = self.t

        self.bint = self.b - 3 * self.t
        self.hint = self.h - 3 * self.t

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

        self.Dint = self.D - 2 * self.t

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
