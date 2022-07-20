
from math import pi, sqrt
from util import PropGeo, NumPositivo
import pandas as pd
from collections import namedtuple
from types import MethodType

from secao import SecaoGenerica
from normas import NBR8800, AISC360
from material import Material, Aco

#  Importando o banco de dados de perfis
perfis_AISC = pd.read_excel('db-aisc-perfis.xlsx')
perfis_vallourec_ret = pd.read_excel('db-vallourec-perfis.xlsx')
perfis_vallourec_cir = pd.read_excel('db-vallourec-perfis.xlsx', 1)

db_perfis = pd.concat([perfis_AISC, perfis_vallourec_cir, perfis_vallourec_ret], sort=False)

NORMAS = {'NBR8800': NBR8800, 'AISC360': AISC360}

Und = {'mm': 1, 'cm': 10, 'm': 10e3}


class PerfilEstrutural(SecaoGenerica):
    """
    Esta classe define uma seção tranversal de barra de formato genérico
    de acordo com suas propriedades geométricas e seu material, e serve como
    base para a implementação das classes de perfis de formas específicos.

    Parameters
    ----------
    A: float
        área total da seção transversal
    Ix: float
        momento de inércia a flexão da seção em relação ao eixo X (horizontal)
        que passa pelo centroide da seção.
    Iy: float
        momento de inércia a flexão da seção em relação ao eixo Y (Vertical)
        que passa pelo centroide da seção.
    J: float
        constante de torção da seção em relação ao centróide da seção
    mat: `.material.Material`, list, dict, str
        material que compõe o perfil.
    Wx: float
        módulo elástico da seção em relação ao eixo X (horizontal)
    Wy: float
        módulo elástico da seção em relação ao eixo Y (Vertical)
    Zx: float,
        módulo plástico da seção em relação ao eixo X (horizontal)
    Zy: float
        módulo plástico da seção em relação ao eixo Y (vertical)
    Awx: float
        Área efetiva de cisalhamento em relação ao eixo X (horizontal)
    Awy: float
        Área efetiva de cisalhamento em relação ao eixo Y (vertical)
    xo: float, opcional
        coordenada x do centro de corte da seção trasnversal em relação ao
        centróide da seção
    yo: float, opcional
        coordenada y do centro de corte da seção trasnversal em relação ao
        centróide da seção
    Cw: float
        constante de empenamento da seção tranversal
    simetria: lista de bool
        indica se a seção apresenta eixos de simetria
    """

    _tipos_validos = None

    Wxs = PropGeo('Wxs', 3)
    Wxi = PropGeo('Wxi', 3)
    Wys = PropGeo('Wys', 3)
    Wyi = PropGeo('Wyi', 3)
    Zx  = PropGeo('Zx', 3)
    Zy  = PropGeo('Zy', 3)
    Awx = NumPositivo('Awx')
    Awy = NumPositivo('Awy')
    Cw  = PropGeo('Cw', 6)

    def __init__(self, A, Ix, Iy, J, Wxs, Wxi, Wys, Wyi, Zx, Zy, Awx, Awy, xo, yo, Cw, mat, simetria, norma=None,
                 tipo='GENERICO'):

        if isinstance(mat, list):
            mat = Aco(*mat)
        if isinstance(mat, dict):
            mat = Aco(**mat)

        super().__init__(A, Ix, Iy, J, mat)

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
        self.mat = mat
        self.simetria = simetria
        self.tipo = tipo

        self.raiz_E_fy = sqrt(self.mat.E / self.mat.fy)
        self.raiz_fy_E = sqrt(self.mat.fy / self.mat.E)

        if norma is not None:
            self.definir_norma(norma)


    @property
    def rx(self):
        """ Raio de giração em relação ao X da seção transversal """
        return sqrt(self.Ix / self.A)

    @property
    def ry(self):
        """ Raio de giração em relação ao Y da seção transversal """
        return sqrt(self.Iy / self.A)

    @property
    def ro(self) -> float:
        """ Raio de giração polar de inércia da seção em relação ao centro de cisalhamento """
        return sqrt((self.Iy + self.Ix) / self.A + self.xo ** 2 + self.yo ** 2)

    # Capacidade resistente do perfil sem consideração dos efeitos de instabilidade

    @property
    def Afy(self) -> float:
        """ Resitência axial ao escoamento"""
        return self.A * self.mat.fy

    @property
    def Afu(self) -> float:
        """ Resitência axial a ruptura"""
        return self.A * self.mat.fu

    @property
    def Vplx(self):
        """ Força cortante resistênte na direção X """
        return self.Awx * self.mat.fy * self.mat.cv

    @property
    def Vply(self):
        """ Força cortante resistênte na direção Y"""
        return self.Awy * self.mat.fy * self.mat.cv

    @property
    def Mplx(self):
        """ Momento de plastificação da seção em relação ao eixo X"""
        return self.Zx * self.mat.fy

    @property
    def Mply(self):
        """ Momento de plastificação da seção em relação ao eixo Y"""
        return self.Zy * self.mat.fy

    @property
    def Mrx(self):
        """ Momento de início de escoamento da seção em relação ao eixo X """
        return min(self.Wxs, self.Wxi) * self.mat.fy

    @property
    def Mry(self):
        """ Momento de início de escoamento da seção em relação ao eixo Y """
        return min(self.Wys, self.Wyi) * self.mat.fy

    @property
    def simetria_x(self):
        return self.simetria[0]

    @property
    def simetria_y(self):
        return self.simetria[1]

    @property
    def bissimetrico(self):
        return True if self.simetria_x and self.simetria_y else False

    # Métodos para definição da carga critica de flambagem de um elemento de barra de comprimentos
    # de flambagem klx, kly e klz

    def indice_esbeltez(self, Lx, Ly):
        """
        Retorna o indice de esbeltez de uma barra de comprimento destravado Lb
        formado pelo perfil em relação aos eixos X e Y, respectivamente.

        Parameter
        ---------
        Lb: 'float'
            comprimento destravado da barra

        Return
        ------

        """
        return Lx / self.rx, Ly / self.ry

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
        par_estabilidade = namedtuple('par_estabilidade', 'Ne fe Nex Ney Nez Nexz Neyz Me')

        Nex = pi ** 2 * self.EIx / (klx ** 2)
        Ney = pi ** 2 * self.EIy / (kly ** 2)
        Nez = (pi ** 2 * self.mat.E * self.Cw / (klz ** 2) + self.GJ) / self.ro ** 2
        Nexz = None
        Neyz = None
        Me = None

        if self.bissimetrico:
            Ne = min(Nex, Ney, Nez)
            Me = self.ro * sqrt(Ney * Nez) if self.Ix > self.Iy else self.ro * sqrt(Nex * Nez)

        elif not self.bissimetrico and self.simetria_y:

            Neyz = (Ney + Nez) / (2 * (1 - (self.yo / self.ro) ** 2)) * \
                   (1 - sqrt(1 - 4 * Ney * Nez * (1 - (self.yo / self.ro) ** 2) / (Ney + Nez) ** 2))

            Ne = min(Nex, Neyz)
            Me = self.ro * sqrt(Nex * Nez)

        elif not self.bissimetrico and self.simetria_x:

            Nexz = (Nex + Nez) / (2 * (1 - (self.xo / self.ro) ** 2)) * \
                   (1 - sqrt(1 - 4 * Nex * Nez * (1 - (self.xo / self.ro) ** 2) / (Nex + Nez) ** 2))

            Ne = min(Ney, Nexz)
            Me = self.ro * sqrt(Nex * Nez)

        else:

            a = 1 - (self.yo ** 2 + self.xo ** 2) / self.ro ** 2
            b = -Nex - Ney - Nez + (self.xo / self.ro) ** 2 * Ney + (self.yo / self.ro) ** 2 * Nex
            c = Ney * Nex + Nex * Nez
            d = - Nex * Ney * Nez

            det_0 = b ** 2 - 3 * a * c
            det_1 = 2 * b ** 3 - 9 * a * b * c + 27 * a ** 2 * d
            det = (det_1 ** 2 - 4 * det_0 ** 3) / (-27 * a ** 2)

            C = ((sqrt(det_1 ** 2 - 4 * det_0 ** 3) + det_1) / 2) ** (1 / 3)

            u = -1 + sqrt(-3) / 2

            raizes = list()

            for n in range(1, 4):
                raizes.append((- (b + (u ** n * C) + det_0) / (u ** n * C)) / (3 * a))

            Ne = min(raizes)

        fe = Ne / self.A

        return par_estabilidade(Ne, fe, Nex, Ney, Nez, Nexz, Neyz, Me)

    def prop(self):
        """ Retorna as propriedades geométricas do perfil """
        propriedades_perfil = namedtuple('propriedades_perfil', 'A Ix Iy J Wxs Wxi Wys Wyi Zx Zy Awx Awy xo yo Cw')
        return propriedades_perfil(self.A, self.Ix, self.Iy, self.J, self.Wxs, self.Wxi, self.Wys, self.Wyi, self.Zx,
                                   self.Zy, self.Awx, self.Awy, self.xo, self.yo, self.Cw)

    def definir_norma(self, norma):
        """ Método responsável por definir os métodos de verificação de acordo com a norma requerida"""

        if norma not in NORMAS.keys():
            raise ValueError('A norma {} não está implementada'.format(norma))

        else:
            metodos = ('Ntrd_brt', 'Ncrd', 'Vrdx', 'Vrdy', 'Mrdx', 'Mrdy', 'Trd')
            for metodo in metodos:
                self.__setattr__(metodo, MethodType(getattr(NORMAS[norma], metodo), self))

    def _validar_nome(self, nome):
        perfil = db_perfis[db_perfis['Nomes'] == nome]

        str_tipo = perfil.Tipo.values[0]

        if str_tipo not in self._tipos_validos:
            raise ValueError('{} não é um nome válido para o perfil do tipo {}'.format(nome, self.tipo))

        else:
            self._dados_perfil = perfil


class PerfilI(PerfilEstrutural):
    """
    Está classe define um perfil I soldado mono ou duplamente simétrico. O eixo Y é tomado como eixo de vertical,
    e o eixo X como o eixo horizontal.
    Perfis I soldados monossimétricos apresentam sua simetria em relação ao eixo Y.

    Parameter
    ---------
    d: float
       Altura total do perfil
    bfs: float
       Largura da mesa superior
    bfi: float
        Largura da mesa inferior
    tw: float
        Espessura da alma
    tfs: float
        Espessura da mesa superior
    tfi: float
        Espessura da mesa inferior
    mat: `.material.Material`, list, dict, str
        material que compõe o perfil.
    norma: str, opcional
        Nome na norma na qual se deseja fazer as verificações de capacidade resistente de perfil
        O nome da norma dever ser entre os apresentados abaixo

            -'NBR8800' :ref: `(veja mais) <metalpy.norma._nbr8800>`
            -'AISC360' :ref: `(veja mais) <metalpy.norma._aisc360>`

    Examples
    --------
    Exercicio 6.5.1 do livro **Estrutura de aço - dimensionamento prático** `Pfeil, Michèle; Pfeil, Walter`, pag.197:

    >>>from perfis import PerfilI
    >>>from material import Aco

    >>>#Definindo o aço do tipo MR250 com as propriedades em kN/cm²
    >>>MR250 = Aco(20000, 0.3, 25, 30, 0.6)

    >>>#Dados do perfil VS 500 X 86 em cm
    >>>d = 50
    >>>bf = 25
    >>>tw = 0.63
    >>>tf = 1.6

    >>>VS500X86 = PerfilI(d, bf, bf, tw, tf, tf, MR250, 'NBR8800')

    Attribute
    ----------
    d: float
       Altura total do perfil
    dl: float
        Distância interna entre as mesas
    h: float
        Largura plana da alma
    bfs: float
       Largura da mesa superior
    bfi: float
        Largura da mesa inferior
    tw: float
        Espessura da alma
    tfs: float
        Espessura da mesa superior
    tfi: float
        Espessura da mesa inferior
    esb_alma: float
        Esbeltez da alma
    esb_mesa_s: float
        Esbeltez da mesa superior
    esb_mesa_i: float
        Esbeltez da mesa inferior
    Iys: float
        Momento de inércia da mesa superior em relação ao eixo y
    Iyi: float
        Momento de inércia da mesa inferior em relação ao eixo y
    hcg: float
        Altura do centro geométrico
    hpl: float
        Altura da linha neutra na plastificação da seção
    A: float
        Área do perfil
    Ix: float
        Momento de inércia do perfil em relação ao eixo X
    Iy: float
        Momento de inércia do perfil em relação ao eixo Y
    Wxs: float
        Módulo elástico superior do perfil em relação ao eixo X, para perfis monossimétricos
    Wxi: float
        Módulo elástico inferior do perfil em relação ao eixo X, para perfis monossimétricos
    Wx: float
        Módulo elástico em relação ao eixo X, para perfis bissimétricos.
    Wy: float
        Módulo elástico do perfil em relação ao eixo Y
    Zx: float
        Módulo plástico do perfil em relação ao eixo X
    Zy: float
        Módulo plástico do perfil em relação ao eixo Y
    Awx: float
        Área efetiva de cisalhamento paralela ao eixo X
    Awy: float
        Área efetiva de cisalhamento paralela ao eixo Y
    xo: float
        coordenada X do centro de cisalhamento em relação ao centróide da seção
    yo: float
        coordenada Y do centro de cisalhamento em relação ao centróide da seção
    Cw: float
        Costante de empenamento
    mat: objeto `.material.Material`
        Material que compõe o perfil
    simetria: list de bool
        Indica se existe simetria em relação aos eixos X e Y.

        `True` indica a existência de simetria e `False` indica assimetria,
        o primeiro termo da lista indica a simetria em relação ao eixo X (eixo horizontal),
        e o segundo indica a simetria em relação ao eixo Y (eixo vertical)
    tipo: str, default: 'I SOLDADO'
        tipo de perfil
    """

    tipo = 'I SOLDADO'

    d = NumPositivo('d')
    bfs = NumPositivo('bfs')
    bfi = NumPositivo('bfi')
    tw = NumPositivo('tw')
    tfs = NumPositivo('tfs')
    tfi = NumPositivo('tfi')

    def __init__(self, d, bfs, bfi, tw, tfs, tfi, mat, norma=None):

        self.d = d
        self.dl = d - tfi - tfs
        self.h = self.dl
        self.tw = tw
        self.bfs = bfs
        self.bfi = bfi
        self.tfs = tfs
        self.tfi = tfi
        self.esb_mesa_s = bfs / (2 * tfs)
        self.esb_mesa_i = bfi / (2 * tfi)
        self.esb_alma = self.dl / tw
        self.Iys = None
        self.Iyi = None
        self.hcg = None
        self.hpl = None

        simetria = [True, True]
        simetria[0] = False if bfs != bfi or tfi != tfs else True

        prop_geo = self._prop_geo()
        super().__init__(**prop_geo, mat=mat, simetria=simetria, norma=norma, tipo=self.tipo)

        self.Wy = self.Wys

        if self.bissimetrico:
            self.Wx = self.Wxi
            self.esb_mesa = self.esb_mesa_i

    def _prop_geo(self):
        """
        Método que calcula as propriedades geométricas de um perfil I soldado bissimetrico ou
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

        ycg = (Am_inf * self.tfi / 2 + Aalma * (self.tfi + self.dl / 2)
               + Am_sup * (self.d - self.tfs / 2)) / A

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

            # Altura da linha neutra na plastificação da seção
            ypl = self.d - (Aalma + Am_inf + Am_sup) / (2 * self.bfs)

            ys = (ypl + self.d) / 2

            ycg_mi = self.tfi / 2

            yi = (Am_inf * ycg_mi + (self.tfi + self.dl / 2) + self.bfs * (self.tfs - (self.d - ypl)) *
                  (ypl + self.tfi + self.dl / 2)) / (Am_inf + Aalma + self.bfs * (self.tfs - (self.d - ypl)))

        elif Am_inf > Aalma + Am_sup:

            # Altura da linha neutra na plastificação da seção
            ypl = (Aalma + Am_sup + self.bfi * self.tfi) / (2 * self.bfi)

            ys = (Am_sup * (self.d - self.tfs / 2) + Aalma * (self.tfi + self.dl / 2) + self.bfi * (
                    ypl + self.tfi) / 2) / \
                 (Am_sup + Aalma + self.bfi * (self.tfi - ypl))

            yi = (self.tfi - ypl) / 2

        else:

            # Altura da linha neutra na plastificação da seção
            ypl = self.dl + self.tfi - (Am_inf - Am_sup + self.tw * self.dl) / (2 * self.tw)

            Aalma_sup = (Aalma - self.tw * (ypl - self.tfi))
            Aalma_inf = Aalma - Aalma_sup
            ycg_mesa_sup = (self.d - self.tfs / 2)
            ycg_mesa_inf = self.tfi / 2
            ycg_alma_sup = (ypl + self.dl + self.tfi) / 2
            ycg_alma_inf = (self.tfi + ypl) / 2

            ys = (Am_sup * ycg_mesa_sup + Aalma_sup * ycg_alma_sup) / (Am_sup + Aalma_sup)
            yi = (Am_inf * ycg_mesa_inf + Aalma_inf * ycg_alma_inf) / (Am_inf + Aalma_inf)

        Zx = A * (ys - yi) / 2
        Zy = (Am_inf * self.bfi + Aalma * self.tw + Am_sup * self.bfs) / 4

        # Áreas de cisalhamento (Awx e Awy)
        # ---------------------------------

        Awx = self.dl * self.tw
        Awy = self.tfi * self.bfi + self.tfs * self.bfs

        # Coordenadas do centro de cisalhamento em relação ao centro geométrico (xo e yo)
        # ---------------------------------------------------------------------------

        xo = 0

        h = self.d - self.tfi / 2 - self.tfs / 2
        ycc = (self.d - self.tfs / 2) - h * I2y / (I1y + I2y)
        yo = ycc - ycg

        # Constante de empenamento (Cw)
        # -----------------------------

        C = (self.bfi * self.bfs) ** 3 / (self.bfi ** 3 + self.bfs ** 3)
        Cw = C * (self.tfi + self.tfs) * (h ** 2) / 24

        self.Iys = I1y + Am_sup * (self.dl + self.tfs) ** 2 / 4
        self.Iyi = I2y + Am_inf * (self.dl + self.tfi) ** 2 / 4
        self.hpl = ypl
        self.hcg = ycg

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wxs, 'Wxi': Wxi, 'Wyi': Wyi, 'Wys': Wys,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}


class PerfilILam(PerfilEstrutural):
    """
    Está classe define um perfil do I laminado, podendo ser dos tipos W, S, H e HP.

    Parameter
    ---------
    nome: str
        Nome do perfil, definido da seguinte forma '(TIPO)(ALTURA DO PERFIL EM mm)X(PESO DO PERFIL EM kg/m)'
        sendo 'tipo' um dos seguintes valores: 'W', 'HP' e 'S'

        Exemplos:
            - 'W530X85'
            - 'HP200X53'
            - 'S380X81.4'
    mat: `.material.Material`, list, dict, str
        material que compõe o perfil.
    und: str, default:'mm'
        unidade de medida das propriedades do perfil.
        O nome das unidades deve ser um entre os apresentados abaixo:
            - 'mm'
            - 'cm'
            - 'm'
        Por definição os valores da base dados estão em mm.
    norma: str, opcional
        Nome na norma na qual se deseja fazer as verificações de capacidade resistente de perfil
        O nome da norma dever ser entre os apresentados abaixo

            - 'NBR8800' :ref: `(veja mais) <metalpy.norma._nbr8800>`
            - 'AISC360' :ref: `(veja mais) <metalpy.norma._aisc360>`

    Examples
    --------
    Exercicio 6.5.1 do livro **Estrutura de aço - dimensionamento prático** `Pfeil, Michèle; Pfeil, Walter`, pag.197:

    >>> from perfis import PerfilILam
    >>> from material import Aco

    >>> #Definindo o aço do tipo MR250 com as propriedades em kN/cm²
    >>>  MR250 = Aco(20000, 0.3, 25, 30, 0.6)

    >>> #Definindo um perfil I laminado do tipo W com as propriedades em cm
    >>>  W530X85 = PerfilILam('W530X85', MR250, und='cm', norma='NBR8800')

    Attribute
    ----------
    d: float
       Altura total do perfil
    dl: float
        Distância interna as mesas dos perfis
    h: float
        Largura plana da alma
    bf: float
       Largura da mesa
    tw: float
        Espessura da alma
    tf: float
        Espessura da mesa
    r: float
        Raio de concordância
    esb_mesa: float
        Esbeltez da mesa
    esb_alma: float
        Esbeltez da alma
    A: float
        Área do perfil
    Ix: float
        Momento de inércia do perfil em relação ao eixo X
    Iy: float
        Momento de inércia do perfil em relação ao eixo Y
    Wx: float
        Módulo elástico em relação ao eixo X, para perfis bissimétricos.
    Wy: float
        Módulo elástico do perfil em relação ao eixo Y
    Zx: float
        Módulo plástico do perfil em relação ao eixo X
    Zy: float
        Módulo plástico do perfil em relação ao eixo Y
    Awx: float
        Área efetiva de cisalhamento paralela ao eixo X
    Awy: float
        Área efetiva de cisalhamento paralela ao eixo Y
    xo: float
        coordenada X do centro de cisalhamento em relação ao centróide da seção
    yo: float
        coordenada Y do centro de cisalhamento em relação ao centróide da seção
    Cw: float
        Costante de empenamento
    mat: objeto `.material.Material`
        Material que compõe o perfil
    simetria: list de bool, default: [True, True]
        Indica se existe simetria em relação aos eixos X e Y.
    tipo: str, default: 'I LAMINADO'
        tipo de perfil
    """

    _tipos_validos = ('W', 'HP', 'S')
    tipo = 'I LAMINADO'

    d = PropGeo('d', 1)
    bf = PropGeo('bf', 1)
    tf = PropGeo('tf', 1)
    tw = PropGeo('tw', 1)
    kdes = PropGeo('kdes', 1)

    def __init__(self, nome, mat, und='mm', norma=None):

        self._validar_nome(nome)

        self.und = und

        self.d = float(self._dados_perfil['d'])
        self.bf = float(self._dados_perfil['bf'])
        self.tf = float(self._dados_perfil['tf'])
        self.tw = float(self._dados_perfil['tw'])
        self.kdes = float(self._dados_perfil['kdes'])

        self.r = self.kdes - self.tf
        self.h = self.d - 2 * self.tf
        self.dl = self.h - 2 * self.r

        self.esb_alma = float(self._dados_perfil['h/tw'])
        self.esb_mesa = float(self._dados_perfil['bf/2tf'])

        simetria = [True, True]

        super().__init__(**self._prop_geo(), mat=mat, simetria=simetria, norma=norma, tipo=self.tipo)

        self.Wx = self.Wxs
        self.Wy = self.Wys

    def _prop_geo(self):
        """ Método que obtém as propriedades do perfil do banco de dados """
        A = float(self._dados_perfil['A'])
        Ix = float(self._dados_perfil['Ix']) * 1E6
        Iy = float(self._dados_perfil['Iy']) * 1E6
        J = float(self._dados_perfil['J']) * 1E3
        Wx = float(self._dados_perfil['Wx']) * 1E3
        Wy = float(self._dados_perfil['Wy']) * 1E3
        Zx = float(self._dados_perfil['Zx']) * 1E3
        Zy = float(self._dados_perfil['Zy']) * 1E3
        Cw = float(self._dados_perfil['Cw']) * 1E9

        Awy = self.d * self.tw
        Awx = 2 * self.bf * self.tf

        xo = 0
        yo = 0

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wx, 'Wxi': Wx, 'Wyi': Wy, 'Wys': Wy,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}


class Caixao(PerfilEstrutural):
    """
    Está classe define um perfil soldado do tipo seção caixão retangular duplamente simétrico,
    tendo como eixo X o eixo paralelo as mesas e eixo Y o eixo perpendicular as mesas.

    Parameter
    ---------
    h: float
        Altura total do perfil
    b: float
        Largura total do perfil
    tw: float
        Espessura das almas
    tf: float
        Espessura das mesas
    mat: `.material.Material`, list, dict, str
        material que compõe o perfil.
    norma: str, opcional
        Nome na norma na qual se deseja fazer as verificações de capacidade resistente de perfil
        O nome da norma dever ser entre os apresentados abaixo

            - 'NBR8800' :ref: `(veja mais) <metalpy.norma._nbr8800>`
            - 'AISC360' :ref: `(veja mais) <metalpy.norma._aisc360>`
    Examples
    --------
    >>> from perfis import Caixao
    >>> from material import Aco

    >>> #Definindo o aço do tipo MR250 com as propriedades em kN/cm²
    >>> MR250 = Aco(20000, 0.3, 25, 30, 0.6)

    >>> #Definindo uma seção caixão com 20cm de altura, 20cm de largura e espessura de 3.5cm
    >>> C200x15 = Caixao(20, 20, 3.5, 3.5)

    Attribute
    ---------
    h: float
        Altura total do perfil
    b: float
        Largura total do perfil
    tw: float
        Espessura das almas
    tf: float
        Espessura das mesas
    hint: float
        Distância interna entre as mesas do perfil
    bint: float
        Distância interna entre as almas do perfil
    esb_alma: float
        Esbeltez das almas
    esb_mesa: float
        Esbeltez das mesas
    A: float
        Área do perfil
    Ix: float
        Momento de inércia do perfil em relação ao eixo X
    Iy: float
        Momento de inércia do perfil em relação ao eixo Y
    Wx: float
        Módulo elástico em relação ao eixo X.
    Wy: float
        Módulo elástico do perfil em relação ao eixo Y
    Zx: float
        Módulo plástico do perfil em relação ao eixo X
    Zy: float
        Módulo plástico do perfil em relação ao eixo Y
    Awx: float
        Área efetiva de cisalhamento paralela ao eixo X
    Awy: float
        Área efetiva de cisalhamento paralela ao eixo Y
    xo: float
        coordenada X do centro de cisalhamento em relação ao centróide da seção
    yo: float
        coordenada Y do centro de cisalhamento em relação ao centróide da seção
    Cw: float
        Costante de empenamento
    mat: objeto `.material.Material`
        Material que compõe o perfil
    simetria: list de bool, default: [True, True]
        Indica se existe simetria em relação aos eixos X e Y.
    tipo: str, default: 'CAIXAO'
        tipo de perfil
    """

    tipo = 'CAIXAO'

    h = NumPositivo('h')
    b = NumPositivo('b')
    tw = NumPositivo('tw')
    tf = NumPositivo('tf')

    def __init__(self, h, b, tw, tf, mat, norma=None):
        self.h = h
        self.b = b
        self.tw = tw
        self.tf = tf
        self.hint = self.h - 2 * tf
        self.bint = self.b - 2 * tw

        self.esb_alma = self.hint / tw
        self.esb_mesa = self.bint / tf

        self.Wt = (self.b - self.tf) * (self.h - self.tw) * (self.tf + self.tw)

        simetria = [True, True]

        super().__init__(**self._prop_geo(), mat=mat, simetria=simetria, norma=norma, tipo='CAIXAO')

        self.Wx = self.Wxs
        self.Wy = self.Wys

    def _prop_geo(self):
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

    """
    Esta classe define um perfil laminado do tipo tubo retangular

    Parameter
    ---------
    nome: str
        Nome do perfil, definido da seguinte forma:
         '(TIPO)(ALTURA TOTAL DO PERFIL mm)X(LARGURA TOTAL EM mm)X(ESPESSURA EM mm)'

        sendo 'tipo' um dos seguintes valores: 'TQ', 'TR' e 'HSS', 'TQ' e 'TR' são perfis do
        catalogo da Vallourec e 'HSS' são perfis padronizados pela AISC.

        Exemplos:
            - 'TQ110X110X12.5'
            - 'TR360X210X11'
            - 'HSS508X508X22.2'
    mat: `.material.Material`, list, dict, str
        material que compõe o perfil.
    und: str, default:'mm'
        unidade de medida das propriedades do perfil.
        O nome das unidades deve ser um entre os apresentados abaixo:
            - 'mm'
            - 'cm'
            - 'm'
        Por definição os valores da base dados estão em mm.
    norma: str, opcional
        Nome na norma na qual se deseja fazer as verificações de capacidade resistente de perfil
        O nome da norma dever ser entre os apresentados abaixo

            - 'NBR8800' :ref: `(veja mais) <metalpy.norma._nbr8800>`
            - 'AISC360' :ref: `(veja mais) <metalpy.norma._aisc360>`
    Examples
    --------

    Attribute
    ---------
    h: float
        Altura total do perfil
    b: float
        Largura total do perfil
    tw: float
        Espessura das almas
    tf: float
        Espessura das mesas
    hint: float
        Altura plana da alma do perfil
    bint: float
        Largura plana das mesas do perfil
    esb_alma: float
        Esbeltez das almas
    esb_mesa: float
        Esbeltez das mesas
    A: float
        Área do perfil
    Ix: float
        Momento de inércia do perfil em relação ao eixo X
    Iy: float
        Momento de inércia do perfil em relação ao eixo Y
    Wx: float
        Módulo elástico em relação ao eixo X.
    Wy: float
        Módulo elástico do perfil em relação ao eixo Y
    Wt: float
        Módulo de torção
    Zx: float
        Módulo plástico do perfil em relação ao eixo X
    Zy: float
        Módulo plástico do perfil em relação ao eixo Y
    Awx: float
        Área efetiva de cisalhamento paralela ao eixo X
    Awy: float
        Área efetiva de cisalhamento paralela ao eixo Y
    xo: float
        coordenada X do centro de cisalhamento em relação ao centróide da seção
    yo: float
        coordenada Y do centro de cisalhamento em relação ao centróide da seção
    Cw: float
        Costante de empenamento
    mat: objeto `.material.Material`
        Material que compõe o perfil
    simetria: list de bool, default: [True, True]
        Indica se existe simetria em relação aos eixos X e Y.
    tipo: str, default: 'TUBO RET'
        tipo de perfil
    """

    tipo = 'TUBO RET'

    _tipos_validos = ('TC', 'TR', 'HSS')

    h = PropGeo('h', 1)
    b = PropGeo('b', 1)
    t = PropGeo('t', 1)
    Wt = PropGeo('Wt', 3)

    def __init__(self, nome, mat, und='mm', norma=None):

        self._validar_nome(nome)

        self.und = und

        self.h = float(self._dados_perfil['Ht'])
        self.b = float(self._dados_perfil['B'])
        self.esb_alma = float(self._dados_perfil['h/tdes'])
        self.esb_mesa = float(self._dados_perfil['b/tdes'])
        self.t = float(self._dados_perfil['tdes'])
        self.tw = self.t
        self.tf = self.t
        self.Wt = float(self._dados_perfil['Wt']) * 1E3

        self.bint = self.esb_mesa * self.t
        self.hint = self.esb_alma * self.t

        simetria = [True, True]

        super().__init__(**self.prop_geo(), mat=mat, simetria=simetria, norma=norma,  tipo=self.tipo)

        self.Wx = self.Wxs
        self.Wy = self.Wys

    def prop_geo(self):

        A = float(self._dados_perfil['A'])
        Ix = float(self._dados_perfil['Ix']) * 1E6
        Iy = float(self._dados_perfil['Iy']) * 1E6
        J = float(self._dados_perfil['J']) * 1E3
        Wx = float(self._dados_perfil['Wx']) * 1E3
        Wy = float(self._dados_perfil['Wy']) * 1E3
        Zx = float(self._dados_perfil['Zx']) * 1E3
        Zy = float(self._dados_perfil['Zy']) * 1E3

        Awy = 2 * self.hint * self.t
        Awx = 2 * self.bint * self.t

        xo = 0
        yo = 0

        Cw = 0

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wx, 'Wxi': Wx, 'Wyi': Wy, 'Wys': Wy,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}


class TuboCir(PerfilEstrutural):

    """
    Está classe define um perfil tubo circular

    Parameter
    ---------
    nome: str
        Nome do perfil, definido da seguinte forma:
         '(TIPO)(DIÂMETRO EXTERNO DO PERFIL mm)X(ESPESSURA EM mm)', para perfis dos tipo 'TC' e 'HSS',
         para perfis do tipo 'Pipe' consultar nomes em ...

        Exemplos:
            - 'HSS508X12.7'
            - 'TC88.9X8'
            - 'Pipe200XS'
    mat: `.material.Material`, list, dict, str
        material que compõe o perfil.
    und: str, default:'mm'
        unidade de medida das propriedades do perfil.
        O nome das unidades deve ser um entre os apresentados abaixo:
            - 'mm'
            - 'cm'
            - 'm'
        Por definição os valores da base dados estão em mm.
    norma: str, opcional
        Nome na norma na qual se deseja fazer as verificações de capacidade resistente de perfil
        O nome da norma dever ser entre os apresentados abaixo

            - 'NBR8800' :ref: `(veja mais) <metalpy.norma._nbr8800>`
            - 'AISC360' :ref: `(veja mais) <metalpy.norma._aisc360>`

    Examples
    --------
    Exercicio 3.7.7 adaptado do livro ***Projeto de estruturas de Edificações com perfis Tubulates de Aço*** pag.190

    >>>from perfis import TuboCir
    >>>from material import Aco

    >>>#Definição do Aço com as propriedades em kN/cm²
    >>>A572 = Aco(20000, 0.3, 34.5, 45.0, 0.6)
    >>>L = 800 #cm

    >>>#Definindo um perfil TC323.8X6.4 (perfil catalogado mais próximo do apresentado no exercicio)
    >>>TC320X6 = TuboCir('TC323.8X6.4', A572, 'cm', 'NBR8800')

    >>>#Obtendo o momento resistênte de perfil
    >>> TC320X6.Mrdx()

    >>>#Obtendo a força cortante resistente
    >>>TC320X6.Vrdy(400), TC320X6.Awx

    Attribute
    ---------
    D: float
        Diâmetro externo do perfil
    t: float
        Espessura do perfil
    Dint: float
        Diâmetro interno do perfil
    esb: float
        Esbeltez do perfil
    A: float
        Área do perfil
    Ix: float
        Momento de inércia do perfil em relação ao eixo X
    Iy: float
        Momento de inércia do perfil em relação ao eixo Y
    Wx: float
        Módulo elástico em relação ao eixo X.
    Wy: float
        Módulo elástico do perfil em relação ao eixo Y
    Wt: float
        Módulo de torção
    Zx: float
        Módulo plástico do perfil em relação ao eixo X
    Zy: float
        Módulo plástico do perfil em relação ao eixo Y
    Awx: float
        Área efetiva de cisalhamento paralela ao eixo X
    Awy: float
        Área efetiva de cisalhamento paralela ao eixo Y
    xo: float
        coordenada X do centro de cisalhamento em relação ao centróide da seção
    yo: float
        coordenada Y do centro de cisalhamento em relação ao centróide da seção
    Cw: float
        Costante de empenamento
    mat: objeto `.material.Material`
        Material que compõe o perfil
    simetria: list de bool, default: [True, True]
        Indica se existe simetria em relação aos eixos X e Y.
    tipo: str, default: 'TUBO CIR'
        tipo de perfil
    """

    tipo = 'TUBO CIR'

    _tipos_validos = ('TC', 'HSS', 'PIPE')

    D = PropGeo('D', 1)
    t = PropGeo('t', 1)
    Awx = PropGeo('Awx', 2)
    Awy = PropGeo('Awy', 2)
    Wt = PropGeo('Wt', 3)

    esb = NumPositivo('esb')

    def __init__(self, nome, mat, und='mm', norma=None):
        self._validar_nome(nome)

        self.und = und

        self.D = float(self._dados_perfil['D'])
        self.t = float(self._dados_perfil['tdes'])
        self.esb = float(self._dados_perfil['D/t'])

        try:
            self.Wt = float(self._dados_perfil['Wt']) * 1E3
        except ValueError:
            self.Wt = 0

        self.Dint = self.D - 2 * self.t

        simetria = [True, True]

        super().__init__(**self.prop_geo(), mat=mat, simetria=simetria, norma=norma, tipo=self.tipo)

        self.W = self.Wxs

    def prop_geo(self):

        A = float(self._dados_perfil['A'])
        Ix = float(self._dados_perfil['Ix']) * 1E6
        Iy = float(self._dados_perfil['Iy']) * 1E6
        J = float(self._dados_perfil['J']) * 1E3
        Wx = float(self._dados_perfil['Wx']) * 1E3
        Wy = float(self._dados_perfil['Wy']) * 1E3
        Zx = float(self._dados_perfil['Zx']) * 1E3
        Zy = float(self._dados_perfil['Zy']) * 1E3

        Awy = 0.5 * A
        Awx = 0.5 * A

        xo = 0
        yo = 0

        Cw = 0

        return {'A': A, 'Ix': Ix, 'Iy': Iy, 'J': J, 'Wxs': Wx, 'Wxi': Wx, 'Wyi': Wy, 'Wys': Wy,
                'Zx': Zx, 'Zy': Zy, 'Awx': Awx, 'Awy': Awy, 'xo': xo, 'yo': yo, 'Cw': Cw}