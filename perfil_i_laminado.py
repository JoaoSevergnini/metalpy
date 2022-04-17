from math import sqrt, pi
from perfil_de_aco import PerfilDeAço
from material import Material
import pandas as pd


class PerfilILaminado(PerfilDeAço):
    """
    Esta classe define um perfil de aço do tipo I laminado (W, H, HP).

    Parameter
    ----------
    nome: 'str'
        nome da seção transversal de acordo com a tabela da IASC

    material: colocar o nome do aço ou uma lista com os parametros ou um objeto da propria classe material

    """

    perfis_W = pd.read_excel('aisc-shapes-database-v15.0.xlsx', 1).iloc[1:283, 84:]
    perfis_Hp = db2 = pd.read_excel('aisc-shapes-database-v15.0.xlsx', 1).iloc[329:351, 84:]
    perfis = pd.concat([perfis_W, perfis_Hp])

    def __init__(self, nome, material):

        perfil = self.perfis[self.perfis['EDI_Std_Nomenclature.1'] == nome]

        self.ht = float(perfil['d.1'])
        self.bf = float(perfil['bf.1'])
        self.tf = float(perfil['tf.1'])
        self.tw = float(perfil['tw.1'])
        self.r = float(perfil['kdes.1']) - self.tf

        A = float(perfil['A.1'])
        Ix = float(perfil['Ix.1']) * 1E6
        Iy = float(perfil['Iy.1']) * 1E6
        J = float(perfil['J.1']) * 1E3
        Wx = float(perfil['Sx.1']) * 1E3
        Wy = float(perfil['Sy.1']) * 1E3
        Zx = float(perfil['Zx.1']) * 1E3
        Zy = float(perfil['Zy.1']) * 1E3
        Cw = float(perfil['Cw.1']) * 1E9

        xo = 0
        yo = 0

        simetria = [True, True]

        super().__init__(A, Ix, Iy, J, Wx, Wy, Zx, Zy, xo, yo, Cw, material, simetria)

        self.esb_alma = self.h / self.tw
        self.esb_mesa = self.bf / (2 * self.tf)

    @property
    def hw(self):
        return self.ht - 2 * self.tf

    @property
    def h(self):
        return self.ht - 2 * self.tf - 2 * self.r

    # -------------------------------------------------------------------------------------
    # --------------------------Verificações de resistência--------------------------------
    # -------------------------------------------------------------------------------------

    # ----------------------------------NBR 8800-------------------------------------------
    # COMPRESSÃO
    # -----------

    def par_esbeltez_limites_AL_Ncrd(self):
        return 0.56 * self.raiz_E_fy, 1.03 * self.raiz_E_fy

    def fator_Qs(self):

        # elp = esbeltez limite para plastificação
        # elr = esbeltez limite para início de escoamento
        elp, elr = self.par_esbeltez_limites_AL_Ncrd()

        if elp > self.esb_mesa:
            return 1
        if elp < self.esb_mesa <= elr:
            return 1.415 - 0.74 * self.esb_alma * sqrt(self.material.fy / self.material.E)
        elif self.esb_mesa > elr:
            return 0.69 * self.material.e / (self.material.fy * self.esb_mesa ** 2)

    def fator_Qa(self, frc):

        tensao = self.material.fy * frc
        ca = 0.34

        bef = 1.92 * self.tw * sqrt(self.material.E / tensao) * \
              (1 - ca / self.esb_alma * sqrt(self.material.E / tensao))

        bef = bef if bef < self.h else self.h

        Aef = self.A - (self.h - bef) * self.tw

        return Aef / self.A

    # CORTANTE
    # -----------

    @property
    def Awx(self):
        return 2 * self.bf * self.tf

    @property
    def Awy(self):
        return self.ht * self.tw

    # CORTANTE EM X
    def kv_Vrdx(self, a=None):
        return 1.2

    # CORTANTE EM Y
    def kv_Vrdy(self, a=None):
        if a is None or a / self.h > 3 or a / self.h > (260 / (self.h / self.tw)) ** 2:
            return 5
        else:
            return 5 + 5 / (a / self.h) ** 2

    # MOMENTO EM X
    # ------------
    # Estado Limite FLT

    def par_esbeltez_limite_Mrdx_FLT(self):

        # parâmetro de esbeltez limite de plastificação (elp)
        elp = 1.76 * self.raiz_E_fy

        # parâmetro de esbeltez limite de inicio de escoamento (elr)
        beta_1 = (0.7 * self.material.fy * self.Wx) / (self.material.E * self.J)

        elr = (1.38 * sqrt(self.Iy * self.J)) / (self.ry * self.J * beta_1) * \
              sqrt(1 + sqrt(1 + 27 * self.Cw * beta_1 ** 2 / self.Iy))

        return elp, elr

    def Mrx_FLT(self):
        return 0.7 * self.material.fy * self.Wx

    def Mcrx_FLT(self, Cb, Lb):

        Mcr = (Cb * pi ** 2 * self.material.E * self.Iy / Lb ** 2) \
              * sqrt(self.Cw / self.Iy * (1 + 0.039 * self.J * Lb ** 2 / self.Cw))
        return Mcr

    def Mnx_FLT(self, Cb, Lb):

        esbeltez = self.indice_esbeltez_X(Lb)
        elp, elr = self.par_esbeltez_limite_Mrdx_FLT()

        if esbeltez < elp:
            return self.Mplx
        if elp < esbeltez < elr:
            return Cb * (self.Mplx - (self.Mplx - self.Mrx_FLT()) * (esbeltez - elp) / (elr - elp))
        elif esbeltez > elp:
            return self.Mcrx_FLT(Cb, Lb)

    # Estado Limite FLM

    def par_esbeltez_limite_Mrdx_FLM(self):
        return 0.38 * self.raiz_E_fy, 0.83 * sqrt(1 / 0.7) * self.raiz_E_fy

    def Mrx_FLM(self):
        return 0.7 * self.material.fy * self.Wx

    def Mcrx_FLM(self):
        return 0.69 * self.material.E * self.Wx / self.esb_mesa ** 2

    def Mnx_FLM(self):

        elp, elr = self.par_esbeltez_limite_Mrdx_FLM()

        if self.esb_mesa < elp:
            return self.Mplx
        if elp < self.esb_mesa < elr:
            return self.Mplx - (self.Mplx - self.Mrx_FLM()) * (self.esb_mesa - elp) / (elr - elp)
        elif self.esb_mesa > elr:
            return self.Mcrx_FLM()

    # Estado Limite FLA

    def par_esbeltez_limite_Mrdx_FLA(self):
        return 3.76 * self.raiz_E_fy, 5.7 * self.raiz_E_fy

    def Mrx_FLA(self):
        return self.material.fy * self.Wx

    def Mnx_FLA(self):

        elp, elr = self.par_esbeltez_limite_Mrdx_FLA()

        if self.esb_alma < elp:
            return self.Mplx
        elif elp < self.esb_alma < elr:
            return self.Mplx - (self.Mplx - self.Mrx_FLA()) * (self.esb_alma - elp) / (elr - elp)

    # MOMENTO EM Y
    # ------------
    # Estado Limite FLT
    def Mny_FLT(self, Cb=None, Lb=None):
        return self.Mply

    # Estado Limite FLM
    def par_esbeltez_limite_Mrdy_FLM(self):
        return 0.38 * self.raiz_E_fy, 0.83 * sqrt(1 / 0.7) * self.raiz_E_fy

    def Mry_FLM(self):
        return 0.70 * self.material.fy * self.Wy

    def Mcry_FLM(self):
        return 0.69 * self.material.E / self.esb_mesa ** 2 * self.Wy

    def Mny_FLM(self):

        elp, elr = self.par_esbeltez_limite_Mrdy_FLM()

        if self.esb_mesa < elp:
            return self.Mply
        elif elp < self.esb_mesa < elr:
            return self.Mply - (self.Mply - self.Mry_FLM()) * (self.esb_mesa - elp) / (elr - elp)
        elif self.esb_mesa > elr:
            return self.Mcry_FLM

    # Estado Limite FLA
    def Mny_FLA(self):
        return self.Mply
