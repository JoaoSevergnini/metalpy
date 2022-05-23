from math import sqrt, pi
from perfil_de_aco import PerfilAco, PerfilAcoNBR8800
from material import Material
import pandas as pd


class PerfilILaminadoNBR8800(PerfilAcoNBR8800):
    """
    Esta classe define um perfil de aço do tipo I laminado (W, H, HP).

    Parameter
    ----------
    nome: 'str'
        nome da seção transversal de acordo com a tabela da IASC

    material: Material, list, dict
        material ou propriedades do aço que compõe o perfil
    """

    def __init__(self, nome, material):

        perfil = self.perfis_laminados[self.perfis_laminados['EDI_Std_Nomenclature.1'] == nome]

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
        tipo = 'I Laminado'

        super().__init__(A, Ix, Iy, J, Wx, Wy, Zx, Zy, xo, yo, Cw, material, simetria, tipo)

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

    def par_esbeltez_limite_AL_Ncrd(self):
        return 0.56 * self.raiz_E_fy, 1.03 * self.raiz_E_fy

    def par_esbeltez_limite_AA_Ncrd(self):
        return 1.49 * self.raiz_E_fy

    def fator_Qs(self):

        # elp = esbeltez limite para plastificação
        # elr = esbeltez limite para início de escoamento
        elp, elr = self.par_esbeltez_limite_AL_Ncrd()

        if elp > self.esb_mesa:
            return 1
        if elp < self.esb_mesa <= elr:
            return 1.415 - 0.74 * self.esb_alma * sqrt(self.material.fy / self.material.E)
        elif self.esb_mesa > elr:
            return 0.69 * self.material.e / (self.material.fy * self.esb_mesa ** 2)

    def fator_Qa(self, frc):

        if self.esb_alma > self.par_esbeltez_limite_AA_Ncrd():
            return 1
        else:
            tensao = self.material.fy * frc
            ca = 0.34

            bef = 1.92 * self.tw * sqrt(self.material.E / tensao) * \
                (1 - ca / self.esb_alma * sqrt(self.material.E / tensao))

            bef = bef if bef < self.h else self.h

            Aef = self.A - (self.h - bef) * self.tw

            return Aef / self.A

    # CORTANTE
    # -----------

    def _Awx(self):
        return 2 * self.bf * self.tf

    def _Awy(self):
        return self.ht * self.tw

    # CORTANTE EM X
    def kv_Vrdx(self):
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

    # Estado Limite FLM
    def par_esbeltez_limite_Mrdx_FLM(self):
        return 0.38 * self.raiz_E_fy, 0.83 * sqrt(1 / 0.7) * self.raiz_E_fy

    def Mcrx_FLM(self):
        return 0.69 * self.material.E * self.Wx / self.esb_mesa ** 2

    # Estado Limite FLA
    def par_esbeltez_limite_Mrdx_FLA(self):
        return 3.76 * self.raiz_E_fy, 5.7 * self.raiz_E_fy

    # MOMENTO EM Y
    # ------------
    # Estado Limite FLT
    def Mny_FLT(self, Cb=None, Lb=None):
        return self.Mply

    # Estado Limite FLM
    def par_esbeltez_limite_Mrdy_FLM(self):
        return 0.38 * self.raiz_E_fy, 0.83 * sqrt(1 / 0.7) * self.raiz_E_fy

    def Mcry_FLM(self):
        return 0.69 * self.material.E / self.esb_mesa ** 2 * self.Wy

    # Estado Limite FLA
    def Mny_FLA(self):
        return self.Mply
