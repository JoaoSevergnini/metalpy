from math import sqrt, pi
from perfil_de_aco import PerfilDeAço
from material import Material
import pandas as pd

# ==================
# RENOMEAR A CLASSE
# ==================
class NomeDoTipoDePerfil(PerfilDeAço):
    """
    Esta classe define um perfil de aço laminado do tipo  XXXX.

    Parameter
    ----------
    nome: 'str'
        nome da seção transversal de acordo com a tabela da IASC

    material: 'Material', 'list', 'dict'
             material que compõe a seção.

    """

    # =================================================================================================
    # li = NÚMERO DA LINHA ONDE COMEÇA O TIPO DE SEÇÃO ESPECIFICA DE ACORDO COM A TABELA DA AISC -1
    # lf = NÚMERO DA LINHA ONDE TERMINA O TIPO DE SEÇÃO ESPECIFICA DE ACORDO COM A TABELA DA AISC + 1
    li = XX
    lf = XX
    # =================================================================================================

    perfis = pd.read_excel('aisc-shapes-database-v15.0.xlsx', 1).iloc[li:lf, 84:]

    def __init__(self, nome, material):

        perfil = self.perfis[self.perfis['EDI_Std_Nomenclature.1'] == nome]

        # ===============================================================
        # COLOCAR AQUI AS PROPRIEDADES ESPECIFICAS DO PERFIL A SER CRIADO
        # ===============================================================

        A = float(perfil['A.1'])
        Ix = float(perfil['Ix.1']) * 1E6
        Iy = float(perfil['Iy.1']) * 1E6
        J = float(perfil['J.1']) * 1E3
        Wx = float(perfil['Sx.1']) * 1E3
        Wy = float(perfil['Sy.1']) * 1E3
        Zx = float(perfil['Zx.1']) * 1E3
        Zy = float(perfil['Zy.1']) * 1E3
        Cw = float(perfil['Cw.1']) * 1E9

        # =========================================================
        # INDICAR A POSIÇÃO DO CENTRO DE CORTE (xo, yo) EM RELAÇÃO
        # AO CENTRO GEOMÉTRICO DA SEÇÃO
        # =========================================================
        xo = #FÓRMULA PARA O CALCULO DE XO
        yo = #FÓRMULA PARA O CALCULO DE YO

        # =======================================
        # INDICAR A SIMETRIA DO TIPO DE PERFIL
        #
        # True = simetrico
        # False = Não simetrico
        #
        # simetria[0] - indica simetria no eixo X
        # simetria[1] - indica simetria no eixo Y
        # =======================================
        simetria = [True, True]

        super().__init__(A, Ix, Iy, J, Wx, Wy, Zx, Zy, xo, yo, Cw, material, simetria)

        self.esb_alma = #FÓRMULA PARA O CALCULO DA ESBELTEZ DA ALMA DE ACORDO COM AS PROPRIEDADES DO PERFIL
        self.esb_mesa = #FÓRMULA PARA O CALCULO DA ESBELTEZ DA MESA DE ACORDO COM AS PROPRIEDADES DO PERFIL

    # ====================================
    # COLOCAR AQUI OS MÉTODOS PROPERTIES
    # ====================================


    # -------------------------------------------------------------------------------------
    # --------------------------Verificações de resistência--------------------------------
    # -------------------------------------------------------------------------------------

    # ----------------------------------NBR 8800-------------------------------------------
    # COMPRESSÃO
    # -----------

    def par_esbeltez_limites_AL_Ncrd(self):

        # =================================================================================================
        # DEVE SER IMPLEMENTADO APENAS NOS PERFIL QUE APRESENTEM ELEMENTOS
        # APOIADOS LIVRES ( seções do tipo I, C, L e T )

        # CASO NÃO EXISTA ELEMENTOS DESSE TIPO ESSE MÉTODO PODE SER EXCLUÍDO

        # coef1 = COEFICIENTE PARA O CÁLCULO DO LIMITE DE PLASTIFICAÇÃO DE ACORDO COM A SEÇÃO F.2 DA NBR8800
        # coef2 = COEFICIENTE PARA O CÁLCULO DO LIMITE DE ESCOAMENTO DE ACORDO COM A SEÇÃO F.2 DA NBR8800
        # =================================================================================================

        return coef1 * self.raiz_E_fy, coef2 * self.raiz_E_fy

    def fator_Qs(self):

        # =============================================================================
        # DEVE SER IMPLEMENTADO APENAS NOS PERFIL QUE APRESENTEM ELEMENTOS
        # APOIADOS LIVRES ( seções do tipo I, C, L e T )

        # CASO NÃO EXISTA ELEMENTOS DESSE TIPO ESSE MÉTODO DEVE RETORNAR O VALOR 1
        # =============================================================================

        # elp = esbeltez limite para plastificação
        # elr = esbeltez limite para início de escoamento
        elp, elr = self.par_esbeltez_limites_AL_Ncrd()

        # =============================================================================
        # coef3 = COEFICIENTE PARA O CÁLCULO DE QS DE ACORDO COM A SEÇÃO F.2 DA NBR8800
        # coef4 = COEFICIENTE PARA O CÁLCULO DE QS DE ACORDO COM A SEÇÃO F.2 DA NBR8800
        # coef5 = COEFICIENTE PARA O CÁLCULO DE QS DE ACORDO COM A SEÇÃO F.2 DA NBR8800
        # =============================================================================

        if elp > self.esb_mesa:
            return 1
        if elp < self.esb_mesa <= elr:
            return coef3 - coef4 * self.esb_alma * sqrt(self.material.fy / self.material.E)
        elif self.esb_mesa > elr:
            return coef5 * self.material.e / (self.material.fy * self.esb_mesa ** 2)

    def fator_Qa(self, frc):
        # ========================================================================
        # DEVE SER IMPLEMENTADO APENAS NOS PERFIL QUE APRESENTAM ELEMENTOS
        # APOIADOS APOIADOS ( seções do tipo I, C, e TUBULARES )

        # CASO NÃO EXISTA ELEMENTOS DESSE TIPO ESSE MÉTODO DEVE RETORNAR O VALOR 1
        # ========================================================================

        tensao = self.material.fy * frc

        # =================================
        # INDICAR O VALOR DO COEFICIENTE ca
        #
        # ca = 0.38 - TUBOS RETANGULARES
        # ca = 0.34 - OUTROS PERFIS
        # =================================
        # =====================================
        # ADAPTAR O VALOR DE ca ABAIXO CASO NECESSÁRIO
        #======================================
        ca = 0.34

        # ========================================================================
        # PARA CADA ELEMENTO APOIADO-APOIADO DA SEÇÃO TRANSVERSAL A LARGURA EFETIVA
        # DEVE SER CALCULADO COM A FÓRMULA DA LARGURA EFETIVA bef ABAIXO
        #
        # bef = 1.92 . t . sqrt(E / tensao) . [ 1 - ca / (b/t) . sqrt( E / tensao)]
        #
        # SENDO QUE O VALOR NÃO PODE ULTRAPASSAR A LARGURA
        # =========================================================================

        # =====================================================================================
        # ADAPTAR AS EQUAÇÕES ABAIXO PARA ATENDER TODOS OS ELEMENTOS APOIADOS APOIADOS DA SEÇÃO
        # E SUBSTITUIR  'b' E 't' PELA PROPRIEDADE DA CLASSE QUE REPRESENTAM A LARGURA E ESPESSURA
        # DO ELEMENTO AA CONSIDERADO
        # =======================================================================================

        bef = 1.92 * 't' * sqrt(self.material.E / tensao) * \
              (1 - ca / 'b/t' * sqrt(self.material.E / tensao))

        bef = bef if bef < 'b' else 'b'

        Aef = self.A - ('b' - bef) * 't'

        return Aef / self.A



    # CORTANTE
    # -----------

    @property
    def Awx(self):
        # ===================================================================================
        # RETORNAR O CÁLCULO Awx ESPECÍFICO DA SEÇÃO DE ACORDO COM A SEÇÃO 5.4.3.1 DA NBR8800
        # ====================================================================================
        return

    @property
    def Awy(self):
        # ===================================================================================
        # RETORNAR O CÁLCULO Awy ESPECÍFICO DA SEÇÃO DE ACORDO COM A SEÇÃO 5.4.3.1 DA NBR8800
        # ====================================================================================
        return #FORMULA DO Aw

    # CORTANTE EM X
    def kv_Vrdx(self, a=None):
        # =================================================================================================
        # RETORNAR O CÁLCULO kv ESPECÍFICO EM RELAÇÃO A X DA SEÇÃO DE ACORDO COM A SEÇÃO 5.4.3.1 DA NBR8800
        # =================================================================================================
        return # VALOR DO kv

    # CORTANTE EM Y
    def kv_Vrdy(self, a=None):
        # =================================================================================================
        # RETORNAR O CÁLCULO kv ESPECÍFICO EM RELAÇÃO A Y DA SEÇÃO DE ACORDO COM A SEÇÃO 5.4.3.1 DA NBR8800
        # =================================================================================================
            return # VALOR DO kv



    # MOMENTO EM X
    # ------------
    # Estado Limite FLT

    def par_esbeltez_limite_Mrdx_FLT(self):
        # =============================================================================
        # DEVE SER IMPLEMENTADO APENAS SE A SEÇÃO APRESENTAR O ESTADO LIMITE FLT QUANDO
        # FLETIDO NO EIXO X
        #
        # CASO NÃO EXISTA ELEMENTOS DESSE TIPO ESSE MÉTODO DEVE SER EXCLUÍDO
        # =============================================================================

        # ===========================================================================
        # OS PARAMETROS DE ESBELTEZ LIMITE DE PLASTIFICAÇÃO (elp) E O PARAMETRO DE
        # DE ESBELTEZ LIMITE DE INICIO DE ESCOAMENTO (elr) DEVEM SER IMPLEMENTADOS DE
        # ACORDO COM A TABELA G.1 DO ANEXO G DA NBR8800:2008
        #============================================================================

        # parâmetro de esbeltez limite de plastificação (elp)
        elp = #FÓRMULA PARA O CÁLCULO DO PARÂMETRO DE ESBELTEZ


        # parâmetro de esbeltez limite de início de escoamento (elr)
        elr = #FÓRMULA PARA O CÁLCULO DO PARÂMETRO DE ESBELTEZ

        return elp, elr

    def Mrx_FLT(self):
        # =============================================================================
        # DEVE SER IMPLEMENTADO APENAS SE A SEÇÃO APRESENTAR O ESTADO LIMITE FLT QUANDO
        # FLETIDO NO EIXO X
        #
        # CASO NÃO EXISTA ELEMENTOS DESSE TIPO ESSE MÉTODO DEVE SER EXCLUÍDO
        # =============================================================================

        # ================================================================================
        # O MOMENTO FLETOR DE INICIO DE ESCOAMENTO DE ACORDO COM A TABELA G.1 E SEÇÃO G.2
        # ANEXO G DA NBR8800:2008, É DADO PELO PRODUTO DA TENSÃO DE ESCOAMENTO (fy),
        # REDUZIDO PELA TENSÃO RESIDUAL, E O MÓDULO DE RESISTÊNCIA ELÁSTICO DA SEÇÃO
        # Mr = (fy - tensao residual) . W
        #
        # CASO EXISTA TENSÃO RESIDUAL, ELA É TOMADA COMO SENDO 30% DA TENSÃO DE ESCOAMENTO
        # SENDO ASSIM Mr = 0.7 . fy . W
        #
        # VER TABELA G1 PARA OBTER O CÁLCULO DA TENSÃO RESIDUAL ESPECIFICA DO PERFIL
        #=================================================================================

        return 0.7 * self.material.fy * self.Wx

    def Mcrx_FLT(self, Cb, Lb):

        # =============================================================================
        # DEVE SER IMPLEMENTADO APENAS SE A SEÇÃO APRESENTAR O ESTADO LIMITE FLT QUANDO
        # FLETIDO NO EIXO X
        #
        # CASO NÃO EXISTA ELEMENTOS DESSE TIPO ESSE MÉTODO DEVE SER EXCLUÍDO
        # =============================================================================

        # ======================================================================================
        # O MOMENTO FLETOR DE FLAMBAGEM ELÁSTICO DE ACORDO COM A TABELA G.1 DO
        # ANEXO G DA NBR8800:2008
        #
        # VER TABELA G1 E SEÇÃO G.2 PARA OBTER O CÁLCULO DA TENSÃO RESIDUAL ESPECÍFICA DO PERFIL
        #========================================================================================

        Mcr = # FÓRMULA PARA O CÁLCULO DO Mcr VER TABELA G.1 E SEÇÃO G.2 DO ANEXO G DA NBR8800:2008
        return Mcr

    def Mnx_FLT(self, Cb, Lb):

        # =============================================================================
        # DEVE SER IMPLEMENTADO APENAS SE A SEÇÃO APRESENTAR O ESTADO LIMITE FLT QUANDO
        # FLETIDO NO EIXO X
        #
        # CASO NÃO EXISTA ELEMENTOS DESSE TIPO ESSE MÉTODO DEVE RETORNAR O MOMENTO DE
        # PLASTIFICAÇÃO Mplx
        # ==============================================================================

        # ==========================================================================
        # PARA OS TIPOS DE PERFIS DA TABELA G.1 DO ANEXO G O PROCEDIMENTO DE CÁLCULO
        # É REALIZADO DE ACORDO COM O QUE ESTÁ IMPLEMENTADO ABAIXO, PARA OUTROS CASOS
        # VER ITENS DA SEÇÃO G.2 DO ANEXO G DA NBR8800:2008
        # =============================================================================

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
