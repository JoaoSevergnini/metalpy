from math import sqrt
from material import Material
from perfil_estrutural import PerfilEstrutural
import pandas as pd
import abc


class PerfilAco(PerfilEstrutural):
    """

    Esta classe define um seção tranversal de barra de formato genérico
    de acordo com suas propriedades geométricas e seu material.

    Parameter
    ----------

    A: 'float'
        área total da seção transversal

    Ix: 'float'
        momento de inércia a flexão do perfil em relação ao eixo X (horizontal)
        que passa pelo centroide da seção.

    Iy: 'float'
        momento de inércia a flexão da perfil em relação ao eixo Y (Vertical)
        que passa pelo centroide da seção.

    J: 'float'
        constante de torção da seção em relação ao centróide da seção

    aco: 'Material', 'list', 'dict'
        material que compõe a seção.

    Wx: 'float'
        módulo elástico do perfil em relação ao eixo X (horizontal)

    Wy: 'float'
        módulo elástico do perfil em relação ao eixo Y (Vertical)

    Zx: 'float'
        módulo plástico do perfil em relação ao eixo X (horizontal)

    xo: 'float'
        coordenada x do centro de corte da seção trasnversal em relação ao
        centróide da seção

    yo: 'float'
        coordenada y do centro de corte da seção trasnversal em relação ao
        centróide da seção

    Cw: 'float'
        constante de empenamento do perfil

    simetria: 'list'
        indica se o perfil apresenta eixos de simetria

    tipo: 'string'
        tipo do perfil
    """

    # Data frame com todos os perfis laminados e suas propriedades, catalogados pela AISC
    perfis_laminados = pd.read_excel('aisc-shapes-database-v15.0.xlsx', 1).iloc[:, 84:]

    def __init__(self, A, Ix, Iy, J, Wx, Wy, Zx, Zy, xo, yo, Cw, aco, simetria, tipo):
        self.tipo = tipo

        super().__init__(A, Ix, Iy, J, Wx, Wy, Zx, Zy, self._Awx, self._Awy, xo, yo, Cw, aco, simetria)

        self.esb_alma = None
        self.esb_mesa = None

        self.raiz_E_fy = sqrt(self.material.E / self.material.fy)
        self.raiz_fy_E = sqrt(self.material.fy / self.material.E)

    # -------------------------------------------------------------------------------------
    # --------------------------Verificações de resistência--------------------------------
    # -------------------------------------------------------------------------------------

    # TRAÇÂO
    # --------
    def Ntrd_ESB(self, coef_seg):
        """
        Método que determina a resistência ao escoamento da seção bruta de perfil métálico
        """
        return NotImplementedError('Método não implementado')

    # COMPRESSÃO
    # -----------
    @property
    def par_esbeltez_limite_AL_Ncrd(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e
        de início de escoamento de elementos apoiados livres em compressão
        """

        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil, caso o perfil
        # a ser criado não apresente elementos apoiados livres este
        # método não precisa ser implementado.

        raise NotImplementedError('Método não implementado')

    @property
    def par_esbeltez_limite_AA_Ncrd(self):
        """
        Retorna o parâmetro de esbeltez de início de escoamento de
        elementos apoiados apoiados em compressão
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil, caso o perfil
        # a ser criado não apresente elementos apoiados apoiados este método
        # não precisa ser implementado.

        raise NotImplementedError('Método não implementado')

    def Ncrd(self, klx, kly, klz, coef_seg):
        """
        Método que determina a resistência a compressão de cálculo de uma barra de aço
        """
        raise NotImplementedError('Método não implementado')

    # CORTANTE
    # -----------
    def par_esbeltez_limite_Vrd(self, kv):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de início de escoamento
        dos elementos solicitados por esforço cortante
        """
        raise NotImplementedError('Método não implementado')

    # cortante em X
    # -------------
    def _Awx(self):
        """ Área efetiva de cisalhamento na direção X"""

        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def kv_Vrdx(self):
        """
        Retorna o valor do coeficiente kv do perfil para determinação da força
        resistênte ao corte na direção x.

        Coeficiente que leva em consideração a existência de enrijecedores
        e a forma do perfil na resistência ao cisalhamento
        """

        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Vrdx(self, coef_seg):
        """
        Deve retornar a força cortante resistente de cálculo do perfil de aço,
        com carga aplicada paralela ao eixo y do perfil
        """
        raise NotImplementedError('Método não implementado')

    # cortante em Y
    # -------------
    def _Awy(self):
        """ Área efetiva de cisalhamento na direção Y"""

        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def kv_Vrdy(self, a=None):
        """
        Retorna o valor do coeficiente kv do perfil para determinação da
        força resistênte ao corte na direção y.

        Coeficiente que leva em consideração a existência de enrijecedores
        e a forma do perfil na resistência ao cisalhamento

        Patameter
        ---------
        a: 'float'
            distância entre os centros dos enrijecedores.

        Return
        ------
        kv: 'float'
            coeficiente kv
        """

        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Vrdy(self, coef_seg, a=None):
        """
        Deve retornar a força cortante resistente de cálculo do perfil de aço,
        com carga aplicada paralela ao eixo y do perfil
        """

        raise NotImplementedError('Método não implementado')

    # MOMENTO FLETOR EM X
    # ------------
    # Estado Limite FLT

    def par_esbeltez_limite_Mrdx_FLT(self):
        """
        Deve retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem lateral com
        torção para barras fletidas em relação ao eixo x.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mrx_FLT(self):
        """
        Retorna o momento fletor em X correspondente ao início de
        escoamento da seção, para o estado limite de flambagem lateral
        com torção.
        """
        raise NotImplementedError('Método não implementado')

    def Mcrx_FLT(self, cb, lb):
        """
        Retorna o momento fletor em X correspondente a flambagem elástica,
        para o estado limite de flambagem lateral com torção.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mnx_FLT(self, Cb, Lb):
        """
        Determina o momento fletor resistente nominal de uma barra para
        o estado limite último de flambagem lateral com torção.

        Parameter
        --------
        Cb: 'float'
            coeficiente Cb

        Lb: 'float'
            comprimento destravado da barra

        Return
        ------
        """

        # Para perfis que não apresentam a flambagem lateral com torção
        # com um estado limite este método deve ser implementado no perfil
        # especifico retornando o momento de plastificação do perfil (Mpl)

        raise NotImplementedError('Método não implementado')

    # Estado Limite FLM
    def par_esbeltez_limite_Mrdx_FLM(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem local da
        mesa para barras fletidas em relação ao eixo x.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mrx_FLM(self):
        """
        Retorna o momento fletor em X correspondente ao inicio de escoamento da seção,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mcrx_FLM(self):
        """
        Retorna o momento fletor em X correspondente a flambagem elástica,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mnx_FLM(self):
        """
        Determina o momento fletor resistente nominal em X de uma barra para
        o estado limite último de flambagem local da mesa.

        Return
        ------
        """
        # Para perfis que não apresentam a flambagem local da mesa
        # como um estado limite este método deve ser implementado no perfil
        # especifico retornando o momento de plastificação do perfil (Mpl)

        raise NotImplementedError('Método não implementado')

    # Estado Limite FLA
    def par_esbeltez_limite_Mrdx_FLA(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem local da
        alma para barras fletidas em relação ao eixo x.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mrx_FLA(self):
        """
        Retorna o momento fletor em X correspondente ao inicio de escoamento da seção,
        para o estado limite de flambagem local da alma.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mcrx_FLA(self):
        """
        Retorna o momento fletor em X correspondente a flambagem elástica,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mnx_FLA(self):
        """
        Determina o momento fletor resistente nominal em X de uma barra para
        o estado limite último de flambagem local da alma.

        Return
        ------
        """
        # Para perfis que não apresentam a flambagem local da alma
        # com um estado limite este método deve ser implementado no perfil
        # especifico retornando o momento de plastificação do perfil (Mpl)

        raise NotImplementedError('Método não implementado')

    def Mrdx(self, Lb, coef_seg, Cb=1):
        """
        Método responsável por calcular o momento fletor resitente de cálculo para uma
        barra de comprimento Lb em relação ao eixo X do perfil.

        Parameter
        --------
        Cb: 'float'
            coeficiente Cb

        Lb: 'float'
            comprimento destravado da barra

        coef_seg: 'float'
            coeficiente de minoração da resistência

        Return
        ------

        """
        raise NotImplementedError('Método não implementado')

    # MOMENTO EM Y
    # ------------
    # Estado Limite FLT

    def par_esbeltez_limite_Mrdy_FLT(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem lateral com
        torção para barras fletidas em relação ao eixo y.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.
        raise NotImplementedError('Método não implementado')

    def Mry_FLT(self):
        """
        Retorna o momento fletor em Y correspondente ao início de escoamento da seção,
        para o estado limite de flambagem lateral com torção.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mcry_FLT(self, Cb, Lb):
        """
        Retorna o momento fletor em Y correspondente a flambagem elástica,
        para o estado limite de flambagem lateral com torção.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mny_FLT(self, Cb, Lb):
        """
        Determina o momento fletor resistente nominal em Y de uma barra para
        o estado limite último de flambagem lateral com torção.

        Return
        ------
        """
        # Para perfis que não apresentam a flambagem lateral com torção
        # com um estado limite este método deve ser implementado no perfil
        # especifico retornando o momento de plastificação do perfil (Mpl)

        raise NotImplementedError

    # Estado Limite FLM
    def par_esbeltez_limite_Mrdy_FLM(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem local
        da mesa para barras fletidas em relação ao eixo y conforme o que
        consta na seção G2 do anexo G da NBR8800:2008.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mry_FLM(self):
        """
        Retorna o momento fletor em Y correspondente ao início de escoamento
        da seção, para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.
        raise NotImplementedError

    def Mcry_FLM(self):
        """
        Retorna o momento fletor em Y correspondente a flambagem elástica,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mny_FLM(self):
        """
        Determina o momento fletor resistente nominal em Y de uma barra para
        o estado limite último de flambagem local da mesa.

        Return
        ------
        """
        # Para perfis que não apresentam a flambagem local da mesa
        # com um estado limite este método deve ser implementado no perfil
        # especifico retornando o momento de plastificação do perfil (Mpl)

        raise NotImplementedError

    # Estado Limite FLA

    def par_esbeltez_limite_Mrdy_FLA(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem local
        da alma para barras fletidas em relação ao eixo y.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mry_FLA(self):
        """
        Retorna o momento fletor em Y correspondente ao início de escoamento
        da seção, para o estado limite de flambagem local da alma.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mcry_FLA(self):
        """
        Retorna o momento fletor em Y correspondente a flambagem elástica,
        para o estado limite de flambagem local da alma.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mny_FLA(self):
        """
        Determina o momento fletor resistente nominal em Y de uma barra para
        o estado limite último de flambagem local da alma.

        Return
        ------
        """
        # Para perfis que não apresentam a flambagem local da alma
        # com um estado limite este método deve ser implementado no perfil
        # especifico retornando o momento de plastificação do perfil (Mpl)
        raise NotImplementedError

    def Mrdy(self, Lb, coef_seg, Cb=1):
        """
        Método responsável por calcular o momento fletor resitente de cálculo para uma
        barra de comprimento Lb em relação ao eixo Y do perfil.

        Parameter
        --------
        Cb: 'float'
            coeficiente Cb determinado.

        Lb: 'float'
            comprimento destravado da barra

        coef_seg: 'float'
            coeficiente de minoração da resistência

        Return
        ------

        """

        raise NotImplementedError


class PerfilAcoNBR8800(PerfilAco):
    """

    Esta classe define um seção tranversal de barra de formato genérico
    de acordo com suas propriedades geométricas e seu material.

    Parameter
    ----------

    A: 'float'
        área total da seção transversal

    Ix: 'float'
        momento de inércia a flexão do perfil em relação ao eixo X (horizontal)
        que passa pelo centroide da seção.

    Iy: 'float'
        momento de inércia a flexão da perfil em relação ao eixo Y (Vertical)
        que passa pelo centroide da seção.

    J: 'float'
        constante de torção da seção em relação ao centróide da seção

    material: 'Material', 'list', 'dict'
        material que compõe a seção.

    Wx: 'float'
        módulo elástico do perfil em relação ao eixo X (horizontal)

    Wy: 'float'
        módulo elástico do perfil em relação ao eixo Y (Vertical)

    Zx: 'float'
        módulo plástico do perfil em relação ao eixo X (horizontal)

    xo: 'float'
        coordenada x do centro de corte da seção trasnversal em relação ao
        centróide da seção

    yo: 'float'
        coordenada y do centro de corte da seção trasnversal em relação ao
        centróide da seção

    Cw: 'float'
        constante de empenamento do perfil

    simetria: 'list'
        indica se o perfil apresenta eixos de simetria
    """

    def __init__(self, A, Ix, Iy, J, Wx, Wy, Zx, Zy, xo, yo, Cw, aco, simetria, tipo):

        super().__init__(A, Ix, Iy, J, Wx, Wy, Zx, Zy, xo, yo, Cw, aco, simetria, tipo)

        self.c_tensao_res = 0.7

    # -------------------------------------------------------------------------------------
    # --------------------------Verificações de resistência--------------------------------
    # -------------------------------------------------------------------------------------

    # --------------------------------NBR8800/2008-----------------------------------------

    # TRAÇÂO
    # --------

    def Ntrd_ESB(self, gama_a1=1.1):
        """
        Método que determina a resistência ao escoamento da seção bruta de perfil métálico

        A resitência ao escoamento da seção bruta é determinada de acordo com o item a) da
        seção 5.2.2 da NBR8800:2008, que é obtida pelo produto da área total do perfil pelo
        tensão de escoamento, sendo o valor obtido deste produto minorado pelo coeficiente de
        segurança gama a1.

        Parameter
        ---------
        gama_a1: 'float' (default=1,1)
                coeficiente de segurança gama_a1

        Return
        ------

        """
        return self.Afy / gama_a1

    # COMPRESSÃO
    # -----------
    def ind_esbeltez_reduzido(self, klx, kly, klz, Q=1):
        """
        Método que determina o indice de elbeltez reduzido de uma barra de aço de determinado perfil.

        O indice de esbeltez reduzido é determinado de acordo com o item 5.3.3.2 da NBR8800:2008,
        que apresenta um limite superior para o seu valor de 200.


        Parameter
        ---------
        klx:'float'
            comprimento de flambagem por flexão em relação ao eixo x

        kly:'float'
            comprimento de flambagem por flexão em relação ao eixo Y

        klz:'float'
            comprimento de flambagem por torção em relação ao eixo longitudinal Z

        Q:'float' (default = 1)
          fator de redução total associado a flambagem local

        Return
        ------
        ier: 'float'
            indice de esbeltez reduzido
        """
        Ne = self.Ne(klx, kly, klz)
        ier = sqrt(Q * self.Afy / Ne)
        return ier

    def fator_reducao_compressao(self, klx, kly, klz, Q=1):
        """
        Método que determina o fator de redução da resistência a compressão do perfil (Fator Chi).

        O fator de redução da resistência a compressão (Fator Chi) é determinado de acordo
        com o item 5.3.3 da NBR8800:2008 em função do indíce de esbeltez reduzido da barra.

        Parameter
        --------

        ier: 'float'
            indice de esbeltez reduzido

        Return:
        -------
        frc: 'float'
              fator de redução de compressão

        """
        ier = self.ind_esbeltez_reduzido(klx, kly, klz, Q)
        if ier <= 1.5:
            frc = 0.658 ** (ier ** 2)
        else:
            frc = 0.877 / ier ** 2
        return frc

    def fator_Qs(self):
        """
        Método que determina o fator de redução de resistência a compressão,
        associado a flambagem local de elementos apoiados livres.

        O fator Qs é determinado de acordo com o item F.2  do anexo F da NBR8800:2008.
        """

        # Como o cálculo deste fator é em função do tipo de seção, este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil, caso o perfil o perfil a ser criado não apresente elementos
        # apoiados livres este método deve ser implementado retornando o valor 1.

        raise NotImplementedError

    def fator_Qa(self, frc):
        """
        Método que determina o fator de redução de resistência a compressão,
        associado a flambagem local de elementos apoiados apoiados.

        O fator Qa é determinado de acordo com o item F.3 do anexo F da NBR8800:2008.

        Como o cálculo deste fator é em função do tipo de seção, este método deve
        ser implementado em cada uma das classes especificas de cada tipo de perfil.
        """

        # Como o cálculo deste fator é em função do tipo de seção, este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil, caso o perfil o perfil a ser criado não apresente elementos
        # apoiados apoiados este método deve ser implementado retornando o valor 1.

        raise NotImplementedError

    def fator_Q(self, frc):
        """
        Método que determina o fator de redução total associado a flambagem local.

        O fator de redução Q associado a flambagem local é determinado de acordo com o anexo F da NBR8800:2008, sendo
        determinada pela multiplicação dos fatores Qa e Qs que estão associados, respectivamente, a flambagem local
        dos elementos apoiados-apoiados(AA), como almas de perfis I por exemplo, e elementos apoiados-livres(AL), como
        mesas de perfis I.

        Parameter
        --------

        frc: 'float'
            fator de redução da compressão

        Return:
            'float'
        """
        return self.fator_Qa(frc) * self.fator_Qs()

    def Ncrd(self, klx, kly, klz, gama_a1=1.1):
        """
        Método que determina a resistência a compressão de cálculo
        de uma barra de aço de acordo com a NBR8800

        Parameter
        ---------

        klx:'float'
            comprimento de flambagem por flexão em relação ao eixo x

        kly:'float'
            comprimento de flambagem por flexão em relação ao eixo Y

        klz:'float'
            comprimento de flambagem por torção em relação ao eixo
            longitudinal Z

        Return
        ------
        Ncrd = 'float'
            resistência a compressão do perfil
        """

        frc = self.fator_reducao_compressao(klx, kly, klz)

        Q = self.fator_Q(frc)
        frc = self.fator_reducao_compressao(klx, kly, klz, Q)

        Ncrd = frc * Q * self.Afy / gama_a1
        return Ncrd

    # CORTANTE
    # -----------
    def par_esbeltez_limite_Vrd(self, kv):
        return 1.1 * self.raiz_E_fy * sqrt(kv)

    # EM X
    # -------------
    def Vrdx(self, gama_a1=1.1, a=None):
        """
        Método que determina a força cortante resistente de cálculo na
        direção X do perfil de acordo com a NBR8800:2008.

        O procedimento para a determinação da capacidade resistênte ao corte da
        seção transversal é realizado conforme a seção 5.4.3 da NBR8800:2008.

        Parameter
        --------
        a: 'float'
            distância entre eixos de enrijecedores

        gama_a1: 'float' (default = 1.1)
            coeficiente de minoração da resistência.

        Return
        ------
        Vrdx: 'float'
            Força cortante resistênte de cálculo na direção x.
        """

        kv = self.kv_Vrdx()

        # elp = esbeltez limite para plastificação
        # elr = esbeltez limite para início de escoamento
        elp, elr = self.par_esbeltez_limite_Vrd(kv)

        if self.esb_mesa <= elp:
            return self.Vplx / gama_a1

        elif elp < self.esb_mesa <= elr:
            return (elp / self.esb_mesa) * (self.Vplx / gama_a1)

        else:
            return 1.24 * (elp / self.esb_mesa) ** 2 * (self.Vplx / gama_a1)

    # CORTANTE EM Y
    # -------------
    def Vrdy(self, a=None, gama_a1=1.1):

        """
        Método que determina a força cortante resistente de cálculo na
        direção Y do perfil de acordo com a NBR8800:2008.

        O procedimento para a determinação da capacidade resistente ao corte
        da seção transversal é realizado conforme a seção 5.4.3 da NBR8800:2008.

        Parameter
        --------
        a: 'float'
            distância entre eixos de enrijecedores

        gama_a1: 'float' (default = 1.1)
            coeficiente de minoração da resistência.

        Return
        ------
        Vrdy: 'float'
            Força cortante resistênte de cálculo na direção y
        """

        kv = self.kv_Vrdy(a)

        # elp = esbeltez limite para plastificação
        # elr = esbeltez limite para início de escoamento
        elp, elr = self.par_esbeltez_limite_Vrd(kv)

        if self.esb_alma <= elp:
            return self.Vply / gama_a1

        elif elp < self.esb_alma <= elr:
            return (elp / self.esb_alma) * (self.Vply / gama_a1)

        else:
            return 1.24 * (elp / self.esb_alma) ** 2 * (self.Vply / gama_a1)

    # MOMENTO FLETOR EM X
    # ------------
    # Estado Limite FLT
    def Mrx_FLT(self):
        """
        Retorna o momento fletor em X correspondente ao início de escoamento da seção,
        para o estado limite de flambagem lateral com torção.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        return self.Mrx * self.c_tensao_res

    def Mcrx_FLT(self, Cb, Lb):
        """
        Retorna o momento fletor em X correspondente a flambagem elástica,
        para o estado limite de flambagem lateral com torção.

        Cb: 'float'
            coeficiente Cb determinado conforme item 5.4.2.3 da NBR8800:2008

        Lb: 'float'
            comprimento destravado da barra

        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        return Cb * self.Mex(Lb, Lb)

    def Mnx_FLT(self, Cb, Lb):
        """
        Determina o momento fletor resistente nominal de uma barra para
        o estado limite último de flambagem lateral com torção.

        Parameter
        --------
        Cb: 'float'
            coeficiente Cb determinado conforme item 5.4.2.3 da NBR8800:2008

        Lb: 'float'
            comprimento destravado da barra

        Return
        ------

        """

        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        esb = self.indice_esbeltez_X(Lb)
        elp, elr = self.par_esbeltez_limite_Mrdx_FLT()

        if esb < elp:
            return self.Mplx
        elif elp < esb < elr:
            return self.Mplx - (self.Mplx - self.Mrx_FLT) * (esb - elp) / (elr - elp)
        elif elr < esb:
            Mcrx = self.Mcrx_FLT(Cb, Lb)
            return Mcrx if Mcrx < self.Mplx else self.Mplx

    # Estado Limite FLM
    def Mrx_FLM(self):
        """
        Retorna o momento fletor em X correspondente ao inicio de escoamento da seção,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        return self.Mrx * self.c_tensao_res

    def Mcrx_FLM(self):
        """
        Retorna o momento fletor em X correspondente a flambagem elástica,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mnx_FLM(self):
        """
        Determina o momento fletor resistente nominal em X de uma barra para
        o estado limite último de flambagem local da mesa.

        Return
        ------
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        elp, elr = self.par_esbeltez_limite_Mrdx_FLM()

        if self.esb_alma < elp:
            return self.Mplx
        elif elp < self.esb_alma < elr:
            return self.Mplx - self.Mplx - (self.Mplx - self.Mrx_FLM) * (self.esb_alma - elp) / (elr - elp)
        elif self.esb_alma > elr:
            return self.Mcrx_FLM()

    # Estado Limite FLA
    def Mrx_FLA(self):
        """
        Retorna o momento fletor em X correspondente ao inicio de escoamento da seção,
        para o estado limite de flambagem local da alma.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.
        return self.Mrx

    def Mcrx_FLA(self):
        """
        Retorna o momento fletor em X correspondente a flambagem elástica,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError('Método não implementado')

    def Mnx_FLA(self):
        """
        Determina o momento fletor resistente nominal em X de uma barra para
        o estado limite último de flambagem local da alma.

        Return
        ------
        """
        elp, elr = self.par_esbeltez_limite_Mrdx_FLA()

        if self.esb_alma < elp:
            return self.Mplx
        elif elp < self.esb_alma < elr:
            return self.Mplx - self.Mplx - (self.Mplx - self.Mrx_FLM) * (self.esb_alma - elp) / (elr - elp)

    def Mrdx(self, Lb, gama_a1=1.1, Cb=1):
        """
        Método responsável por calcular o momento fletor resitente de cálculo para uma
        barra de comprimento Lb em relação ao eixo X do perfil, de acordo com a NBR8800.

        Parameter
        --------
        Cb: 'float'
            coeficiente Cb determinado conforme item 5.4.2.3 da NBR8800:2008

        Lb: 'float'
            comprimento destravado da barra

        gama_a1: 'float'
            coeficiente de minoração da resistência

        Return
        ------

        """
        return min(self.Mnx_FLA(), self.Mnx_FLM(), self.Mnx_FLT(Cb, Lb)) / gama_a1

    # MOMENTO EM Y
    # ------------

    # Estado Limite FLT
    def Mry_FLT(self):
        """
        Retorna o momento fletor em Y correspondente ao início de escoamento da seção,
        para o estado limite de flambagem lateral com torção.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        return self.Mry * self.c_tensao_res

    def Mcry_FLT(self, Cb, Lb):
        """
        Retorna o momento fletor em Y correspondente a flambagem elástica,
        para o estado limite de flambagem lateral com torção.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        return Cb * self.Mey(Lb, Lb)

    def Mny_FLT(self, Cb, Lb):
        """
        Determina o momento fletor resistente nominal em Y de uma barra para
        o estado limite último de flambagem lateral com torção.

        Return
        ------
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.
        esb = self.indice_esbeltez_Y(Lb)
        elp, elr = self.par_esbeltez_limite_Mrdy_FLT()

        if esb > elp:
            return self.Mply
        elif elp < esb < elr:
            return Cb * (self.Mply - (self.Mply - self.Mry_FLT()) * (esb - elp) / (elr - elp))
        else:
            return self.Mcry_FLT(Cb, Lb)

    # Estado Limite FLM
    def Mry_FLM(self):
        """
        Retorna o momento fletor em Y correspondente ao início de escoamento da seção,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.
        return self.c_tensao_res * self.Mry

    def Mcry_FLM(self):
        """
        Retorna o momento fletor em Y correspondente a flambagem elástica,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mny_FLM(self):
        """
        Determina o momento fletor resistente nominal em Y de uma barra para
        o estado limite último de flambagem local da mesa.

        Return
        ------
        """
        elp, elr = self.par_esbeltez_limite_Mrdy_FLM()

        if self.esb_mesa < elp:
            return self.Mply
        elif elp < self.esb_mesa < elr:
            return self.Mply - (self.Mply - self.Mry_FLM()) * (self.esb_mesa - elp) / (elr - elp)
        elif elr < self.esb_mesa:
            return self.Mcry_FLM()

    # Estado Limite FLA
    def Mry_FLA(self):
        """
        Retorna o momento fletor em Y correspondente ao início de escoamento da seção,
        para o estado limite de flambagem local da alma.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mcry_FLA(self):
        """
        Retorna o momento fletor em Y correspondente a flambagem elástica,
        para o estado limite de flambagem local da alma.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mny_FLA(self):
        """
        Determina o momento fletor resistente nominal em Y de uma barra para
        o estado limite último de flambagem local da alma.

        Return
        ------
        """
        elp, elr = self.par_esbeltez_limite_Mrdy_FLA()

        if self.esb_alma > elp:
            return self.Mply
        elif elp > self.esb_alma > elr:
            return self.Mply - (self.Mply - self.Mry_FLA()) * (self.esb_mesa - elp) / (elr - elp)

    def Mrdy(self, Lb, Cb=1, gama_a1=1.1):
        """
        Método responsável por calcular o momento fletor resitente de cálculo para uma
        barra de comprimento Lb em relação ao eixo Y do perfil, de acordo com a NBR8800.

        Parameter
        --------
        Cb: 'float'
            coeficiente Cb determinado conforme item 5.4.2.3 da NBR8800:2008

        Lb: 'float'
            comprimento destravado da barra

        gama_a1: 'float'
            coeficiente de minoração da resistência

        Return
        ------

        """

        return min(self.Mny_FLA(), self.Mny_FLM(), self.Mny_FLT(Cb, Lb)) / gama_a1


class PerfilAcoAISC360(PerfilAco):

    def __init__(self, A, Ix, Iy, J, Wx, Wy, Zx, Zy, xo, yo, Cw, aco, simetria, tipo):
        super().__init__(A, Ix, Iy, J, Wx, Wy, Zx, Zy, xo, yo, Cw, aco, simetria, tipo)

    # TRAÇÃO
    def Ntrd_ESB(self, phi_s=0.9):
        return self.Afy * phi_s

    # COMPRESSÃO
    def Aef(self):
        raise NotImplementedError('Método não implementado')

    def Fcr(self, klx, kly, klz):

        Fe = self.fe(klx, kly, klz)

        if self.material.fy / Fe <= 2.25:
            return 0.658 ** (self.material.fy / Fe) * self.material.fy
        else:
            return 0.877 * Fe

    def Ncrd(self, klx, kly, klz, phi_c=0.9):
        return self.Fcr(klx, kly, klz) * self.Aef() * phi_c

    # CORTANTE
    def par_esbeltez_limite_Vrd(self, kv):
        return 1.1 * sqrt(kv) * self.raiz_E_fy, 1.37 * sqrt(kv) * self.raiz_E_fy

    def Cv2_Vrdx(self):

        kv = self.kv_Vrdx()
        elp, elr = self.par_esbeltez_limite_Vrd(kv)

        if self.esb_mesa <= elp:
            return 1
        elif elp < self.esb_mesa <= elr:
            return elp/self.esb_mesa
        else:
            return 1.51 * kv * (self.raiz_E_fy / self.esb_mesa) ** 2

    def Vrdx(self, phi_v=0.90):
        return self.Vplx * self.Cv2_Vrdx() * phi_v

    def Cv2_Vrdy(self, a):

        kv = self.kv_Vrdy(a)
        elp, elr = self.par_esbeltez_limite_Vrd(kv)

        if self.esb_alma <= elp:
            return 1
        elif elp < self.esb_alma <= elr:
            return elp/self.esb_alma
        else:
            return 1.51 * kv * (self.raiz_E_fy / self.esb_alma) ** 2

    def Vrdy(self, phi_v=0.90, a=None):
        return self.Vply * self.Cv2_Vrdy(a) * phi_v

    # MOMENTO FLETOR
