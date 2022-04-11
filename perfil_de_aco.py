from math import sqrt
from secao import SecaoGenerica


class PerfilDeAço(SecaoGenerica):
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

    material: 'Material'
        material que compõe a seção em relação ao eixo X (horizontal)
        que passa pelo centroide da seção.

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
        constante de empenamento do pefil

    simetria:
        indica se o perfil apresenta eixos de simetria
    """

    def __init__(self, A, Ix, Iy, J, Wx, Wy, Zx, Zy, xo, yo, Cw, material, simetria):

        super().__init__(A, Ix, Iy, J, material, Wx, Wy, xo, yo, Cw, simetria)

        self.Zx = Zx
        self.Zy = Zy

        self.esb_alma = None
        self.esb_mesa = None

        self.raiz_E_fy = sqrt(self.material.E / self.material.fy)
        self.raiz_fy_E = sqrt(self.material.fy / self.material.E)

    # -------------------------------------------------------------------------------------
    # --------------------------Verificações de resistência--------------------------------
    # -------------------------------------------------------------------------------------

    # --------------------------------NBR8800/2008-----------------------------------------

    # TRAÇÂO
    # --------

    def resist_esc_secao_bruta_NBR8800(self, gama_a1=1.1):
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
        return self.A * self.material.fy / gama_a1

    # COMPRESSÃO
    # -----------
    def par_esbeltez_limites_AL_Ncrd(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de início de escoamento,
        respectivamente dos elementos apoiado-livre que compõe o perfil, de acordo com os itens
        F.2.a) à F.2.d) do anexo F da NBR8800:2008.
        """
        # Como o cálculo deste fator é em função do tipo de seção, este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil, caso o perfil o perfil a ser criado não apresente elementos
        # apoiados livres este método não precisa ser implementado.

        raise NotImplementedError

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
        ier = sqrt(Q * self.A * self.material.fy) / Ne
        return ier

    def fator_reducao_compressao(self, ier):
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
        if ier <= 1.5:
            frc = 0.658 ** (ier ** 2)
        else:
            frc = 0.877 / ier ** 2
        return frc

    def fator_Qs(self):
        """ Método que determina o fator de redução de resistência a compressão, associado a flambagem local de
        de elementos apoiados livres.

        O fator Qs é determinado de acordo com o item F.2  do anexo F da NBR8800:2008.
        """

        # Como o cálculo deste fator é em função do tipo de seção, este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil, caso o perfil o perfil a ser criado não apresente elementos
        # apoiados livres este método deve ser implementado retornando o valor 1.

        raise NotImplementedError

    def fator_Qa(self, frc):
        """
        Método que determina o fator de redução de resistência a compressão, associado a flambagem local de
        de elementos apoiados apoiados.

        O fator Qa é determinado de acordo com o item F.3 do anexo F da NBR8800:2008.

        Como o cálculo deste fator é em função do tipo de seção, este método deve ser implementado em cada uma das
        classes especificas de cada tipo de perfil.
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

    def Ncrd_NBR8800(self, klx, kly, klz, gama_a1=1.1):
        """
        Método que determina a resistência a compressão de cálculo de uma barra de aço de acordo com a NBR8800

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
        Ncrd = 'float'
            resistência a compressão do perfil
        """

        ier = self.ind_esbeltez_reduzido(klx, kly, klz)
        frc = self.fator_reducao_compressao(ier)
        Q = self.fator_Q(frc)

        ier = self.ind_esbeltez_reduzido(klx, kly, klz, Q)
        frc = self.fator_reducao_compressao(ier)

        Ncrd = frc * Q * self.A * self.material.fy / gama_a1
        return Ncrd

    # CORTANTE
    # -----------
    @property
    def Awx(self):
        """ Área efetiva de cisalhamento na direção X"""

        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    @property
    def Awy(self):
        """ Área efetiva de cisalhamento na direção Y"""

        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def par_esbeltez_limites_Vrd(self, kv):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de início de escoamento,
        respectivamente, conforme o item 5.4.3.1.1 da NBR8800:2008.

        Parameter
        --------
        kv: 'float'
            coeficiente kv do perfil

        Returns: 'float'
        """
        return 1.1 * sqrt(kv) * self.raiz_E_fy, 1.37 * sqrt(kv) * self.raiz_E_fy

    def Vpl(self, Aw):
        """
        Método que determina a força cortante de plastificação da alma.

        O valor de Vpl é determinado conforme o item 5.4.3.1.2 na NBR8800:2008.

        Parameter
        --------
        Aw: 'float'
            área efetiva de cisalhamento

        Return
        ------
        vpl: 'float'
            força cortante de plastificação
        """
        vpl = 0.60 * Aw * self.material.fy
        return vpl

    # EM X
    # -------------
    def kv_Vrdx(self, a):
        """
        Retorna o valor do coeficiente kv do perfil determinação da força
        resistênte ao corte na direção x.

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
        raise NotImplementedError

    def Vrdx_NBR8800(self, a, gama_a1=1.1):
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

        kv = self.kv_Vrdx(a)

        # elp = esbeltez limite para plastificação
        # elr = esbeltez limite para início de escoamento
        elp, elr = self.par_esbeltez_limites_Vrd(kv)

        if self.esb_mesa <= elp:
            return self.Vpl(self.Awx) / gama_a1

        elif elp < self.esb_mesa <= elr:
            return (elp / self.esb_mesa) * (self.Vpl(self.Awx) / gama_a1)

        else:
            return 1.24 * (elp / self.esb_mesa) ** 2 * (self.Vpl(self.Awx) / gama_a1)

    # CORTANTE EM Y
    # -------------
    def kv_Vrdy(self, a):
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
        raise NotImplementedError

    def Vrdy_NBR8800(self, a, gama_a1=1.1):

        """
        Método que determina a força cortante resistente de cálculo na
        direção X do perfil de acordo com a NBR8800:2008.

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
        elp, elr = self.par_esbeltez_limites_Vrd(kv)
        self.Vpl(self.Awy)

        if self.esb_alma <= elp:
            return self.Vpl(self.Awy) / gama_a1

        elif elp < self.esb_alma <= elr:
            return (elp / self.esb_alma) * (self.Vpl(self.Awy) / gama_a1)

        else:
            return 1.24 * (elp / self.esb_alma) ** 2 * (self.Vpl(self.Awy) / gama_a1)

    # MOMENTO FLETOR EM X
    # ------------
    @property
    def Mplx(self):
        """ Momento de plastificação da seção em relação ao eixo X"""
        return self.Zx * self.material.fy

    # Estado Limite FLT
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

    def par_esbeltez_limite_Mrdx_FLT(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem lateral com
        torção para barras fletidas em relação ao eixo x conforme o que consta na
        seção G2 do anexo G da NBR8800:2008.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mrx_FLT(self):
        """
        Retorna o momento fletor em X correspondente ao início de escoamento da seção,
        para o estado limite de flambagem lateral com torção.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mcrx_FLT(self, cb, lb):
        """
        Retorna o momento fletor em X correspondente a flambagem elástica,
        para o estado limite de flambagem lateral com torção.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

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

        raise NotImplementedError

    # Estado Limite FLM
    def par_esbeltez_limite_Mrdx_FLM(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem local da
        mesa para barras fletidas em relação ao eixo x conforme o que consta na
        seção G2 do anexo G da NBR8800:2008.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mrx_FLM(self):
        """
        Retorna o momento fletor em X correspondente ao inicio de escoamento da seção,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

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

        raise NotImplementedError

    # Estado Limite FLA
    def par_esbeltez_limite_Mrdx_FLA(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem local da
        alma para barras fletidas em relação ao eixo x conforme o que consta na
        seção G2 do anexo G da NBR8800:2008.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mrx_FLA(self):
        """
        Retorna o momento fletor em X correspondente ao inicio de escoamento da seção,
        para o estado limite de flambagem local da alma.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mcrx_FLA(self):
        """
        Retorna o momento fletor em X correspondente a flambagem elástica,
        para o estado limite de flambagem local da mesa.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mnx_FLA(self):
        """
        Determina o momento fletor resistente nominal em X de uma barra para
        o estado limite último de flambagem local da alma.

        Return
        ------
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mrdx_NBR8800(self, Lb, gama_a1=1.1, Cb=1):
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

    @property
    def Mply(self):
        """ Momento de plastificação da seção em relação ao eixo Y"""
        return self.Zy * self.material.fy

    # Estado Limite FLT
    def indice_esbeltex_Y(self, Lb):
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

    def par_esbeltez_limite_Mrdy_FLT(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem lateral com
        torção para barras fletidas em relação ao eixo y conforme o que consta na
        seção G2 do anexo G da NBR8800:2008.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.
        raise NotImplementedError

    def Mry_FLT(self):
        """
        Retorna o momento fletor em Y correspondente ao início de escoamento da seção,
        para o estado limite de flambagem lateral com torção.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    def Mcry_FLT(self, Cb):
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
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

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
        Retorna o momento fletor em Y correspondente ao início de escoamento da seção,
        para o estado limite de flambagem local da mesa.
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
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

    # Estado Limite FLA

    def par_esbeltez_limite_Mrdy_FLA(self):
        """
        Retorna os parâmetros de esbeltez limite de plastificação e de
        início de escoamento, respectivamente, relativos a flambagem local
        da alma para barras fletidas em relação ao eixo y conforme o que
        consta na seção G2 do anexo G da NBR8800:2008.
        """
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.

        raise NotImplementedError

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
        # Como o cálculo deste fator é em função do tipo de seção,
        # este método deve ser implementado em cada uma das
        # classes especificas de cada tipo de perfil.
        raise NotImplementedError

    def Mrdy_NBR8800(self, Lb, Cb=1, gama_a1=1.1):
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
