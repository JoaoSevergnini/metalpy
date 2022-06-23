from math import sqrt, pi
from collections import namedtuple
from warnings import warn


class NBR8800:
    """
    Está classe apresenta os **métodos de verificação da capacidade resistênte** de perfis de aço
    fornecidos pela norma brasileira **NBR8800:2008**: `Projeto de estruturas de aço e de estruturas
    mistas de aço e concreto de edifícios`, de acordo com o método dos estados limites últimos **(ELU)**.

    Static method
    -------------
    Ntrd_brt(perfil, gama_a1=1.1)
        Determina a força resistênte de tração de cálculo ao escoamento da seção bruta do perfil.
    Ncrd(perfil, klx, kly, klz, gama_a1=1.1, data='False')
        Determina a força resistênte de compressão de cálculo de uma barra de aço.
    Vrdx(perfil, gama_a1=1.1, data=False)
        Determina a força resistênte de cisalhamento de cálculo na direção X (Horizontal)
    Vrdy(perfil, gama_a1=1.1, data=False)
        Determina a força resistênte de cisalhamento de cálculo na direção Y (Vertical)
    Mrdx(perfil, Lb=None, Cb=1.0, gama_a1=1.1, data=False)
        Determina o momento resistênte de cálculo do perfil em relação ao eixo X (Horizontal)
    Mrdy(perfil, Lb=None, Cb=1.0, gama_a1=1.1, data=False)
         Determina o momento resistênte de cálculo do perfil em relação ao eixo Y (Vertical)
    """
    c_tensao_res = 0.7

    # -> Métodos para determinação da resistência a tração
    # ---------------------------------------------------------

    @staticmethod
    def Ntrd_brt(perfil, gama_a1=1.1, data=False):
        """
        Método que determina a força axial resistênte de cálculo ao escoamento da seção bruta do perfil
        de acordo com a **NBR8800:2008**.

        Ver seção 5.2.2 da NBR8800:2008

        Parameter
        ---------
        perfil: objeto PerfilEstrutural
            perfil estrutural.
            podendo ser um objeto de uma das seguintes classes:
                - PerfilI
                - PerfilILam
                - TuboRet
                - TuboCir
                - Caixao
        gama_a1: float, default=True
                coeficiente de segurança gama_a1
        data: bool, default=False
              Se data=True o método deve retornar os dados utilizados na obtenção de Ntrd.

        Examples
        --------

        Return
        ------
        Ntrd: float
            Força axial resistênte de cálculo ao escoamento da seção bruta
        Ntrd, Ntrd_data: float, objeto Ntrd_data
            Força axial resistênte de cálculo ao escoamento da seção bruta e dados de cálculo.
            Caso data=True
        """

        warn('Item 5.2.8.1:A NBR8800:2008 recomenda que índice de esbeltez de barras tracionadas não '
             'supere o valor de 300')

        Ntrd_dados = namedtuple('Ntrd_dados', 'Ntrk A fy')
        Ntrd = perfil.Afy / gama_a1
        return Ntrd if not data else (Ntrd, Ntrd_dados(perfil.Afy, perfil.A, perfil.mat.fy))

    # -> Métodos para determinação da resistência a compressão
    # ---------------------------------------------------------

    @staticmethod
    def _kc(esb):
        """ Conforme item (c) da seção F.2 do anexo F da NBR8800 """

        kc = 4 / sqrt(esb)
        kc = 0.35 if kc < 0.35 else kc
        kc = 0.76 if kc > 0.76 else kc
        return kc

    # Métodos auxiliares para a determinação do fator Qs de acordo com a seção F.2 anexo F da NBR8800
    @staticmethod
    def _Qs_g3(perfil):
        """ Fator Qs para mesas de perfis do grupo 3 (cantoneiras) """

        elp = 0.45 * perfil.raiz_E_fy
        elr = 0.91 * perfil.raiz_E_fy

        if elp > perfil.esb:
            return 1
        elif elp < perfil.esb_mesa <= elr:
            return 1.340 - 0.76 * perfil.esb * perfil.raiz_fy_E
        else:
            return 0.53 * (perfil.raiz_E_fy / perfil.esb) ** 2

    @staticmethod
    def _Qs_g4(perfil):
        """ Fator Qs para mesas de perfis do grupo 4 (laminados I, H, U e T) """

        elp = 0.56 * perfil.raiz_E_fy
        elr = 1.03 * perfil.raiz_E_fy

        if elp > perfil.esb_mesa:
            return 1
        elif elp < perfil.esb_mesa <= elr:
            return 1.415 - 0.74 * perfil.esb_mesa * perfil.raiz_fy_E
        else:
            return (0.69 / perfil.esb_mesa ** 2) * perfil.raiz_E_fy ** 2

    @staticmethod
    def _Qs_g5(perfil):
        """ Fator Qs para mesas de perfis do grupo 5 (soldados I, U e T) """

        kc = NBR8800._kc(perfil.esb_alma)

        elp = 0.64 * perfil.raiz_E_fy * sqrt(kc)
        elr = 1.17 * perfil.raiz_E_fy * sqrt(kc)

        if perfil.tipo == 'I SOLDADO' and not perfil.simetria_x:
            esb_mesa = max(perfil.esb_mesa_s, perfil.esb_mesa_i)
        else:
            esb_mesa = perfil.esb_mesa

        if elp > esb_mesa:
            return 1
        elif elp < esb_mesa <= elr:
            return 1.415 - 0.65 * (esb_mesa / sqrt(kc)) * perfil.raiz_fy_E
        else:
            return 0.90 * kc * perfil.raiz_E_fy ** 2 / (esb_mesa ** 2)

    @staticmethod
    def _Qs_g6(perfil):
        """ Fator Qs para alma de perfis T"""

        elp = 0.75 * perfil.raiz_E_fy
        elr = 1.03 * perfil.raiz_E_fy

        if elp > perfil.esb_alma:
            return 1
        elif elp < perfil.esb_alma <= elr:
            return 1.908 - 1.22 * perfil.esb_alma * perfil.raiz_fy_E
        else:
            return 0.69 * perfil.raiz_E_fy ** 2 / perfil.esb_alma ** 2

    @staticmethod
    def _bef(b, t, elr, E, fy, ca, chi):
        """
        Método auxiliar para o cálculo da largura efetiva de um elemento AA (apoiado apoiado) em compressão.

        Ver item F.3.2 do anexo F da NBR8880:2008.
        """

        if b / t <= elr:
            return b
        else:
            raiz_E_fcr = sqrt(E / (fy * chi))
            bef = 1.92 * t * raiz_E_fcr * (1 - ca / (b / t) * raiz_E_fcr)
            return min(bef, b)

    @staticmethod
    def _Qa(perfil, chi):
        """
        Fator Qa para perfis com elementos apoiado-poiado (AA) comprimidos.

        Ver seção F.3 do anexo F da NBR8800:2008.
        """
        elr = 1.40 * perfil.raiz_E_fy if perfil.tipo == "TUBO RET" else 1.49 * perfil.raiz_E_fy

        # Para perfis tubo retangulares e caixão
        if perfil.tipo in ('TUBO RET', 'CAIXAO'):
            ca = 0.38 if perfil.tipo == 'TUBO RET' else 0.34
            bef_alma = NBR8800._bef(perfil.hint, perfil.tw, elr, perfil.mat.E, perfil.mat.fy, ca, chi)
            bef_mesa = NBR8800._bef(perfil.bint, perfil.tf, elr, perfil.mat.E, perfil.mat.fy, ca, chi)

            Aef = perfil.A - 2 * (perfil.bint - bef_mesa) * perfil.tf - 2 * (perfil.hint - bef_alma) * perfil.tw
            return Aef / perfil.A

        # Para outros perfis com elementos AA
        else:
            bef_alma = NBR8800._bef(perfil.h, perfil.tw, elr, perfil.mat.E, perfil.mat.fy, 0.34, chi)
            Aef = perfil.A - (perfil.h - bef_alma) * perfil.tw
            return Aef / perfil.A

    @staticmethod
    def _Q(perfil, chi):
        """ Fator de redução da capacidade resistênte a compressão devido a flambagem local do perfil
        (Ver anexo F da NBR8800:2008) """

        if perfil.tipo in ('I LAMINADO', 'U LAMINADO'):
            return NBR8800._Qa(perfil, chi) * NBR8800._Qs_g4(perfil)

        elif perfil.tipo in ('I SOLDADO', 'U SOLDADO'):
            return NBR8800._Qa(perfil, chi) * NBR8800._Qs_g5(perfil)

        elif perfil.tipo == 'T SOLDADO':
            return NBR8800._Qa(perfil, chi) * min(NBR8800._Qs_g4(perfil), NBR8800._Qs_g6(perfil))

        elif perfil.tipo == 'T LAMINADO':
            return NBR8800._Qa(perfil, chi) * min(NBR8800._Qs_g5(perfil), NBR8800._Qs_g6(perfil))

        elif perfil.tipo in ('CAIXAO', 'TUBO RET'):
            return NBR8800._Qa(perfil, chi)

        elif perfil.tipo == 'CANTONEIRA':
            return NBR8800._Qs_g4(perfil)

        elif perfil.tipo == 'TUBO CIR':

            elp = 0.11 * perfil.raiz_E_fy ** 2
            elr = 0.45 * perfil.raiz_E_fy ** 2

            if perfil.esb <= elp:
                return 1
            elif elp < perfil.esb <= elr:
                return 0.038 * perfil.raiz_E_fy ** 2 / perfil.esb + 2 / 3
            else:
                raise ValueError('Item F.4.2: Esbeltez do perfil supera o valor permitido pela NBR8800 para perfis '
                                 'tubulares circulares ')
        else:
            raise NotImplementedError('Cálculo do fator Q não implementado para perfis do tipo {}'.format(perfil.tipo))

    @staticmethod
    def Ncrd(perfil, klx, kly, klz, gama_a1=1.1, data=False):
        """
        Método que determina a força axial de compressão resistênte de cálculo de uma
        barra de aço de acordo com a **NBR8800:2008**.

        Ver seção 5.3 da NBR8800:2008.

        Parameter
        ----------
        perfil: objeto PerfilEstrutural
            perfil estrutural.
            podendo ser um objeto de uma das seguintes classes:
                - PerfilI
                - PerfilILam
                - TuboRet
                - TuboCir
                - Caixao
        klx: float
            comprimento de flambagem por flexão em relação ao eixo x
        kly: float
            comprimento de flambagem por flexão em relação ao eixo Y
        klz: float
            comprimento de flambagem por torção em relação ao eixo
            longitudinal Z
        gama_a1: float, default=1.1
                coeficiente de segurança
        data: bool, default=False
              Se data=True o método deve retornar os dados utilizados na obtenção de Ncrd.
        Examples
        --------

        Return
        ------
        Ncrd: float
            Força axial de compressão resistênte de cálculo.
        Ncrd, Ncrd_dados: float, objeto Ncrd_dados
            Força axial de compressão resistênte de cálculo e dados de cálculo.
            Caso data=True
        """
        if max(perfil.indice_esbeltez(klx, kly)) > 200:
            raise ValueError('Item 5.3.4.1: O índice de esbeltez de uma barra comprimida não deve ser superior a 200')

        Ncrd_dados = namedtuple('Ncrd_dados', 'Ncrk A fy Ne ier1 Chi1 ier2 Chi2 Q')

        # Fator Chi sem a consideração de flambagem local (Q=1)
        # -----------------------------------------------------

        Ne = perfil.par_estabilidade(klx, kly, klz).Ne
        ier1 = sqrt(perfil.Afy / Ne)  # índice de esbeltez reduzido Q=1

        if ier1 <= 1.5:
            Chi1 = 0.658 ** (ier1 ** 2)
        else:
            Chi1 = 0.877 / ier1 ** 2

        Q = NBR8800._Q(perfil, Chi1)

        # Fator Chi com a consideração de flambagem local
        # -----------------------------------------------

        ier2 = sqrt(Q * perfil.Afy / Ne)  # índice de esbeltez reduzido
        if ier2 <= 1.5:
            Chi2 = 0.658 ** (ier2 ** 2)
        else:
            Chi2 = 0.877 / ier2 ** 2

        Ncrk = Chi2 * Q * perfil.Afy
        Ncrd = Ncrk / gama_a1

        return Ncrd if not data else (Ncrd, Ncrd_dados(Ncrk, perfil.A, perfil.mat.fy, Ne, ier1, Chi1, ier2, Chi2, Q))

    # -> Métodos para determinação da resistência ao esforço cortante
    # ---------------------------------------------------------------

    # EM X
    # -------------
    @staticmethod
    def Vrdx(perfil, Lv=None, gama_a1=1.1, data=False):
        """
        Método que determina a força cortante resistente de cálculo do perfil para cargas aplicadas na direção X
        de acordo com a **NBR8800:2008**.

        ver seção 5.4.3 da NBR8800:2008.

        Parameter
        ---------
        perfil: objeto PerfilEstrutural
            perfil estrutural.
            podendo ser um objeto de uma das seguintes classes:
                - PerfilI
                - PerfilILam
                - TuboRet
                - TuboCir
                - Caixao
        gama_a1: float, default=1.1
                coeficiente de segurança gama_a1
        Lv: float, default=None
            distância entre as seções de forças cortantes máxima e nula.
            (só é necessário caso o perfil seja uma instância da classe TuboCir)
        data: bool, default=False
              Se data=True o método deve retornar os dados utilizados na obtenção de Vrdx.
        Examples
        --------

        Return
        ------
        Vrdx: float
            Força cortante resistênte de cálculo na direção x.
        Vrdx, Vrdx_dados: float, Vrdx_dados
            Força cortante resistênte de cálculo na direção x e dados de cálculo.
            Caso data=True
        """

        if perfil.tipo in ('I SOLDADO', 'I LAMINADO', 'U SOLDADO', 'U LAMINADO', 'TUBO RET'):

            Vrdx_dados = namedtuple('Vrdx_dados', 'Vpl kv elp elr')

            kv = 1.2 if perfil.tipo in ('I SOLDADO', 'I LAMINADO', 'U SOLDADO', 'U LAMINADO') else 5

            elp = 1.1 * sqrt(kv) * perfil.raiz_E_fy
            elr = 1.37 * sqrt(kv) * perfil.raiz_E_fy

            Vpl = perfil.Vplx

            if perfil.esb_mesa <= elp:
                Vrdx = Vpl / gama_a1
                return Vrdx if not data else (Vrdx, Vrdx_dados(Vpl, kv, elp, elr))

            elif elp < perfil.esb_mesa <= elr:
                Vrdx = (elp / perfil.esb_mesa) * (perfil.Vplx / gama_a1)
                return Vrdx if not data else (Vrdx, Vrdx_dados(Vpl, kv, elp, elr))

            else:
                Vrdx = 1.24 * (elp / perfil.esb_mesa) ** 2 * (perfil.Vplx / gama_a1)
                return Vrdx if not data else (Vrdx, Vrdx_dados(Vpl, kv, elp, elr))

        elif perfil.tipo == 'TUBO CIR':
            return NBR8800._Vrd_tubo(perfil, Lv, gama_a1, data)

        else:
            return NotImplementedError('Vrdx não implementado para perfis do tipo {}'.format(perfil.tipo))

    # CORTANTE EM Y
    # -------------
    @staticmethod
    def Vrdy(perfil, a=None, Lv=None, gama_a1=1.1, data=False):

        """
        Método que determina a força cortante resistente de cálculo do perfil para cargas aplicadas na direção Y
        de acordo com a **NBR8800:2008**.

        ver seção 5.4.3 da NBR8800:2008.

        Parameter
        ---------
        perfil: objeto PerfilEstrutural
            perfil estrutural.
            podendo ser um objeto de uma das seguintes classes:
                - PerfilI
                - PerfilILam
                - TuboRet
                - TuboCir
                - Caixao
        gama_a1: float, default=1.1
                coeficiente de segurança gama_a1
        a: float, default=None
            distância entre enrijecedores.
            (só é necessário caso o perfil seja uma instância das classes PerfilI, PerfilILam)
        Lv: float, default=None
            distância entre as seções de forças cortantes máxima e nula.
            (só é necessário caso o perfil seja uma instância da classe TuboCir)
        data: bool, default=False
              Se data=True o método deve retornar os dados utilizados na obtenção de Vrdy.
        Examples
        --------

        Return
        ------
        Vrdy: float
            Força cortante resistênte de cálculo na direção y.
        Vrdy, Vrdy_dados: float, Vrdy_dados
            Força cortante resistênte de cálculo na direção y e dados de cálculo.
            Caso data=True
        """

        if perfil.tipo != 'TUBO CIR':

            Vrdy_dados = namedtuple('Vrdy_dados', 'Vpl kv elp elr')

            if perfil.tipo in ('I SOLDADO', 'I LAMINADO', 'U SOLDADO', 'U LAMINADO'):

                if a is None or a / perfil.h > 3 or a / perfil.h > (260 / perfil.esb_alma) ** 2:
                    kv = 5
                else:
                    kv = 5 + 5 / (a / perfil.h) ** 2

            elif perfil.tipo == 'TUBO RET':
                kv = 5

            else:
                kv = 1.2

            elp = 1.1 * sqrt(kv) * perfil.raiz_E_fy  # Parâmetro de esbeltez limite de plastificação
            elr = 1.37 * sqrt(kv) * perfil.raiz_E_fy  # Parâmetro de esbeltez limite de início de escoamento

            Vpl = perfil.Vply

            if perfil.esb_alma <= elp:
                Vrdy = Vpl / gama_a1
                return Vrdy if not data else (Vrdy, Vrdy_dados(Vpl, kv, elp, elr))

            elif elp < perfil.esb_mesa <= elr:
                Vrdy = (elp / perfil.esb_mesa) * (Vpl / gama_a1)
                return Vrdy if not data else (Vrdy, Vrdy_dados(Vpl, kv, elp, elr))

            else:
                Vrdy = 1.24 * (elp / perfil.esb_mesa) ** 2 * (Vpl / gama_a1)
                return Vrdy if not data else (Vrdy, Vrdy_dados(Vpl, kv, elp, elr))

        elif perfil.tipo == 'TUBO CIR':
            return NBR8800._Vrd_tubo(perfil, Lv, gama_a1, data)

        else:
            return NotImplementedError('Vrdy não implementado para perfis do tipo {}'.format(perfil.tipo))

    @staticmethod
    def _Vrd_tubo(perfil, Lv, gama_a1, data):
        """ Determina a força cortante resistênte de cálculo para tubos circulares"""

        if Lv is None:
            raise ValueError('Lv não fornecido')

        Vrdy_dados = namedtuple('Vrdy_dados', 'fcr')

        fcr1 = 1.60 * perfil.mat.E / (sqrt(Lv / perfil.D) * perfil.esb ** (5 / 4))
        fcr2 = 0.78 * perfil.mat.E / perfil.esb ** 1.5

        fcr = max(fcr1, fcr2) if max(fcr1, fcr2) < 0.60 * perfil.mat.fy else 0.6 * perfil.mat.fy

        Vrdy = fcr * perfil.Awy / gama_a1
        return Vrdy if not data else (Vrdy, Vrdy_dados(fcr))

    # MOMENTO FLETOR EM X
    # ------------
    # Estado Limite FLT
    @staticmethod
    def _Mnx_FLT(perfil, Cb, Lb):

        """ Determina o momento fletor resistente nominal de uma barra para o estado limite último
            de flambagem lateral com torção em relação ao eixo X"""

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'TUBO RET', 'CAIXAO', 'U SOLDADO', 'U LAMINADO'):

            esb = perfil.indice_esbeltez(Lb, Lb)[0]

            # Determinação dos parâmetros necessários para determinação do momento fletor
            if perfil.tipo in ('TUBO RET', 'CAIXAO'):
                elp = 0.13 * perfil.mat.E * sqrt(perfil.J * perfil.A) / perfil.Mplx
                elr = 0.13 * perfil.mat.E * sqrt(perfil.J * perfil.A) / perfil.Mplx

            else:
                elp = 1.76 * perfil.raiz_E_fy

                beta_1 = perfil.Wxs * NBR8800.c_tensao_res / (perfil.mat.E * perfil.J)
                beta_2 = 1

                if perfil.tipo == 'I SOLDADO' and not perfil.bissimetrico:
                    alfa_y = perfil.Iys / perfil.Iyi
                    beta_3 = 0.45 * (perfil.d - (perfil.tfs + perfil.tfi) / 2) * (alfa_y - 1) / (alfa_y + 1)
                    beta_2 = 5.2 * beta_1 * beta_3 + 1

                # primeiro termo da equação
                t1 = 1.38 * sqrt(perfil.Iy * perfil.J) / (perfil.rx + perfil.J + beta_1)
                # segundo termo da equação
                t2 = sqrt(beta_2 + sqrt(beta_2 ** 2 + 27 * perfil.Cw * beta_1) / perfil.Iy)
                elr = t1 * t2

            if esb < elp:
                return perfil.Mplx

            elif elp < esb < elr:
                Mr = perfil.Mrx * 0.7 if not perfil.tipo == 'I LAMINADO' and not perfil.bissimetrico \
                    else min(0.7 * perfil.Wxs, 0.7 * perfil.Mr)
                return perfil.Mplx - (perfil.Mplx - Mr) * (esb - elp) / (elr - elp)

            elif elr < esb:
                Me = perfil.par_estabilidade(Lb, Lb, Lb).Me
                return Cb * Me if Me < perfil.Mplx else perfil.Mplx

        elif perfil.tipo in ('T LAMINADO', 'T SOLDADO'):
            B = 2.3 * perfil.d * sqrt(perfil.Iy / perfil.J) / Lb
            return pi * sqrt(perfil.EIy * perfil.GJ) * (B + sqrt(1 + B ** 2)) / Lb

    # Estado Limite FLM
    @staticmethod
    def _Mnx_FLM(perfil):
        """
        Determina o momento fletor resistente nominal em X de uma barra para
        o estado limite último de flambagem local da mesa.
        """

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'TUBO RET', 'CAIXAO', 'U SOLDADO', 'U LAMINADO'):

            # Determinação dos parâmetros necessários para determinação do momento fletor

            if perfil.tipo in ('TUBO RET', 'CAIXAO'):
                elp = 1.12 * perfil.raiz_E_fy
                elr = 1.40 * perfil.raiz_E_fy

                Mrx = perfil.mat.fy  # calcular Wef
                Mcrx = perfil.Wx * perfil.mat.fy  # calcular Wef e fazer Wef/Wx * fy

            else:
                elp = 0.38 * perfil.raiz_E_fy

                if perfil.tipo in ('I SOLDADO', 'U SOLDADO'):
                    kc = 4 / sqrt(perfil.esb_alma)
                    elr = 0.95 * sqrt(kc / 0.7) * perfil.raiz_E_fy
                    Mcrx = 0.90 * perfil.mat.E * kc * perfil.Wxs / perfil.esb_mesa ** 2
                else:
                    elr = 0.83 * sqrt(1 / 0.7) * perfil.raiz_E_fy
                    Mcrx = 0.69 * perfil.mat.E * perfil.Wx / perfil.esb_mesa ** 2

                Mrx = 0.7 * perfil.Wxs

            if perfil.esb_mesa < elp:
                return perfil.Mplx
            elif elp < perfil.esb_mesa < elr:
                return perfil.Mplx - (perfil.Mplx - Mrx) * (perfil.esb_alma - elp) / (elr - elp)
            elif perfil.esb_mesa > elr:
                return Mcrx

        elif perfil.tipo == ('T SOLDADO', 'T LAMINADO'):

            elp = 0.38 * perfil.raiz_E_fy
            elr = 1 * perfil.raiz_E_fy

            if perfil.esb_mesa < elp:
                return perfil.Mplx
            elif elp < perfil.esb_mesa < elr:
                return (1.19 - 0.5 * perfil.esb_mesa * perfil.raiz_fy_E) * perfil.Ws * perfil.mat.fy
            elif perfil.esb_mesa > elr:
                return 0.69 * perfil.Ws * perfil.mat.fy / perfil.esb_mesa ** 2

    # Estado Limite FLA

    @staticmethod
    def _Mnx_FLA(perfil):
        """
        Determina o momento fletor resistente nominal em X de uma barra para
        o estado limite último de flambagem local da alma.
        """

        elr = 5.7 * perfil.raiz_E_fy

        if perfil.tipo in ('I LAMINADO', 'U SOLDADO', 'U LAMINADO', 'CAIXAO') or (perfil.tipo == 'I SOLDADO'
                                                                                  and perfil.bissimetrico):
            elp = 3.76 * perfil.raiz_E_fy

        elif perfil.tipo == 'I SOLDADO' and not perfil.bissimetrico:

            hc = 2 * (perfil.d - perfil.tfs - perfil.hcg)
            hp = 2 * (perfil.d - perfil.tfs - perfil.hpl)

            elp = (hc / hp) * perfil.raiz_E_fy / (0.54 * perfil.Mpl / perfil.Mrx - 0.09) ** 2

        else:
            elp = 2.42 * perfil.raiz_E_fy

        if perfil.esb_alma < elp:
            return perfil.Mplx
        elif elp < perfil.esb_alma < elr:
            return perfil.Mplx - (perfil.Mplx - perfil.Mrx) * (perfil.esb_alma - elp) / (elr - elp)
        else:
            print("ALMA ESBELTA")

    @staticmethod
    def _Mn_Tubo(perfil, data=False):

        Mrd_dados = namedtuple("Mrd_dados", "Mn elp elr")

        elp = 0.07 * perfil.mat.E / perfil.mat.fy
        elr = 0.31 * perfil.mat.E / perfil.mat.fy

        if perfil.esb < elp:
            Mn = perfil.Mplx
            return Mn if not data else (Mn, Mrd_dados(Mn, elp, elr))

        elif elp < perfil.esb <= elr:
            Mn = (0.021 * perfil.mat.E / perfil.esb + perfil.mat.fy) * perfil.W
            return Mn if not data else (Mn, Mrd_dados(Mn, elp, elr))

        else:
            Mn = 0.33 * perfil.mat.E * perfil.Wx / perfil.esb
            return Mn if not data else (Mn, Mrd_dados(Mn, elp, elr))

    @staticmethod
    def Mrdx(perfil, Lb=None, gama_a1=1.1, Cb=1):
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

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO', 'U SOLDADO'):
            return min(NBR8800._Mnx_FLT(perfil, Lb, Cb), NBR8800._Mnx_FLM(perfil), NBR8800._Mnx_FLA(perfil),
                       1.5 * perfil.Mrx) / gama_a1

        elif perfil.tipo in ('CAIXAO', 'TUBO RET') and perfil.Wx > perfil.Wy:
            return min(NBR8800._Mnx_FLT(perfil, Lb, Cb), NBR8800._Mnx_FLM(perfil), NBR8800._Mnx_FLA(perfil),
                       1.5 * perfil.Mrx) / gama_a1

        elif perfil.tipo in ('CAIXAO', 'TUBO RET') and perfil.Wx <= perfil.Wy:
            return min(NBR8800._Mnx_FLM(perfil), NBR8800._Mnx_FLA(perfil), 1.5 * perfil.Mrx) / gama_a1

        elif perfil.tipo in ('T LAMINADO', 'T SOLDADO'):
            return min(NBR8800._Mnx_FLT(perfil, Cb, Lb), NBR8800._Mnx_FLA(perfil), 1.5 * perfil.Mrx) / gama_a1

        elif perfil.tipo == 'TUBO CIR':
            return min(NBR8800._Mn_Tubo(perfil), 1.5 * perfil.Mrx) / gama_a1

    # MOMENTO FLETOR EM Y
    # ------------

    # Estado Limite FLT
    @staticmethod
    def _Mny_FLT(perfil, Cb, Lb):

        elp = 0.13 * perfil.mat.E * sqrt(perfil.J * perfil.A) / perfil.Mpl
        elr = 0.13 * perfil.mat.E * sqrt(perfil.J * perfil.A) / perfil.Mpl

        esb = perfil.indice_esbeltez(Lb, Lb)[1]

        if esb < elp:
            return perfil.Mply
        elif elp < esb < elr:
            return perfil.Mply - (perfil.Mply - 0.7 * perfil.Mry) * (esb - elp) / (elr - elp)
        else:
            Me = perfil.par_estabilidade(Lb, Lb, Lb).Me
            return Cb * Me if Me < perfil.Mplx else perfil.Mply

    @staticmethod
    def _Mny_FLM(perfil):

        # Determinação dos parâmetros necessários para determinação do momento fletor

        if perfil.tipo in ('TUBO RET', 'CAIXAO'):
            elp = 1.12 * perfil.raiz_E_fy
            elr = 1.40 * perfil.raiz_E_fy

            Mry = Wef * perfil.mat.fy  # calcular Wef
            Mcry = perfil.Wy * perfil.mat.fy  # calcular Wef e fazer Wef/Wy * fy

        else:
            elp = 0.38 * perfil.raiz_E_fy

            if perfil.tipo in ('I SOLDADO', 'U SOLDADO'):
                kc = 4 / sqrt(perfil.esb_alma)
                elr = 0.95 * sqrt(kc / 0.7) * perfil.raiz_E_fy
                Mcry = 0.90 * perfil.mat.E * kc * perfil.Wys / perfil.esb_mesa ** 2

            else:
                elr = 0.83 * sqrt(1 / 0.7) * perfil.raiz_E_fy
                Mcry = 0.69 * perfil.mat.E * perfil.Wys / perfil.esb_mesa ** 2

            Mry = 0.7 * perfil.Wy

        if perfil.esb_mesa < elp:
            return perfil.Mply
        elif elp < perfil.esb_mesa < elr:
            return perfil.Mply - (perfil.Mply - Mry) * (perfil.esb_alma - elp) / (elr - elp)
        elif perfil.esb_mesa > elr:
            return Mcry

    @staticmethod
    def _Mny_FLA(perfil):

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO', 'U SOLDADO'):

            elp = 1.12 * perfil.raiz_E_fy
            elr = 1.40 * perfil.raiz_E_fy

            Wef = perfil.W  # IMPLEMENTAR Wef

            Mry = Wef * perfil.mat.fy
            Mcry = Wef ** 2 * perfil.mat.fy / perfil.Wy

        else:

            elp = 3.72 * perfil.raiz_E_fy if perfil.tipo == 'CAIXAO' else 2.42 * perfil.raiz_E_fy
            elr = 5.7 * perfil.raiz_E_fy

            Mry = perfil.Wy
            Mcry = 'Perfis tubulares não apresentam'

            # PARA PERFIS TUBULARES NÂO EXISTE O Mcry PARA FLAMBAGEM LOCAL DA ALMA??

        if perfil.esb_mesa < elp:
            return perfil.Mply
        elif elp < perfil.esb_mesa < elr:
            return perfil.Mply - (perfil.Mply - Mry) * (perfil.esb_alma - elp) / (elr - elp)
        elif perfil.esb_mesa > elr:
            return Mcry

    @staticmethod
    def Mrdy(perfil, Lb, gama_a1=1.1, Cb=1):

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO', 'U LAMINADO') and perfil.bissimetrico:
            return min(NBR8800._Mny_FLM(perfil), NBR8800._Mny_FLA(perfil)) / gama_a1

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Iy > perfil.Ix:
            return min(NBR8800._Mny_FLT(perfil, Cb, Lb), NBR8800._Mny_FLM(perfil), NBR8800._Mny_FLA(perfil)) / gama_a1

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Iy <= perfil.Ix:
            return min(NBR8800._Mny_FLM(perfil), NBR8800._Mny_FLA(perfil)) / gama_a1

        elif perfil.tipo == 'TUBO CIR':
            return NBR8800._Mn_Tubo(perfil) / gama_a1

        else:
            print('Perfil não apresenta resistência ao momento em relação ao momento em torno do eixo Y')
