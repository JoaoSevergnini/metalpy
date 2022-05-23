from perfis import *


class Nbr8800:
    c_tensao_res = 0.7

    # TRAÇÂO
    # --------

    @staticmethod
    def Ntrd_ESB(perfil, gama_a1=1.1):
        """
        Método que determina a resistência ao escoamento da seção bruta de perfil métálico

        ver seção 5.2.2 da NBR8800:2008

        Parameter
        ---------
        gama_a1: 'float' (default=1,1)
                coeficiente de segurança gama_a1

        Return
        ------

        """
        return perfil.Afy / gama_a1

    # COMPRESSÃO
    # -----------

    @staticmethod
    def _Qs_g3(perfil):
        pass

    @staticmethod
    def _Qs_g4(perfil):
        """ Fator Qs para mesas de perfis laminados do tipo I, U e T """

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
        """ Fator Qs para mesas de perfis soldados do tipo I, U e T"""

        kc = 4 / sqrt(perfil.esb_alma)
        kc = 0.35 if kc < 0.35 else kc
        kc = 0.76 if kc > 0.76 else kc

        elp = 0.64 * perfil.raiz_E_fy * sqrt(kc)
        elr = 1.17 * perfil.raiz_E_fy * sqrt(kc)

        if perfil.tipo == 'I SOLDADO' and not perfil.simetria_x:
            esb_mesa = max(perfil.esb_mesa_s, perfil.esb_mesa_sup)
        else:
            esb_mesa = perfil.esb_mesa

        if elp > esb_mesa:
            return 1
        elif elp < esb_mesa <= elr:
            return 1.415 - 0.65 * esb_mesa / sqrt(kc) * perfil.raiz_fy_E
        else:
            return 0.90 * kc / esb_mesa ** 2 * perfil.raiz_E_fy ** 2

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
    def _Qa(perfil, chi):
        """ Fator Qa para perfis com elementos apoiado-poiado"""
        elr = 1.40 * perfil.raiz_E_fy if perfil.tipo == "TUBO RET" else 1.49 * perfil.raiz_E_fy

        if perfil.tipo in ('TUBO RET', 'CAIXAO'):

            if perfil.esb_alma < elr:
                bef_alma = perfil.hint
            else:
                ca = 0.38 if perfil.tipo == 'TUBO RET' else 0.34
                raiz_E_fcr = perfil.raiz_E_fy / sqrt(chi)
                bef = 1.92 * perfil.tw * raiz_E_fcr * (1 - ca * raiz_E_fcr / perfil.esb_alma)
                bef_alma = bef if bef < perfil.hint else perfil.hint

            if perfil.esb_mesa < elr:
                bef_mesa = perfil.bint

            else:
                ca = 0.38 if perfil.tipo == 'TUBO RET' else 0.34
                raiz_E_fcr = perfil.raiz_E_fy / sqrt(chi)
                bef = 1.92 * perfil.tf * raiz_E_fcr * (1 - ca * raiz_E_fcr / perfil.esb_mesa)
                bef_mesa = bef if bef < perfil.bint else perfil.bint

            return 2 * (bef_mesa * bef_alma) / perfil.A

        else:
            if perfil.esb_alma < elr:
                bef_alma = perfil.hint
            else:
                raiz_E_fcr = perfil.raiz_E_fy / sqrt(chi)
                bef = 1.92 * perfil.tw * raiz_E_fcr * (1 - 0.34 * raiz_E_fcr / perfil.esb_alma)
                bef_alma = bef if bef < perfil.h else perfil.h
            return 1 - bef_alma * perfil.tw / perfil.A

    @staticmethod
    def _Q(perfil, chi):
        """ Fator de redução da capacidade resistênte a compressão devido a flambagem local do perfil
        (Ver anexo F da NBR8800:2008) """

        if perfil.tipo in ('I LAMINADO', 'U LAMINADO'):
            return Nbr8800._Qa(perfil, chi) * Nbr8800._Qs_g4(perfil)

        elif perfil.tipo in ('I SOLDADO', 'U SOLDADO'):
            return Nbr8800._Qa(perfil, chi) * Nbr8800._Qs_g5(perfil)

        elif perfil.tipo == 'T SOLDADO':
            return Nbr8800._Qa(perfil, chi) * min(Nbr8800._Qs_g4(perfil), Nbr8800._Qs_g6(perfil))

        elif perfil.tipo == 'T LAMINADO':
            return Nbr8800._Qa(perfil, chi) * min(Nbr8800._Qs_g5(perfil), Nbr8800._Qs_g6(perfil))

        elif perfil.tipo in ('CAIXAO', 'TUBO RET'):
            return Nbr8800._Qa(perfil, chi)

        elif perfil.tipo == 'TUBO CIR':

            elp = 0.11 * perfil.raiz_E_fy ** 2
            elr = 0.45 * perfil.raiz_E_fy ** 2

            if perfil.esb <= elp:
                return 1
            elif elp < perfil.esb <= elr:
                return 0.038 * perfil.raiz_E_fy ** 2 / perfil.esb + 2 / 3

    @staticmethod
    def Ncrd(perfil, klx, kly, klz, gama_a1=1.1, data=False):
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

        # Fator Chi sem a consideração de flambagem local
        # -----------------------------------------------

        Ne = perfil.par_estabilidade(klx, kly, klz)['Ne']
        ier = sqrt(perfil.Afy / Ne)

        if ier <= 1.5:
            Chi = 0.658 ** (ier ** 2)
        else:
            Chi = 0.877 / ier ** 2

        Q = Nbr8800._Q(perfil, Chi)

        # Fator Chi sem a consideração de flambagem local
        # -----------------------------------------------

        Ne = perfil.par_estabilidade(klx, kly, klz)['Ne']
        ier = sqrt(Q * perfil.Afy / Ne)
        if ier <= 1.5:
            Chi = 0.658 ** (ier ** 2)
        else:
            Chi = 0.877 / ier ** 2

        Ncrd = Chi * Q * perfil.Afy / gama_a1

        return Ncrd if not data else Ncrd, {'Ne': Ne, 'ier': ier, 'Chi': Chi, 'Q': Q}

    # CORTANTE
    # -----------

    # EM X
    # -------------
    @staticmethod
    def Vrdx(perfil, gama_a1=1.1, a=None, data=False):
        """
        Método que determina a força cortante resistente de cálculo na
        direção X do perfil de acordo com a NBR8800:2008.

        ver a seção 5.4.3 da NBR8800:2008.

        Parameter
        --------
        perfil: 'Perfil'

        a: 'float'
            distância entre eixos de enrijecedores

        gama_a1: 'float' (default = 1.1)
            coeficiente de minoração da resistência.

        Return
        ------
        Vrdx: 'float'
            Força cortante resistênte de cálculo na direção x.
        """

        if perfil.tipo in ('I SOLDADO', 'I LAMINADO', 'U SOLDADO', 'U LAMINADO', 'TUBO RET'):

            kv = 1.2 if perfil.tipo in ('I SOLDADO', 'I LAMINADO', 'U SOLDADO', 'U LAMINADO') else 5

            elp = 1.1 * sqrt(kv) * perfil.raiz_E_fy
            elr = 1.37 * sqrt(kv) * perfil.raiz_E_fy

            if perfil.esb_mesa <= elp:
                Vrdx = perfil.Vplx / gama_a1
                return Vrdx if not data else Vrdx, {'kv': kv, 'elp': elp, 'elr': elr}

            elif elp < perfil.esb_mesa <= elr:
                Vrdx = (elp / perfil.esb_mesa) * (perfil.Vplx / gama_a1)
                return Vrdx if not data else Vrdx, {'kv': kv, 'elp': elp, 'elr': elr}

            else:
                Vrdx = 1.24 * (elp / perfil.esb_mesa) ** 2 * (perfil.Vplx / gama_a1)
                return Vrdx if not data else Vrdx, {'kv': kv, 'elp': elp, 'elr': elr}

        elif perfil.tipo == 'TUBO CIR':

            fcr1 = 1.60 * perfil.material.E / (sqrt(a / perfil.D) * perfil.esb ** (5/4))
            fcr2 = 0.78 * perfil.material.E / perfil.esb ** 1.5

            fcr = max(fcr1, fcr2) if max(fcr1, fcr2) < 0.60 * perfil.material.fy else 0.6 * perfil.material.fy

            return fcr * perfil.Awy / gama_a1

        else:
            return None

    # CORTANTE EM Y
    # -------------
    @staticmethod
    def Vrdy(perfil, a=None, gama_a1=1.1, data=False):

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

        if perfil.tipo != 'TUBO CIR':

            if perfil.tipo in ('I SOLDADO', 'I LAMINADO', 'U SOLDADO', 'U LAMINADO'):

                if a is None or a / perfil.h > 3 or a / perfil.h > (260 / perfil.esb_alma) ** 2:
                    kv = 5
                else:
                    kv = 5 + 5 / (a / perfil.h) ** 2

            elif perfil.tipo == 'TUBO RET':
                kv = 5

            elif perfil.tipo in ('T LAMINADO', 'T SOLDADO'):
                kv = 1.2

            elp = 1.1 * sqrt(kv) * perfil.raiz_E_fy
            elr = 1.37 * sqrt(kv) * perfil.raiz_E_fy

            if perfil.esb_alma <= elp:
                Vrdy = perfil.Vply / gama_a1
                return Vrdy if not data else Vrdy, {'kv': kv, 'elp': elp, 'elr': elr}

            elif elp < perfil.esb_mesa <= elr:
                Vrdy = (elp / perfil.esb_mesa) * (perfil.Vplx / gama_a1)
                return Vrdy if not data else Vrdy, {'kv': kv, 'elp': elp, 'elr': elr}

            else:
                Vrdy = 1.24 * (elp / perfil.esb_mesa) ** 2 * (perfil.Vplx / gama_a1)
                return Vrdy if not data else Vrdy, {'kv': kv, 'elp': elp, 'elr': elr}

        else:
            fcr1 = 1.60 * perfil.material.E / (sqrt(a / perfil.D) * perfil.esb ** (5/4))
            fcr2 = 0.78 * perfil.material.E / perfil.esb ** 1.5

            fcr = max(fcr1, fcr2) if max(fcr1, fcr2) < 0.60 * perfil.material.fy else 0.6 * perfil.material.fy

            return fcr * perfil.Awy / gama_a1

    # MOMENTO FLETOR EM X
    # ------------
    # Estado Limite FLT
    def Mrx_FLT(self):
        """
        Retorna o momento fletor em X correspondente ao início de escoamento da seção,
        para o estado limite de flambagem lateral com torção.
        """
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
        return Cb * self.Mex(Lb, Lb)

    @staticmethod
    def _Mnx_FLT(perfil, Cb, Lb):
        """ Determina o momento fletor resistente nominal de uma barra para o estado limite último
            de flambagem lateral com torção em relação ao eixo X"""

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'TUBO RET', 'CAIXAO', 'U SOLDADO', 'U LAMINADO'):

            esb = perfil.indice_esbeltez(Lb, Lb)[1]

            if perfil.tipo in ('TUBO RET', 'CAIXAO'):
                elp = 0.13 * perfil.material.E * sqrt(perfil.J * perfil.A) / perfil.Mpl
            else:
                elp = 1.76 * perfil.raiz_E_fy

            if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U SOLDADO', 'U LAMINADO'):

                beta_1 = perfil.Wxs * Nbr8800.c_tensao_res / (perfil.material.E * perfil.J)
                beta_2 = 1

                if perfil.tipo('I LAMINADO') and not perfil.bi_simetrica:
                    alfa_y = None
                    beta_3 = 0.45 * (perfil.d - (perfil.tfs + perfil.tfi)/2) * (alfa_y - 1) / (alfa_y + 1)
                    beta_2 = 5.2 * beta_1 * beta_3 + 1

                elr = 1.38 * sqrt(perfil.Iy * perfil.J) / perfil.rx + perfil.J




                if esb < elp:
                    return perfil.Mplx
                elif elp < esb < elr:
                    return perfil.Mplx - (perfil.Mplx - perfil.Mrx_FLT) * (esb - elp) / (elr - elp)
                elif elr < esb:
                    Mcrx = perfil.Mcrx_FLT(Cb, Lb)
                    return Mcrx if Mcrx < perfil.Mplx else perfil.Mplx

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
        # Este método está implementado para os perfis apresentados
        # na tabela G1 da do anexo G da NBR8800 para os perfis que
        # não estão contidos na tabela esse método deve ser sobrescrito,
        # caso o perfil não apresente FLM como um estado limite esse método
        # deve ser sobrescrito retornando o momento de plastificação (Mpl)

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

        # Este método está implementado para os perfis apresentados
        # na tabela G1 da do anexo G da NBR8800 para os perfis que
        # não estão contidos na tabela esse método deve ser sobrescrito,
        # caso o perfil não apresente FLA como um estado limite esse método
        # deve ser sobrescrito retornando o momento de plastificação (Mpl)

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
