from math import sqrt, pi
from collections import namedtuple
from warnings import warn


class AISC360:

    """
    Está classe apresenta os **métodos de verificação da capacidade resistênte** de perfis de aço
    fornecidos pela norma americana **AISC360**: `Specification for structural steel buildings`, de
    acordo com o método dos estados limites últimos **(ELU)**.

    Static method
    -------------
    Ntrd_brt(perfil, phi_s=0.9)
        Determina a força resistênte de tração de cálculo ao escoamento da seção bruta do perfil.
    Ncrd(perfil, klx, kly, klz, phi_c=0.9, data=False)
        Determina a força resistênte de compressão de cálculo de uma barra de aço.
    Vrdx(perfil, phi_v=0.90, data=False)
        Determina a força resistênte de cisalhamento de cálculo na direção X (Horizontal)
    Vrdy(perfil, a=None, phi_v=0.90, data=False)
        Determina a força resistênte de cisalhamento de cálculo na direção Y (Vertical)
    Mrdx(perfil, Lb, Cb, theta_b=0.90, data=False)
        Determina o momento resistênte de cálculo do perfil em relação ao eixo X (Horizontal)
    Mrdy(perfil, Lb, Cb, theta_b=0.90, data=False)
         Determina o momento resistênte de cálculo do perfil em relação ao eixo Y (Vertical)
    """

    # -> Métodos para determinação da resistência a tração
    # ----------------------------------------------------

    @staticmethod
    def Ntrd_brt(perfil, phi_s=0.9, data=False):
        """
        Método que determina a força axial resistênte de cálculo ao escoamento da seção bruta do perfil
        de acordo com a **AISC360-16**.

        ver seção D2.(a) do capítulo D da AISC360-16

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
        phi_s: float, default=True
                coeficiente de segurança gama_a1
        data: bool, default=False
              Se data=True o método deve retornar os dados utilizados na obtenção de Ntrd.
        Examples
        --------

        Returns
        -------
        Ntrd: float
            Força axial resistênte de cálculo ao escoamento da seção bruta
        Ntrd, Ntrd_data: float, objeto Ntrd_data
            Força axial resistênte de cálculo ao escoamento da seção bruta e dados de cálculo.
            Caso data=True

        """

        warn('Item D1: A AISC360-16 recomenda que índice de esbeltez de barras tracionadas não '
             'supere o valor de 300')

        Ntrd_dados = namedtuple('Ntrd_dados', 'Ntrk A fy')
        Ntrd = perfil.Afy * phi_s
        return Ntrd if not data else (Ntrd, Ntrd_dados(perfil.Afy, perfil.A, perfil.mat.fy))

    # -> Métodos para determinação da resistência a compressão
    # ---------------------------------------------------------

    @staticmethod
    def _bef(b, c1, elr, esb, Fy, Fcr):
        """
        Método auxiliar para o cálculo da largura efetiva de um elemento em compressão.

        Ver item seção E7 do capítulo E da AISC360:2008.
        """

        if esb <= elr * sqrt(Fy/Fcr):
            return b
        else:
            c2 = (1 - sqrt(1 - 4 * c1)) / (2 * c1)
            Fel = (c2 * elr / esb) ** 2 * Fy
            bef = b * (1 - c1 * sqrt(Fel / Fcr)) * sqrt(Fel / Fcr)
            return min(bef, b)

    @staticmethod
    def _kc(esb):
        """ Conforme tabela B4.1a do capítulo B da AISC360"""

        kc = 4 / sqrt(esb)
        kc = 0.35 if kc < 0.35 else kc
        kc = 0.76 if kc > 0.76 else kc
        return kc

    @staticmethod
    def _Aef(perfil, Fcr):
        """ Método para a determinação a área efetiva de um perfil em compressão (Ver capítulo E da AISC360-16)"""

        Fy = perfil.mat.fy

        if perfil.tipo in ('I LAMINADO', 'U LAMINADO', 'T LAMINADO'):

            # esbeltez limite da mesa
            elr_ms = 0.56 * perfil.raiz_E_fy

            # esbeltez limite da alma
            elr_alm = 0.75 * perfil.raiz_E_fy if perfil.tipo == 'T LAMINADO' else 1.49 * perfil.raiz_E_fy

            # largura efetiva das mesas
            b = perfil.bf if perfil.tipo == 'U LAMINADO' else perfil.bf / 2
            bef_ms = AISC360._bef(b, 0.22, elr_ms, perfil.esb_mesa, Fy, Fcr)

            # largura efetiva da alma
            c1 = 0.22 if perfil.tipo == 'T LAMINADO' else 0.18
            h = perfil.d if perfil.tipo == 'T LAMINADO' else perfil.h
            bef_al = AISC360._bef(h, c1, elr_alm, perfil.esb_alma, Fy, Fcr)

            # Número de 'abas'
            n = 4 if perfil.tipo == 'I LAMINADO' else 2

            return perfil.A - n * (b - bef_ms) * perfil.tf - (h - bef_al) * perfil.tw

        elif perfil.tipo == 'I SOLDADO':

            kc = AISC360._kc(perfil.esb_alma)
            elr_ms = 0.64 * sqrt(kc) * perfil.raiz_E_fy

            elr_alm = 1.49 * perfil.raiz_E_fy

            # Largura efetiva das mesas
            bef_mss = AISC360._bef(perfil.bfs / 2, 0.22, elr_ms, perfil.esb_mesa_s, Fy, Fcr)
            bef_msi = AISC360._bef(perfil.bfi / 2, 0.22, elr_ms, perfil.esb_mesa_i, Fy, Fcr)

            # Largura efetiva da alma
            bef_alm = AISC360._bef(perfil.dl, 0.18, elr_alm, perfil.esb_alma, Fy, Fcr)

            Aef = perfil.A - (perfil.bfs - 2 * bef_mss) * perfil.tfs  - (perfil.bfi - 2 * bef_msi) * perfil.tfi \
                - (perfil.dl - bef_alm) * perfil.tw

            return Aef

        elif perfil.tipo in ('CAIXAO', 'TUBO RET'):

            elr = 1.40 * perfil.raiz_E_fy if perfil.tipo == 'TUBO RET' else 1.49 * perfil.raiz_E_fy

            bef_ms = AISC360._bef(perfil.bint, 0.2, elr, perfil.esb_mesa, Fy, Fcr)
            bef_alm = AISC360._bef(perfil.hint, 0.2, elr, perfil.esb_alma, Fy, Fcr)

            return perfil.A - 2 * (perfil.bint - bef_ms) * perfil.tf - 2 * (perfil.hint - bef_alm) * perfil.tw

        elif perfil.tipo == 'TUBO CIR':

            elp = 0.11 * perfil.raiz_E_fy ** 2
            elr = 0.45 * perfil.raiz_E_fy ** 2

            if perfil.esb <= elp:
                return perfil.A
            elif elp < perfil.esb <= elr:
                return (0.038 * perfil.raiz_E_fy ** 2 / perfil.esb + 2 / 3) * perfil.A
            else:
                raise ValueError('A AISC360 não prevê o uso de perfis tubulares com esbeltez maior do que 0.45*E/fy')

        else:
            raise NotImplementedError('Cálculo da Aef não implementado para perfil do tipo {}'.format(perfil.tipo))

    @staticmethod
    def Ncrd(perfil, klx, kly, klz, phi_c=0.9, data=False):
        """
        Método que determina a força axial de compressão resistênte de cálculo de uma
        barra de aço de acordo com a **AISC360-16**.

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
        phi_c: float, default=0.9
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
            raise ValueError('Nota item E2: O índice de esbeltez de uma barra comprimida não deve ser superior a 200')

        Ncrd_dados = namedtuple('Ncrd_dados', 'Ncrk A Fy Fe Fy_Fe Fcr Aef')

        Fe = perfil.par_estabilidade(klx, kly, klz).fe  # Tensão crítica de flambagem

        Fy_Fe = perfil.mat.fy / Fe

        if Fy_Fe <= 2.25:
            Fcr = (0.658 ** Fy_Fe) * perfil.mat.fy
        else:
            Fcr = 0.877 * Fe

        Aef = AISC360._Aef(perfil, Fcr)

        Ncrk = Fcr * Aef
        Ncrd = Ncrk * phi_c

        return Ncrd if not data else (Ncrd, Ncrd_dados(Ncrk, perfil.A, perfil.mat.fy, Fe, Fy_Fe, Fcr, Aef))

    # -> Métodos para determinação da resistência ao esforço cortante
    # ---------------------------------------------------------------

    # cortante em X
    # -------------
    @staticmethod
    def Vrdx(perfil, Lv=None, phi_v=0.90, data=False):
        """
        Método que determina a força cortante resistente de cálculo do perfil para cargas aplicadas na direção X
        de acordo com a **AISC360-16**.

        ver capítulo G da AISC360-16.

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
        phi_v: float, default=0.9
                coeficiente de segurança
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
        Vrdx, Vrdx_dados: float, Vrdy_dados
            Força cortante resistênte de cálculo na direção y e dados de cálculo.
            Caso data=True
        """

        if perfil.tipo != 'TUBO CIR':
            Vrdx_dados = namedtuple('Vrdx_dados', 'Vpl Aw kv Cv2 elp elr')

            if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO', 'T LAMINADO'):
                kv = 1.2
            else:
                kv = 5

            elp = 1.1 * sqrt(kv) * perfil.raiz_E_fy  # Parâmetro de esbeltez limite de plastificação
            elr = 1.37 * sqrt(kv) * perfil.raiz_E_fy  # Parâmetro de esbeltez limite de início de escoamento

            if perfil.esb_mesa <= elp:
                Cv2 = 1
            elif elp < perfil.esb_mesa <= elr:
                Cv2 = 1.1 * sqrt(kv) * perfil.raiz_E_fy / perfil.esb_mesa
            else:
                Cv2 = 1.51 * kv * perfil.raiz_E_fy ** 2 / (perfil.esb_mesa ** 2 * perfil.mat.fy)

            Vrdx = perfil.Vplx * Cv2 * phi_v

            return Vrdx if not data else (Vrdx, Vrdx_dados(perfil.Vplx, perfil.Awx, kv, Cv2, elr, elp))

        elif perfil.tipo == 'TUBO CIR':

            return AISC360._Vrd_tubo(perfil, Lv, phi_v, data)

        else:
            return NotImplementedError('Vrdx não implementado para perfis do tipo {}'.format(perfil.tipo))

    @staticmethod
    def _Vrd_tubo(perfil, Lv, phi_v, data):
        """ Determina a força cortante resistênte de cálculo para tubos circulares"""

        if Lv is None:
            raise ValueError('Lv não fornecido')

        Vrd_dados = namedtuple('Vrd_dados', 'Fcr Aw')

        Fcr1 = 1.60 * perfil.mat.E / (sqrt(Lv / perfil.D) * perfil.esb ** (5 / 4))
        Fcr2 = 0.78 * perfil.mat.E / (perfil.esb ** 1.5)

        Fcr = max(Fcr1, Fcr2) if max(Fcr1, Fcr2) < 0.60 * perfil.mat.fy else 0.6 * perfil.mat.fy

        Vrd = Fcr * perfil.Awx * phi_v

        return Vrd if not data else (Vrd, Vrd_dados(Fcr, perfil.Awx))

    @staticmethod
    def _Vn_perfil_IU(perfil, Cv2, a):

        Aw = perfil.Awy
        Afc = perfil.bfs * perfil.tfs if perfil.tipo == 'I SOLDADO' else perfil.bf * perfil.tf
        Aft = perfil.bfs * perfil.tfs if perfil.tipo == 'I SOLDADO' else perfil.bf * perfil.tf

        if 2 * Aw / (Afc + Aft) <= 2.5 and perfil.h / Afc <= 6 and perfil.h / Aft <= 6:
            return Aw * (Cv2 + (1 - Cv2)) / (1.15 * sqrt(1 + (a / perfil.h) ** 2))
        else:
            return Aw * (Cv2 + (1 - Cv2)) / (1.15 * (a / perfil.h + sqrt(1 + (a / perfil.h) ** 2)))

    # cortante em Y
    # -------------
    @staticmethod
    def Vrdy(perfil, a=None, Lv=None, phi_v=0.90, data=False):
        """
        Método que determina a força cortante resistente de cálculo do perfil para cargas aplicadas na direção Y
        de acordo com a **AISC360-16**.

        ver capítulo G da AISC360-16.

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
        phi_v: float, default=0.9
                coeficiente de segurança
        a: float, default=None
            distância entre enrijecedores.
            (só é necessário caso o perfil seja uma instância das classes PerfilI, PerfilILam)
        Lv: float, default=None
            distância entre as seções de forças cortantes máxima e nula.
            (só é necessário caso o perfil seja uma instância da classe TuboCir)
        data: bool, default=False
              Se data=True o método deve retornar os dados utilizados na obtenção de Vrdx.
        Examples
        --------

        Return
        ------
        Vrdy: float
            Força cortante resistênte de cálculo na direção y.
        Vrdy, Vrdx_dados: float, Vrdy_dados
            Força cortante resistênte de cálculo na direção y e dados de cálculo.
            Caso data=True
        """

        if perfil.tipo != 'TUBO CIR':

            Vrdy_dados = namedtuple('Vrdy_dados', 'Vpl kv Cv2 elp elr')

            if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO', 'U SOLDADO'):
                kv = 5.34 if a is None or a / perfil.h > 3 else 5 + 5 / (a / perfil.h) ** 2
            elif perfil.tipo in ('CAIXAO', 'TUBO RET'):
                kv = 5.0
            else:
                kv = 1.2

            elp = 1.1 * sqrt(kv) * perfil.raiz_E_fy
            elr = 1.37 * sqrt(kv) * perfil.raiz_E_fy

            # Cálculo do Cv2
            if perfil.esb_alma <= elp or (perfil.tipo == 'I LAMINADO' and perfil.esb_alma < 2.24 * perfil.raiz_E_fy):
                Cv2 = 1

                if perfil.tipo == 'I LAMINADO' and perfil.esb_alma < 2.24 * perfil.raiz_E_fy:
                    phi_v = 1

                Vrdy = perfil.Vply * Cv2 * phi_v
                return Vrdy if not data else (Vrdy, Vrdy_dados(Vrdy, kv, Cv2, elp, elr))

            elif elp < perfil.esb_alma < elr:
                Cv2 = 1.1 * sqrt(kv) * perfil.raiz_E_fy / perfil.esb_alma
                if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO'):
                    Vrdy = max(perfil.Vply * Cv2, AISC360._Vn_perfil_IU(perfil, Cv2, a)) * phi_v
                    return Vrdy if not data else (Vrdy, Vrdy_dados(Vrdy, kv, Cv2, elp, elr))
                else:
                    Vrdy = perfil.Vply * Cv2 * phi_v
                    return Vrdy if not data else (Vrdy, Vrdy_dados(Vrdy, kv, Cv2, elp, elr))
            else:
                Cv2 = 1.51 * kv * perfil.raiz_E_fy ** 2 / (perfil.esb_alma ** 2 * perfil.mat.fy)
                Vrdy = perfil.Vply * Cv2 * phi_v
                return Vrdy if data else (Vrdy, Vrdy_dados(Vrdy, kv, Cv2, elp, elr))

        elif perfil.tipo == 'TUBO CIR':
            return AISC360._Vrd_tubo(perfil, Lv, phi_v, data)
        else:
            return NotImplementedError('Vrdx não implementado para perfis do tipo {}'.format(perfil.tipo))

    # -> Métodos para determinação da resistência ao momento fletor
    # --------------------------------------------------------------

    # Momento em relação ao eixo X
    # ----------------------------

    # Estado Limite FLT
    @staticmethod
    def _Mnx_LTB(perfil, Lb, Cb, data=False):
        """ Determina a o momento fletor resistente nominal para o estado limite de flambagem lateral com torção para
        perfis I e U de alma compacta (Compact Web), e perfis tubo retangulares, caixão e T."""

        ELU_FLT_dados = namedtuple('ELU_FLT_dados', 'Lb Lp Lr Mn')

        if perfil.tipo in ('I LAMINADO', 'U LAMINADO'):
            return AISC360._Mnx_LTB_IU_CW(perfil, Lb, Cb, data)

        if perfil.tipo == 'I SOLDADO':

            hc = 2 * (perfil.d - perfil.tfs - perfil.hcg)
            elrw = 5.7 * perfil.raiz_E_fy

            if perfil.bissimetrico:
                elpw = 3.76 * perfil.raiz_E_fy
            else:
                hp = 2 * (perfil.d - perfil.tfs - perfil.hpl)
                elpw = (hc / hp) * perfil.raiz_E_fy / (0.54 * perfil.Mpl / perfil.Mrx - 0.09) ** 2
                elpw = min(elpw, elrw)

            if perfil.bissimetrico and perfil.esb_alma <= elpw:
                return AISC360._Mnx_LTB_IU_CW(perfil, Lb, Cb, data)
            elif elpw < perfil.esb_alma <= elrw:
                return AISC360._Mnx_LTB_I_NCW(perfil, Lb, Cb, data)
            else:
                return AISC360._Mnx_LTB_I_SW(perfil, Lb, Cb, data)

        elif perfil.tipo in ('TUBO RET', 'CAIXAO'):

            JA = perfil.J * perfil.A
            Lp = 0.13 * perfil.mat.E * perfil.ry * sqrt(JA) / perfil.Mplx
            Lr = 2 * perfil.mat.E * perfil.ry * sqrt(JA) / (0.7 * perfil.Mrx)
            try:
                Mcrx = 2 * perfil.mat.E * Cb * sqrt(JA) / (perfil.indice_esbeltez(Lb, Lb)[1])
            except ZeroDivisionError:
                Mcrx = perfil.Mplx

        else:
            Lp = 1.76 * perfil.ry * perfil.raiz_E_fy
            E_fy = perfil.raiz_E_fy ** 2

            # OBSERVAÇÃO:
            # ----------
            # para perfis do tipo T laminado e T soldado este cálculo está considerando que o perfil
            # está tendo sua mesa comprimida

            IyJ = perfil.Ix * perfil.J
            Sx_J = perfil.Wx / perfil.J
            Lr = 1.95 * E_fy * sqrt(IyJ) / perfil.Wx * \
                sqrt(2.36 * (1 / E_fy) * perfil.d * Sx_J + 1)

            B = 2.3 * perfil.d * sqrt(perfil.Iy / perfil.J) / Lb
            Mcrx = 1.95 * perfil.mat.E * sqrt(IyJ) * (B + sqrt(1 + B ** 2)) / Lb

        if Lb <= Lp:
            Mn = perfil.Mplx
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))

        elif Lp > Lb >= Lr:
            ctr = 1 if perfil.tipo == 'T LAMINADO' else 0.7  # Coeficiente de tensão residual
            Mnx_flt = Cb * (perfil.Mplx - (perfil.Mplx - ctr * perfil.Mrx) * (Lb - Lp) / (Lr - Lp))
            Mn = min(Mnx_flt, perfil.Mplx)
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))

        else:
            Mn = min(Mcrx, perfil.Mplx)
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))

    @staticmethod
    def _Mnx_LTB_IU_CW(perfil, Lb, Cb, data=False):
        """ Determina a o momento fletor resistente nominal para o estado limite de flambagem lateral com torção para
            perfis I e U de alma compacta (Compact Web)."""

        ELU_FLT_dados = namedtuple('ELU_FLT_dados', 'Lb Lp Lr Mn')

        Lp = 1.76 * perfil.ry * perfil.raiz_E_fy
        E_fy = perfil.raiz_E_fy ** 2

        ho = perfil.d - perfil.tf
        sqrt_IyCw = sqrt(perfil.Iy * perfil.Cw)
        sqrt_Iy_Cw = sqrt(perfil.Iy / perfil.Cw)
        c = 1 if perfil.tipo == 'I LAMINADO' else (ho / 2) * sqrt_Iy_Cw
        rts = sqrt(sqrt_IyCw / perfil.Wx)
        Jc = perfil.J * c
        Sxho = perfil.Wx * ho

        Lr = 1.95 * rts * (E_fy / 0.7) * sqrt(Jc / Sxho + sqrt((Jc / Sxho) ** 2 + 6.76 * (0.7 / E_fy) ** 2))

        if Lb <= Lp:
            Mn = perfil.Mplx
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))
        elif Lp < Lb <= Lr:
            Mnx_flt = Cb * (perfil.Mplx - (perfil.Mplx - 0.7 * perfil.Mrx) * (Lb - Lp) / (Lr - Lp))
            Mn = min(Mnx_flt, perfil.Mplx)
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))
        else:
            Fcr = Cb * pi ** 2 * perfil.mat.E / (Lb / rts) ** 2 * \
                  sqrt(1 + 0.078 * Jc * (Lb / rts) ** 2 / Sxho)
            Mn = Fcr * perfil.Wx
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))

    @staticmethod
    def _Mnx_LTB_I_NCW(perfil, Lb, Cb, data=False):
        """ Determina a o momento fletor resistente nominal para o estado limite de flambagem lateral com torção para
        perfis I bi e monossimétricos de alma compacta (Compact Web) ou não compacta (NonCompact Web) """

        # ---------------------------------------------------------------------------------------
        # IMPORTANTE
        # Para perfis I monossimétricos o método está implementado considerando a mesa superior
        # como a mesa comprimida
        # ---------------------------------------------------------------------------------------

        ELU_FLT_dados = namedtuple('ELU_FLT_dados', 'Lb Lp Lr Mn')

        # Tensão nominal da mesa comprimida Fl
        # ------------------------------------

        Sxc = perfil.Wxs
        Sxi = perfil.Wxi

        Fl = AISC360._Fl(perfil.mat.fy, Sxc, Sxi)

        # Comprimentos limites Lp e Lr
        # ----------------------------

        hc = 2 * (perfil.d - perfil.tfs - perfil.hcg)
        aw = hc * perfil.tw / (perfil.bfs * perfil.tfs)
        rt = perfil.bfs / sqrt(12 * (1 + aw / 6))

        Lp = 1.1 * rt * perfil.raiz_E_fy

        ho = perfil.d - perfil.tfs / 2 - perfil.tfi / 2
        Sxcho = Sxc * ho
        Fl_E = Fl / perfil.mat.E

        Lr = 1.95 * rt * (1 / Fl_E) * sqrt(perfil.J / Sxcho + sqrt((perfil.J / Sxcho) ** 2
                                                                   + 6.76 * Fl_E ** 2))

        # Fator de plastificação Rpc
        # ---------------------------
        Rpc, dados_Rpc = AISC360._Rpc(perfil, True)
        Rpt, dados_Rpt = AISC360._Rpt(perfil, True)

        if Lb >= Lp:
            Mn = min(Rpc * dados_Rpc.Myc, Rpt * dados_Rpt.Myt)
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))
        elif Lp < Lb <= Lr:
            RpcMyc = Rpc * dados_Rpc.Myc
            Mn = min(Cb * (RpcMyc - (RpcMyc - Fl * Sxc) * (Lb - Lp) / (Lr - Lp)), RpcMyc)
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))
        else:
            # Tensão critica Fcr
            Lb_rt2 = (Lb / rt) ** 2
            Fcr = Cb * pi ** 2 * perfil.mat.E / Lb_rt2 * sqrt(1 + 0.078 * perfil.J / Sxcho * Lb_rt2)

            Myc = perfil.mat.fy * perfil.Ws
            Mn = min(Fcr * Sxc, Rpc * Myc)
            Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))

    @staticmethod
    def _Mnx_LTB_I_SW(perfil, Lb, Cb, data=False):
        """ Determina a o momento fletor resistente nominal para o estado limite de flambagem lateral com torção para
        perfis I bi e monossimétricos de alma esbelta (Slender Web) """

        ELU_FLT_dados = namedtuple('ELU_FLT_dados', 'Lb Lp Lr Mn')

        # Comprimentos limites Lp e Lr
        # ----------------------------

        hc = 2 * (perfil.d - perfil.tfs - perfil.hcg)
        aw = min(hc * perfil.tw / (perfil.bfs * perfil.tfs), 10)
        rt = perfil.bfs / sqrt(12 * (1 + aw / 6))

        Lp = 1.1 * rt * perfil.raiz_E_fy
        Lr = pi * rt * perfil.raiz_E_fy * sqrt(1/0.7)

        Rpg = min(AISC360._Rpg(aw, perfil.esb_alma, perfil.mat.E, perfil.mat.fy), 1)

        Sxc = perfil.Wxs
        Sxt = perfil.Wxi

        if Lb <= Lp:
            Mn = min(Rpg * perfil.mat.fy * Sxc, perfil.fy * Sxt)
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))
        elif Lp < Lb <= Lr:
            Fcr = min(Cb * (perfil.mat.fy - (0.3 * perfil.mat.fy) * (Lb - Lp) / (Lr - Lp)), perfil.mat.fy)
            Mn = Fcr * Rpg * Sxc
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))
        else:
            Fcr = (Cb * pi ** 2 * perfil.mat.E / (Lb / rt) ** 2, perfil.mat.fy)
            Mn = Fcr * Rpg * Sxc
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))

    # Estado Limite FLM
    @staticmethod
    def _Mnx_FLB(perfil, data=False):
        """
        Determina o momento fletor resistente nominal em X de uma barra para
        o estado limite último de flambagem local da mesa.
        """

        # Determinando os parametros de esbeltez limites de flambagem da mesa (elp e elr) e
        # momentos nominal critico (Mcr) de acordo com o tipo de perfil

        ELU_FLM_dados = namedtuple('ELU_FLM_dados', 'esb_mesa elpf elrf Mn')

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO', 'T LAMINADO'):
            elpf = 0.38 * perfil.raiz_E_fy
            kc = AISC360._kc(perfil.esb_alma)

            elrf = perfil.raiz_E_fy if not perfil.tipo == 'I SOLDADO' else 0.95 * sqrt(kc / 0.7) * perfil.raiz_E_fy

            if perfil.tipo in ('I LAMINADO', 'T LAMINADO'):
                Mn_ncf = perfil.Mplx - (perfil.Mplx - 0.7 * perfil.Mrx) * (perfil.esb_mesa - elpf) / (elrf - elpf)

                # coeficiente para definir o Mcr
                c_Mcr = 0.7 if perfil.tipo == 'T LAMINADO' else 0.9 * kc

            else:  # Perfis I Soldados

                Rpc, dados = AISC360._Rpc(perfil, dados=True)

                if perfil.bissimetrico and perfil.esb_alma < dados.elpw:
                    Mn_ncf = perfil.Mplx - (perfil.Mplx - 0.7 * perfil.Mrx) * (perfil.esb_mesa - elpf) / (elrf - elpf)
                    c_Mcr = 0.9 * kc
                elif dados.elpw < perfil.esb_alma < dados.elrw:
                    RpcMyc = Rpc * dados.Myc
                    Fl = AISC360._Fl(perfil.mat.fy, perfil.Wxs, perfil.Wxi)
                    Mn_ncf = RpcMyc - (RpcMyc - Fl * perfil.Wxs) * (perfil.esb_mesa_s - elpf) / (elrf - elpf)
                    c_Mcr = 0.9 * kc
                else:
                    hc = 2 * (perfil.d - perfil.tfs - perfil.hcg)
                    aw = min(hc * perfil.tw / (perfil.bfs * perfil.tfs), 10)
                    Rpg = AISC360._Rpg(aw, perfil.esb_alma, perfil.mat.E, perfil.mat.fy)
                    Mn_ncf = (perfil.mat.fy - (0.3 * perfil.mat.fy) * (perfil.esb_mesa_s - elpf) / (elrf - elpf)) * Rpg\
                              * perfil.Wxs
                    c_Mcr = 0.9 * Rpg * kc

            Mcr = c_Mcr * perfil.mat.E * perfil.Wxs / perfil.esb_mesa ** 2

        else:
            # Caso para tubos retangulares e seções caixão
            elpf = 1.12 * perfil.raiz_E_fy
            elrf = 1.4 * perfil.raiz_E_fy if perfil.tipo == 'TUBO RET' else 1.49 * perfil.raiz_E_fy

            ca = 0.38 if perfil.tipo == 'TUBO RET' else 0.34

            bef = min(1.92 * perfil.tf * perfil.raiz_E_fy * (1 - ca * perfil.raiz_E_fy / perfil.esb_mesa), perfil.bint)

            Mn_ncf = min(perfil.Mplx - (perfil.Mplx - perfil.Mrx) * (3.57 * perfil.esb_mesa * perfil.raiz_fy_E - 4),
                         perfil.Mplx)

            Mcr = perfil.mat.fy * AISC360._Sex(perfil, bef)

        # Determinado o momento nominal referente ao estado limite de flambagem local da mesa
        if perfil.esb_mesa <= elpf:
            Mn = perfil.Mplx
            return Mn if not data else (Mn, ELU_FLM_dados(perfil.esb_mesa, elpf, elrf, Mn))
        elif elpf < perfil.esb_mesa <= elrf:
            return Mn_ncf if not data else (Mn_ncf, ELU_FLM_dados(perfil.esb_mesa, elpf, elrf, Mn_ncf))
        else:
            return Mcr if not data else (Mcr, ELU_FLM_dados(perfil.esb_mesa, elpf, elrf, Mcr))

    # Estado Limite FLA
    @staticmethod
    def _Mnx_WLB(perfil, data=False):
        """
        Determina o momento fletor resistente nominal em X de uma barra para
        o estado limite último de flambagem local da alma.
        """

        ELU_FLA_dados = namedtuple('ELU_FLA_dados', 'esb_alma, elpw, elrw, Mn')

        if perfil.tipo in ('TUBO RET', 'CAIXAO'):

            elpw = 2.42 * perfil.raiz_E_fy
            elrw = 5.7 * perfil.raiz_E_fy

            if perfil.esb_alma < elpw:
                Mn = perfil.Mplx
                return Mn if not data else (Mn, ELU_FLA_dados(perfil.esb_alma, elpw, elrw, Mn))

            elif elpw > perfil.esb_alma > elrw:
                Mn = min(perfil.Mplx - (perfil.Mplx - perfil.Mrx) * (0.305 * perfil.esb_alma * perfil.raiz_E_fy - 0.178),
                         perfil.Mplx)
                return Mn if not data else (Mn, ELU_FLA_dados(perfil.esb_alma, elpw, elrw, Mn))

            else:
                hc = 2 * (perfil.h / 2 - perfil.tfs)
                aw = min(hc * perfil.tw / (perfil.b * perfil.tf), 10)
                Rpg = AISC360._Rpg(aw, perfil.esb_alma, perfil.mat.E, perfil.mat.fy)
                Fcr = 0.9 * perfil.E * 4 / perfil.esb_mesa ** 2
                Mn = min(Rpg * perfil.mat.fy * perfil.Wx, Rpg * Fcr * perfil.Wx)
                return Mn if not data else (Mn, ELU_FLA_dados(perfil.esb_alma, elpw, elrw, Mn))

        elif perfil.tipo == 'T LAMINADO':
            elpw = 0.84 * perfil.raiz_E_fy
            elrw = 1.52 * perfil.raiz_E_fy

            if perfil.esb_alma <= elpw:
                Mn = perfil.Mrx
                return Mn if not data else (Mn, ELU_FLA_dados(perfil.esb_alma, elpw, elrw, Mn))

            elif elpw < perfil.esb_alma <= elrw:
                Fcr = (1.43 - 0.515 * perfil.esb_alma * perfil.raiz_fy_E) * perfil.mat.fy
                Mn = Fcr * perfil.Wx
                return Mn if not data else (Mn, ELU_FLA_dados(perfil.esb_alma, elpw, elrw, Mn))

            else:
                Fcr = 1.52 * perfil.mat.E / perfil.esb_alma ** 2
                Mn = Fcr * perfil.Wx
                return Mn if not data else (Mn, ELU_FLA_dados(perfil.esb_alma, elpw, elrw, Mn))

    @staticmethod
    def _Mn_Tubo(perfil, data=False):
        """ Determina o momento fletor resistente nominal para perfis tubo circulares."""

        Mrd_dados = namedtuple("Mrd_dados", "Mn elp elr")

        elp = 0.07 * perfil.mat.E / perfil.mat.fy
        elr = 0.31 * perfil.mat.E / perfil.mat.fy

        if perfil.esb <= elp:
            Mn = perfil.Mplx
            return Mn if not data else (Mn, Mrd_dados(Mn, elp, elr))

        elif elp < perfil.esb <= elr:
            Mn = 0.021 * (perfil.mat.E / perfil.esb + perfil.mat.fy) * perfil.W
            return Mn if not data else (Mn, Mrd_dados(Mn, elp, elr))
        else:
            Fcr = 0.33 * perfil.mat.E / perfil.esb
            Mn = Fcr * perfil.W
            return Mn if not data else (Mn, Mrd_dados(Mn, elp, elr))

    @staticmethod
    def Mrdx(perfil, Lb, Cb, theta_b=0.90, data=False):
        """
        Método responsável por calcular o momento fletor resitente de cálculo em relação ao eixo X
        para uma barra de comprimento destravado  Lb de  acordo com a **AISC360-16**.

        ver capítulo F da AISC360-16.

        Parameters
        ----------
        perfil: objeto PerfilEstrutural
            perfil estrutural.
            podendo ser um objeto de uma das seguintes classes:
                - PerfilI
                - PerfilILam
                - TuboRet
                - TuboCir
                - Caixao
        Lb: float
            comprimento destravado da barra
        Cb: float
            coeficiente Cb determinado conforme item seção F1 do capitulo F da AISC360.
        theta_b: 'float'
            coeficiente de minoração da resistência
        data: bool, default=False
            Se data=True o método deve retornar os dados utilizados na obtenção de Mxrd.

        Examples
        --------

        Returns
        -------
        Mxrd: float
            Momento resistente de cálculo em relação ao eixo X.
        Mxrd, dados: float, objeto Mxrd_dados
            Momento resistente de cálculo em relação ao eixo X e dados de cálculo.
            Caso data=True
        """

        if perfil.tipo == 'U LAMINADO':
            if not data:
                return AISC360._Mnx_LTB(perfil, Lb, Cb) * theta_b
            else:
                Mnx_LTB, dados_LTB = AISC360._Mnx_LTB(perfil, Lb, Cb, data)
                Mrdx = Mnx_LTB * theta_b
                return Mrdx, dados_LTB

        elif perfil.tipo in ('I LAMINADO', 'I SOLDADO'):
            if not data:
                return min(AISC360._Mnx_LTB(perfil, Lb, Cb), AISC360._Mnx_FLB(perfil)) * theta_b
            else:
                Mnx_LTB, dados_LTB = AISC360._Mnx_LTB(perfil, Lb, Cb, data)
                Mnx_FLB, dados_FLB = AISC360._Mnx_FLB(perfil, data)
                Mrdx = min(Mnx_LTB, Mnx_FLB) * theta_b
                return Mrdx, dados_LTB, dados_FLB

        elif perfil.tipo == 'T LAMINADO':
            if not data:
                return min(AISC360._Mnx_LTB(perfil, Lb, Cb), AISC360._Mnx_FLB(perfil), AISC360._Mnx_WLB(perfil)) * \
                   theta_b
            else:
                Mnx_LTB, dados_LTB = AISC360._Mnx_LTB(perfil, Lb, Cb, data)
                Mnx_FLB, dados_FLB = AISC360._Mnx_FLB(perfil, data)
                Mnx_WLB, dados_WLB = AISC360._Mnx_WLB(perfil, data)
                Mrdx = min(Mnx_LTB, Mnx_FLB, Mnx_WLB) * theta_b
                return Mrdx, dados_LTB, dados_FLB, dados_WLB

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Ix >= perfil.Iy:
            if not data:
                return min(AISC360._Mnx_LTB(perfil, Lb, Cb), AISC360._Mnx_FLB(perfil), AISC360._Mnx_WLB(perfil)) * \
                            theta_b
            else:
                Mnx_LTB, dados_LTB = AISC360._Mnx_LTB(perfil, Lb, Cb, data)
                Mnx_FLB, dados_FLB = AISC360._Mnx_FLB(perfil, data)
                Mnx_WLB, dados_WLB = AISC360._Mnx_WLB(perfil, data)
                Mrdx = min(Mnx_LTB, Mnx_FLB, Mnx_WLB) * theta_b
                return Mrdx, dados_LTB, dados_FLB, dados_WLB

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Ix < perfil.Iy:
            if not data:
                return min(AISC360._Mnx_FLB(perfil), AISC360._Mnx_WLB(perfil)) * theta_b
            else:
                Mnx_FLB, dados_FLB = AISC360._Mnx_FLB(perfil, data)
                Mnx_WLB, dados_WLB = AISC360._Mnx_WLB(perfil, data)
                Mrdx = min(Mnx_FLB, Mnx_WLB) * theta_b
                return Mrdx, dados_FLB, dados_WLB

        elif perfil.tipo == 'TUBO CIR':
            if not data:
                return AISC360._Mn_Tubo(perfil) * theta_b
            else:
                Mn, dados = AISC360._Mn_Tubo(perfil, data=True)
                Mrd = Mn * theta_b
                return Mrd, dados

        else:
            NotImplementedError('Método não implementado para perfis do tipo {}'.format(perfil.tipo))

    # Momento em relação ao eixo Y
    # ----------------------------


    @staticmethod
    def _Mny_LTB(perfil, Lb, Cb, data=False):
        """ Determina o momento fletor resistente nominal de uma barra para o estado limite último
            de flambagem lateral com torção em relação ao eixo Y"""

        # Método válido para perfis tubo retangulares e seções caixão

        ELU_FLT_dados = namedtuple('ELU_FLT_dados', 'Lb Lp Lr Mn')

        sqrt_JA = sqrt(perfil.J * perfil.A)

        # Comprimentos limites de plastificação(Lp) e início de escoamento(Lr)
        Lp = 0.13 * perfil.mat.E * perfil.rx * sqrt_JA / perfil.Mply
        Lr = 2 * perfil.mat.E * perfil.rx * sqrt_JA * (0.7 * perfil.Mry)

        if Lb <= Lp:
            Mn = perfil.Mply
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))
        elif Lp < Lb <= Lr:
            Mny_flt = Cb * (perfil.Mply - (perfil.Mply - 0.7 * perfil.Mry) * (Lb - Lp) / (Lr - Lp))
            Mn = min(Mny_flt, perfil.Mply)
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))
        else:
            Mcry = 2 * perfil.E * Cb * sqrt_JA / (perfil.indice_esbeltez(Lb, Lb)[0])
            Mn = min(Mcry, perfil.Mply)
            return Mn if not data else (Mn, ELU_FLT_dados(Lb, Lp, Lr, Mn))

    @staticmethod
    def _Mny_FLB(perfil, data=False):
        """
         Determina o momento fletor resistente nominal em Y de uma barra para
         o estado limite último de flambagem local da mesa.
         """

        ELU_FLM_dados = namedtuple('ELU_FLM_dados', 'esb_mesa elpf elrf Mn')

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO'):

            elpf = 0.38 * perfil.raiz_E_fy
            elrf = perfil.raiz_E_fy if perfil.tipo in ('I LAMINADO', 'U LAMINADO') else 0.95 * perfil.raiz_E_fy * \
                sqrt(AISC360._kc(perfil.esb_alma))

            esb = perfil.esb_mesa

            if esb <= elpf:
                Mn = min(perfil.Mply,  1.6 * perfil.Mry)
                return Mn if not data else (Mn, ELU_FLM_dados(perfil.esb_mesa, elpf, elrf, Mn))

            elif elpf < esb <= elrf:
                Mn = perfil.Mply - (perfil.Mply - 0.7 * perfil.Mry) * (esb - elpf) / (elrf - elpf)
                return Mn if not data else (Mn, ELU_FLM_dados(perfil.esb_mesa, elpf, elrf, Mn))
            else:
                Fcr = 0.69 * perfil.mat.E / perfil.esb_mesa ** 2
                Mn = Fcr * perfil.Wy
                return Mn if not data else (Mn, ELU_FLM_dados(perfil.esb_mesa, elpf, elrf, Mn))

        elif perfil.tipo in ('TUBO RET', 'CAIXAO'):

            elpf = 1.12 * perfil.raiz_E_fy
            elrf = 1.4 * perfil.raiz_E_fy if perfil.tipo == 'TUBO RET' else 1.49 * perfil.raiz_E_fy

            ca = 0.38 if perfil.tipo == 'TUBO RET' else 0.34

            esb = perfil.esb_alma

            bef = min(1.92 * perfil.tw * perfil.raiz_E_fy * (1 - ca * perfil.raiz_E_fy / perfil.esb_alma), perfil.hint)

            if esb <= elpf:
                Mn = perfil.Mply
                return Mn if not data else (Mn, ELU_FLM_dados(perfil.esb_mesa, elpf, elrf, Mn))
            elif elpf < esb <= elrf:
                Mn = min(perfil.Mply - (perfil.Mp - perfil.Mry) * (3.57 * esb * perfil.raiz_E_fy - 4), perfil.Mply)
                return Mn if not data else (Mn, ELU_FLM_dados(perfil.esb_mesa, elpf, elrf, Mn))
            else:
                Mn = perfil.mat.fy * AISC360._Sey(perfil, bef)
                return Mn if not data else (Mn, ELU_FLM_dados(perfil.esb_mesa, elpf, elrf, Mn))

    @staticmethod
    def _Mny_WLB(perfil, data):

        ELU_FLA_dados = namedtuple('ELU_FLA_dados', 'esb_alma, elpw, elrw, Mn')

        elpw = 2.42 * perfil.raiz_E_fy
        elrw = 5.7 * perfil.raiz_E_fy

        if perfil.esb_mesa <= elpw:
            Mn = perfil.Mply
            return Mn if not data else (Mn, ELU_FLA_dados(perfil.esb_mesa, elpw, elrw, Mn))

        elif elpw < perfil.esb_mesa <= elrw:
            Mn_ncw = perfil.Mply - (perfil.Mply - perfil.Mry) * (0.305 * perfil.esb_mesa * perfil.raiz_E_fy - 0.178)
            Mn = min(Mn_ncw, perfil.Mplx)
            return Mn if not data else (Mn, ELU_FLA_dados(perfil.esb_mesa, elpw, elrw, Mn))

        else:
            hc = 2 * (perfil.b / 2 - perfil.tw)
            aw = min(hc * perfil.tf / (perfil.h * perfil.tw), 10)
            Rpg = AISC360._Rpg(aw, perfil.esb_mesa, perfil.mat.E, perfil.mat.fy)
            Fcr = 0.9 * perfil.E * 4 / perfil.esb_alma ** 2
            Mn = min(Rpg * perfil.mat.fy * perfil.Wy, Rpg * Fcr * perfil.Wy)
            return Mn if not data else (Mn, ELU_FLA_dados(perfil.esb_mesa, elpw, elrw, Mn))

    @staticmethod
    def Mrdy(perfil, Lb, Cb, theta_b=0.90, data=False):
        """
        Método responsável por calcular o momento fletor resitente de cálculo em relação ao eixo Y
        para uma barra de comprimento destravado  Lb de  acordo com a **AISC360-16**.

        Ver capítulo F da AISC360-16.

        Parameter
        --------
        perfil: objeto PerfilEstrutural
            perfil estrutural.
            podendo ser um objeto de uma das seguintes classes:
                - PerfilI
                - PerfilILam
                - TuboRet
                - TuboCir
                - Caixao
        Lb: float
            comprimento destravado da barra
        Cb: float
            coeficiente Cb determinado conforme item seção F1 do capitulo F da AISC360.
        theta_b: 'float'
            coeficiente de minoração da resistência
        data: bool, default=False
            Se data=True o método deve retornar os dados utilizados na obtenção de Mrdy.
        Examples
        --------

        Return
        ------
        Mrdy: float
            Momento resistente de cálculo em relação ao eixo y.
        Mrdy, dados: float, objeto Mrdy_dados
            Momento resistente de cálculo em relação ao eixo y e dados de cálculo.
            Caso data=True
        """

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO'):
            if not data:
                return AISC360._Mny_FLB(perfil) * theta_b
            else:
                Mny_FLB, dados_FLB = AISC360._Mny_FLB(perfil, data)
                Mrdy = Mny_FLB * theta_b
                return Mrdy, dados_FLB

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Iy >= perfil.Ix:
            if not data:
                return min(AISC360._Mny_FLB(perfil), AISC360._Mny_LTB(perfil, Cb, Lb), AISC360._Mny_WLB(perfil))*theta_b
            else:
                Mny_LTB, dados_LTB = AISC360._Mny_LTB(perfil, Lb, Cb, data)
                Mny_FLB, dados_FLB = AISC360._Mny_FLB(perfil, data)
                Mny_WLB, dados_WLB = AISC360._Mny_WLB(perfil, data)
                Mrdy = min(Mny_WLB, Mny_FLB, Mny_LTB) * theta_b
                return Mrdy, dados_LTB, dados_FLB, dados_WLB

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Iy < perfil.Ix:
            if not data:
                return min(AISC360._Mny_FLB(perfil), AISC360._Mny_WLB(perfil)) * theta_b
            else:
                Mny_FLB, dados_FLB = AISC360._Mny_FLB(perfil, data)
                Mny_WLB, dados_WLB = AISC360._Mny_WLB(perfil, data)
                Mrdy = min(Mny_WLB, Mny_WLB) * theta_b
                return Mrdy, dados_FLB, dados_WLB
        else:
            NotImplementedError('Método não implementado para perfis do tipo {}'.format(perfil.tipo))

    @staticmethod
    def _Rpc(perfil, dados=False):
        """ Fator de plastificação da alma correspondente a compressão da mesa"""

        rpc_dados = namedtuple('rpc_dados', 'elrw elpw Mp Myc')

        hc = 2 * (perfil.d - perfil.tfs - perfil.hcg)

        elrw = 5.7 * perfil.raiz_E_fy

        if perfil.bissimetrico:
            elpw = 3.76 * perfil.raiz_E_fy
        else:
            hp = 2 * (perfil.d - perfil.tfs - perfil.hpl)
            elpw = (hc / hp) * perfil.raiz_E_fy / (0.54 * perfil.Mpl / perfil.Mrx - 0.09) ** 2
            elpw = min(elpw, elrw)

        Mp = min(perfil.Mplx, 1.6 * perfil.Mrx)
        Myc = perfil.mat.fy * perfil.Wxs

        if perfil.Iys / perfil.Iy > 0.23:

            Rpc = Mp / Myc
            esb = hc / perfil.tw

            if esb > elpw:
                Rpc2 = (Mp / Myc - (Mp / Myc - 1) * (esb - elpw) / (elrw - elpw))
                Rpc = min(Rpc, Rpc2)
            return Rpc if not dados else (Rpc, rpc_dados(elrw, elpw, Mp, Mp, Myc))

        else:
            return 1 if not dados else (1, rpc_dados(elrw, elpw, Mp, Mp, Myc))

    @staticmethod
    def _Rpt(perfil, dados=False):
        """ Fator de plastificação da alma correspondente a tração da mesa"""

        rpt_dados = namedtuple('rpt_dados', 'elrw elpw Mp Myc')

        hc = 2 * (perfil.d - perfil.tfs - perfil.hcg)

        elrw = 5.7 * perfil.raiz_E_fy

        if perfil.bissimetrico:
            elpw = 3.76 * perfil.raiz_E_fy
        else:
            hp = 2 * (perfil.d - perfil.tfs - perfil.hpl)
            elpw = (hc / hp) * perfil.raiz_E_fy / (0.54 * perfil.Mpl / perfil.Mrx - 0.09) ** 2
            elpw = min(elpw, elrw)

        Mp = min(perfil.Mplx, 1.6 * perfil.Mrx)
        Myt = perfil.mat.fy * perfil.Wxi

        if perfil.Iys / perfil.Iy > 0.23:

            Rpc = Mp / Myt
            esb = hc / perfil.tw

            if esb > elpw:
                Rpc2 = (Mp / Myt - (Mp / Myt - 1) * (esb - elpw) / (elrw - elpw))
                Rpc = min(Rpc, Rpc2)
            return Rpc if not dados else (Rpc, rpt_dados(elrw, elpw, Mp, Mp, Myt))

        else:
            return 1 if not dados else (1, rpt_dados(elrw, elpw, Mp, Mp, Myt))

    @staticmethod
    def _Fl(Fy, Sxc, Sxi):

        if Sxi / Sxc >= 0.7:
            return 0.7 * Fy
        else:
            Fl = Fy * Sxi / Sxc
            return Fl if Fl >= 0.5 * Fy else 0.5 * Fy

    @staticmethod
    def _Rpg(aw, esb, E, fy):
        """Fator de redução da resistência a flexão"""
        return 1 - aw / (1200 + 300 * aw) * (esb - 5.7 * sqrt(E / fy))

    @staticmethod
    def _Sex(perfil, bef):
        """ Módulo elástico efetivo, considerando possível flambagem local"""

        # Área(A)
        # -------
        Aefm_sup = bef * perfil.tf
        Am_inf = perfil.bint * perfil.tf
        Aalma = perfil.hint * perfil.tw

        A = Aefm_sup + Am_inf + 2 * Aalma

        # Altura do centro geométrico da seção (ycg)
        # ------------------------------------------

        ycg = (Am_inf * perfil.tf / 2 + 2 * Aalma * (perfil.tf + perfil.h / 2)
               + Aefm_sup * (perfil.h - perfil.tf / 2)) / A

        Imsx = bef * perfil.tf ** 3 / 12
        Imix = perfil.b * perfil.tf ** 3 / 12
        Iax = 2 * perfil.tw * perfil.hint ** 3 / 12

        dmsy = perfil.h - perfil.tf / 2 - ycg
        dmiy = ycg - perfil.tf / 2
        da = abs(perfil.h / 2 - ycg)

        Ix = (Imsx + Aefm_sup * dmsy ** 2) + (Imix + Am_inf * dmiy ** 2) + (Iax + 2 * Aalma * da ** 2)

        Sef = Ix / (perfil.h - ycg)

        return Sef

    @staticmethod
    def _Sey(perfil, bef):
        """ Módulo elástico efetivo, consideranco possível flambagem local"""

        # Área(A)
        # -------
        Aef_ac = bef * perfil.tw  # Area efetiva da alma comprimida
        Amesa = perfil.bint * perfil.tf
        Aalma = perfil.hint * perfil.tw

        A = Aef_ac + 2 * Amesa + Aalma

        # Posição x do centro geométrico da seção (xcg)
        # ------------------------------------------

        xcg = (Aef_ac * perfil.tw / 2 + Aalma * (perfil.b - perfil.tw/2)
               + 2 * Amesa * perfil.b/2) / A

        Iac = bef * perfil.tw ** 3 / 12
        Iat = perfil.h * perfil.tf ** 3 / 12
        Im = 2 * perfil.tf * perfil.bint ** 3 / 12

        dacx = xcg - perfil.tw/2
        datx = perfil.b - perfil.tw / 2 - xcg
        dm = abs(perfil.b / 2 - xcg)

        Iy = (Iac + Aef_ac * dacx ** 2) + (Iat + Aalma * datx ** 2) + (Im + 2 * Amesa * dm ** 2)

        Sef = Iy / xcg

        return Sef

    # -> Método para determinação do momento torsor
    # ----------------------------------------------

    @staticmethod
    def Trd(perfil, L, theta_t=0.9, data=False):

        """
         Método responsável por calcular o momento torsor resitente de cálculo
         para uma barra de comprimento L de  acordo com a **AISC360-**.

         Ver capitulo H da AISC360-16.

         Parameter
         --------
         perfil: objeto PerfilEstrutural
             perfil estrutural.
             podendo ser um objeto de uma das seguintes classes:

                 - TuboRet
                 - TuboCir
                 - Caixao
         L: float
             comprimento da barra (necessário somente para perfis da classe TuboCir)
         gama_a1: 'float'
             coeficiente de minoração da resistência
         data: bool, default=False
             Se data=True o método deve retornar os dados utilizados na obtenção de Trd.
         Examples
         --------

         Return
         ------
         Trd: float
             Momento torsor resistente de cálculo.
         Trd, dados: float, objeto Trd_dados
             Momento torsor resistente de cálculo e dados de cálculo.
             Caso data=True
         """

        if perfil.tipo == 'TUBO CIR':

            Trd_dados = namedtuple('Trd_dados', 'esb Fcr1 Fcr2 Trn')

            Fcr1 = 1.23 * perfil.mat.E / (sqrt(L/perfil.D) * perfil.esb ** (5/4))
            Fcr2 = 0.60 * perfil.mat.E / (perfil.esb ** (3/2))
            Fcr = min(max(Fcr1, Fcr2), 0.6 * perfil.mat.fy)

            Trn = perfil.Wt * Fcr
            Trd = Trn * theta_t

            return Trd if not data else (Trd, Trd_dados(perfil.esb, Fcr1, Fcr2, Trn))

        elif perfil.tipo in ('TUBO RET', 'CAIXAO'):

            Trd_dados = namedtuple('Trd_dados', 'elp elr esb Trn')

            elp = 2.45 * perfil.raiz_E_fy
            elr = 3.07 * perfil.raiz_E_fy

            esb = max(perfil.esb_alma, perfil.esb_mesa)

            if esb <= elp:
                Fcr = 0.6 * perfil.mat.fy
            elif elp < esb <= elr:
                Fcr = 0.6 * perfil.mat.fy * elp / esb
            elif elr < esb <= 260:
                Fcr = 0.458 * pi ** 2 * perfil.mat.E / esb ** 2
            else:
                raise ValueError('A AISC360 não define resistência a torção para perfis com esbeltez > 260')

            Trn = Fcr * perfil.Wt
            Trd = Trn * theta_t

            return Trd if not data else (Trd, Trd_dados(elp, elr, esb, Trn))

        else:
            raise NotImplementedError('Trd não implementado para perfis do tipo {}'.format(perfil.tipo))
