from math import sqrt, pi
from collections import namedtuple


class AISC360:

    # TRAÇÃO
    # --------

    @staticmethod
    def Ntrd_brt(perfil, phi_s=0.9):
        return perfil.Afy * phi_s

    # COMPRESSÃO
    # ----------

    @staticmethod
    def _bef(b, c1, elr, esb, Fy, Fcr):

        if esb <= elr:
            return b
        else:
            c2 = (1 - sqrt(1 - 4 * c1)) / (2 * c1)
            Fel = (c2 * elr / esb) ** 2 * Fy

            bef = b * (1 - c1 * sqrt(Fel / Fcr)) * sqrt(Fel / Fcr)

            return bef if bef < b else b

    @staticmethod
    def _kc(esb):
        kc = 4 / sqrt(esb)
        kc = 0.35 if kc < 0.35 else kc
        kc = 0.76 if kc > 0.76 else kc
        return kc

    @staticmethod
    def _Aef(perfil, Fcr):

        Fy = perfil.mat.fy

        if perfil.tipo in ('I LAMINADO', 'U LAMINADO', 'T LAMINADO'):
            elr_ms = 0.56 * perfil.raiz_E_fy * (Fy / Fcr)

            elr_alm = 0.75 * perfil.raiz_E_fy * (Fy / Fcr) if perfil.tipo == 'T LAMINADO' \
                else 1.49 * perfil.raiz_E_fy * (Fy / Fcr)

            # largura efetiva das mesas
            c1 = 0.22
            b = perfil.bf if perfil.tipo == 'U LAMINADO' else perfil.bf / 2
            bef_ms = AISC360._bef(b, 0.22, elr_ms, perfil.esb_mesa, Fy, Fcr)

            # largura efetiva da alma
            c1 = 0.22 if perfil.tipo == 'T LAMINADO' else 0.18
            h = perfil.d if perfil.tipo == 'T LAMINADO' else perfil.h
            bef_al = AISC360._bef(h, c1, elr_alm, perfil.esb_alma, Fy, Fcr)

            # Número de 'abas'
            n = 4 if perfil.tipo == 'I LAMINADO' else 2

            return perfil.A - n * (b - bef_ms) * perfil.tf - (h - bef_al) * perfil.tw

        if perfil.tipo == 'I SOLDADO':
            kc = AISC360._kc(perfil.esb_alma)
            elr_ms = 0.64 * sqrt(kc) * perfil.raiz_E_fy * (Fy / Fcr)

            elr_alm = 1.49 * perfil.raiz_E_fy * (Fy / Fcr)

            # Largura efetiva das mesas
            bef_mss = AISC360._bef(perfil.bfs / 2, 0.22, elr_ms, perfil.esb_mesa_sup, Fy, Fcr)
            bef_msi = AISC360._bef(perfil.bfi / 2, 0.22, elr_ms, perfil.esb_mesa_inf, Fy, Fcr)

            # Largura efetiva da alma
            bef_alm = AISC360._bef(perfil.dl, 0.18, elr_alm, perfil.esb_alma, Fy, Fcr)

            return perfil.A - (perfil.bfs - 2 * bef_mss) * perfil.tfs \
                   - (perfil.bfi - 2 * bef_msi) * perfil.tfi \
                   - (perfil.dl - bef_alm) * perfil.tw

        if perfil.tipo in ('CAIXAO', 'TUBO RET'):
            elr = 1.40 * perfil.raiz_E_fy * (Fy / Fcr) if perfil.tipo == 'TUBO RET' \
                else 1.49 * perfil.raiz_E_fy * (Fy / Fcr)

            bef_ms = AISC360._bef(perfil.b_int, 0.2, elr, perfil.esb_mesa, Fy, Fcr)
            bef_alm = AISC360._bef(perfil.hint, 0.2, elr, perfil.esb_alma, Fy, Fcr)

            return perfil.A - 2 * (perfil.b_int - bef_ms) * perfil.tf - 2 * (perfil.h_int - bef_alm) * perfil.tw

        if perfil.tipo == 'TUBO CIR':

            elp = 0.11 * perfil.raiz_E_fy ** 2
            elr = 0.45 * perfil.raiz_E_fy ** 2

            if perfil.esb > elp:
                return perfil.A
            elif elp < perfil.esb < elr:
                return (0.038 * perfil.raiz_E_fy ** 2 / perfil.esb + 2 / 3) * perfil.A

    @staticmethod
    def Ncrd(perfil, klx, kly, klz, phi_c=0.9, data=False):

        Ncrd_dados = namedtuple('Ncrd_dados', 'Fe Fy_Fe Fcr Aef')

        Fe = perfil.par_estabilidade(klx, kly, klz).fe

        Fy_Fe = perfil.mat.fy / Fe

        if Fy_Fe <= 2.25:
            Fcr = (0.658 ** Fy_Fe) * perfil.mat.fy
        else:
            Fcr = 0.877 * Fe

        Aef = AISC360._Aef(perfil, Fcr)

        ncrd = Fcr * Aef * phi_c

        return ncrd if not data else (ncrd, Ncrd_dados(Fe, Fy_Fe, Fcr, Aef))

    # CORTANTE
    @staticmethod
    def Vrdx(perfil, phi_v=0.90, data=False):

        Vrdx_dados = namedtuple('Vrdx_dados', 'Vpl kv Cv2 elp elr')

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO', 'T LAMINADO'):
            kv = 1.2
        else:
            kv = 5

        elp = 1.1 * sqrt(kv) * perfil.raiz_E_fy
        elr = 1.37 * sqrt(kv) * perfil.raiz_E_fy

        if perfil.esb_mesa <= elp:
            Cv2 = 1
        elif elp < perfil.esb_mesa <= elr:
            Cv2 = 1.1 * sqrt(kv) * perfil.raiz_E_fy / perfil.esb_mesa
        else:
            Cv2 = 1.51 * kv * perfil.raiz_E_fy ** 2 / (perfil.esb_mesa ** 2 * perfil.mat.fy)

        Vrdx = perfil.Vplx * Cv2 * phi_v

        return Vrdx if not data else (Vrdx, Vrdx_dados(perfil.Vplx, kv, Cv2, elr, elp))

    @staticmethod
    def _Vn_perfil_IU(perfil, Cv2, a):

        Aw = perfil.Awy
        Afc = perfil.bfs * perfil.tfs if perfil.tipo == 'I SOLDADO' else perfil.bf * perfil.tf
        Aft = perfil.bfs * perfil.tfs if perfil.tipo == 'I SOLDADO' else perfil.bf * perfil.tf

        if 2 * Aw / (Afc + Aft) <= 2.5 and perfil.h / Afc <= 6 and perfil.h / Aft <= 6:
            return Aw * (Cv2 + (1 - Cv2)) / (1.15 * sqrt(1 + (a / perfil.h) ** 2))
        else:
            return Aw * (Cv2 + (1 - Cv2)) / (1.15 * (a / perfil.h + sqrt(1 + (a / perfil.h) ** 2)))

    @staticmethod
    def Vrdy(perfil, a=None, phi_v=0.90):

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
            return perfil.Vply * Cv2 * phi_v

        elif elp < perfil.esb_alma < elr:
            Cv2 = 1.1 * sqrt(kv) * perfil.raiz_E_fy / perfil.esb_alma
            if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO'):
                return max(perfil.Vply * Cv2, AISC360._Vn_perfil_IU(perfil, Cv2, a)) * phi_v
            else:
                return perfil.Vply * Cv2 * phi_v
        else:
            Cv2 = 1.51 * kv * perfil.raiz_E_fy ** 2 / (perfil.esb_alma ** 2 * perfil.mat.fy)
            return perfil.Vply * Cv2 * phi_v

    #  MOMENTO FLETOR

    @staticmethod
    def _Mnx_LTB_CW(perfil, Lb, Cb):
        """ Determina a o momento nominal para o estado limite de flambagem lateral com torção para
        perfis I e U de alma compacta, e perfis tubo retangulares, caixão e T """

        if perfil.tipo in ('TUBO RET', 'CAIXAO'):

            sqrt_JA = sqrt(perfil.J * perfil.A)
            Lp = 0.13 * perfil.mat.E * perfil.ry * sqrt_JA / perfil.Mplx
            Lr = 2 * perfil.mat.E * perfil.ry * sqrt_JA * (0.7 * perfil.Mrx)
            Mcrx = 2 * perfil.E * Cb * sqrt_JA / (perfil.indice_esbeltez(Lb, Lb)[0])

        else:
            Lp = 1.76 * perfil.ry * perfil.raiz_E_fy
            E_fy = perfil.raiz_E_fy ** 2

            if perfil.tipo in ('I LAMINADO', 'U LAMINADO'):

                ho = perfil.d - perfil.tf
                sqrt_Iy_Cw = sqrt(perfil.Iy / perfil.Cw)
                c = 1 if perfil.tipo == 'I LAMINADO' else (ho / 2) * sqrt_Iy_Cw
                rts = sqrt(sqrt_Iy_Cw / perfil.Wx)
                Jc = perfil.J * c
                Sxho = perfil.Wx * ho

                Lr = 1.95 * rts * (E_fy / 0.7) * \
                     sqrt(Jc / Sxho + sqrt(Jc / Sxho ** 2 + 6.76 * (0.7 / E_fy) ** 2))

                Fcr = Cb * pi ** 2 * perfil.mat.E / (Lb / rts) ** 2 * \
                      sqrt(1 + 0.078 * Jc * (Lb / rts) ** 2 / Sxho)

                Mcrx = Fcr * perfil.Wx

            else:
                IyJ = perfil.Ix * perfil.J
                Sx_J = perfil.Wx / perfil.J
                Lr = 1.95 * E_fy * sqrt(IyJ) / perfil.Wx * \
                     sqrt(2.36 * (1 / E_fy) * perfil.d * Sx_J + 1)

                B = 2.3 * (perfil.d / Lb) * sqrt(perfil.Iy / perfil.J)
                Mcrx = 1.95 * perfil.mat.E * sqrt(IyJ) * (B + sqrt(1 + B ** 2)) / Lb

        if Lb <= Lp:
            return perfil.Mplx

        elif Lp > Lb >= Lr:
            ctr = 1 if perfil.tipo == 'T LAMINADO' else 0.7  # Coeficiente de tensão residual
            Mnx_flt = Cb * (perfil.Mplx - (perfil.Mplx - ctr * perfil.Mrx) * (Lb - Lp) / (Lr - Lp))
            return Mnx_flt if Mnx_flt < perfil.Mplx else perfil.Mplx

        else:
            return Mcrx if Mcrx < perfil.Mplx else perfil.Mplx

    @staticmethod
    def _Mnx_LTB_NCW(perfil, Cb, Lb):

        # Tensão nominal da mesa comprimida Fl
        # ------------------------------------

        Sxc = perfil.Wxs
        Sxi = perfil.Wxi

        Fl = AISC360._Fl(perfil.mat.fy, Sxc, Sxi)

        # Comprimentos limites Lp e Lr
        # ----------------------------

        hc = 2 * (perfil.d - perfil.tfs - perfil.hcg)
        aw = hc * perfil.tw / (perfil.bfs * perfil.tfc)
        rt = perfil.bfs / sqrt(12 * (1 + aw / 6))

        Lp = 1.1 * rt * perfil.raiz_E_fy

        ho = perfil.d - perfil.tfs / 2 - perfil.tfi / 2
        Sxcho = Sxc * ho
        Fl_E = Fl / perfil.mat.E

        Lr = 1.95 * rt * (1 / Fl_E) * sqrt(perfil.J / Sxcho + sqrt((perfil.J / Sxcho) ** 2
                                                                   + 6.76 * Fl_E ** 2))

        # Fator de plastificação Rpc
        # ---------------------------
        Rpc = AISC360._Rpc(perfil)

        # Tensão critica Fcr
        Lb_rt2 = (Lb / rt) ** 2
        Fcr = Cb * pi ** 2 * perfil.mat.E / Lb_rt2 * sqrt(1 + 0.078 * perfil.J / Sxcho * Lb_rt2)

        if Lb >= Lp:
            return perfil.Mplx
        elif Lp < Lb <= Lr:
            Myc = perfil.mat.fy * perfil.Ws
            RpcMyc = Rpc * Myc
            Mn = Cb * (RpcMyc - (RpcMyc - Fl * Sxc) * (Lb - Lp) / (Lr - Lp))
            return Mn if Mn <= RpcMyc else RpcMyc
        else:
            Mn = Fcr * Sxc
            Myc = perfil.mat.fy * perfil.Ws
            RpcMyc = Rpc * Myc
            return Mn if Mn <= RpcMyc else RpcMyc

    @staticmethod
    def _Mnx_FLB(perfil):

        # Determinando os parametros de esbeltez limites de flambagem da mesa (elp e elr) e
        # momentos nominal critico (Mcr) de acordo com o tipo de perfil

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO', 'T LAMINADO'):
            elpf = 0.38 * perfil.raiz_E_fy
            kc = AISC360._kc(perfil.esb_alma)

            elrf = perfil.raiz_E_fy if not perfil.tipo == 'I SOLDADO' else 0.95 * sqrt(kc / 0.7) * perfil.raiz_E_fy

            if perfil.tipo == 'T LAMINADO':
                Mcr = 0.7 * perfil.mat.E * perfil.Wxs / perfil.esb_mesa ** 2
            else:
                Mcr = 0.9 * perfil.mat.E * kc * perfil.Wx / perfil.esb_mesa ** 2

        elif perfil.tipo in ('TUBO RET', 'CAIXAO'):
            elpf = 1.12 * perfil.raiz_E_fy
            elrf = 1.4 * perfil.raiz_E_fy if perfil.tipo == 'TUBO RET' else 1.49 * perfil.raiz_E_fy

            c1 = 0.38 if perfil.tipo == 'TUBO RET' else 0.34
            Fy = perfil.mat.fy
            bef = AISC360._bef(perfil.bint, c1, elrf, perfil.esb_mesa, Fy, Fy)

            Mcr = perfil.mat.Fy * AISC360._Sex(perfil, bef)

        # Determinado o momento nominal referente ao estado limite de flambagem local da mesa
        if perfil.esb_mesa >= elpf:
            return perfil.Mplx

        elif elpf > perfil.esb_mesa >= elrf:

            if perfil.tipo in ('I SOLDADO', 'T LAMINADO'):
                return perfil.Mplx - (perfil.Mplx - 0.7 * perfil.Mrx) * (perfil.esb_mesa - elpf) / (elpf - elrf)

            if perfil.tipo == 'I LAMINADO':

                Rpc, dados = AISC360._Rpc(perfil, dados=True)

                if perfil.bissimetrico and perfil.esb_alma < dados.elpw:
                    return perfil.Mplx - (perfil.Mplx - 0.7 * perfil.Mrx) * (perfil.esb_mesa - elpf) / (elpf - elrf)
                else:
                    RpcMyc = Rpc * dados.Myc
                    Fl = AISC360._Fl(perfil.mat.fy, perfil.Wxs, perfil.Wxi)
                    return RpcMyc - (RpcMyc - Fl * perfil.Wxs) * (perfil.esb_mesa - elpf) / (elpf - elrf)
            else:
                Mn = perfil.Mplx - (perfil.Mplx - perfil.Mrx) * (3.57 * perfil.esb_mesa * perfil.raiz_fy_E - 4)
                return Mn if Mn <= perfil.Mplx else Mn
        else:
            return Mcr

    @staticmethod
    def _Mnx_WLB(perfil):

        if perfil.tipo in ('TUBO RET', 'CAIXAO'):

            elpw = 2.42 * perfil.raiz_E_fy
            elrw = 1.4 * perfil.raiz_E_fy

            if perfil.esb_alma < elpw:
                return perfil.Mplx

            elif elpw > perfil.esb_alma > elrw:
                Mn = perfil.Mplx - (perfil.Mplx - perfil.Mrx) * (0.305 * perfil.esb_alma * perfil.raiz_E_fy - 0.178)
                return min(Mn, perfil.Mplx)

            else:
                hc = 2 * (perfil.h / 2 - perfil.tfs)
                aw = min(hc * perfil.tw / (perfil.b * perfil.tf), 10)
                Rpg = AISC360._Rpg(aw, perfil.esb_alma, perfil.mat.E, perfil.mat.fy)
                Fcr = 0.9 * perfil.E * 4 / perfil.esb_mesa ** 2
                return min(Rpg * perfil.mat.fy * perfil.Wx, Rpg * Fcr * perfil.Wx)

        if perfil.tipo == 'T LAMINADO':
            elpw = 1.52 * perfil.raiz_E_fy
            elrw = 0.84 * perfil.raiz_E_fy

            if perfil.esb_alma <= elpw:
                return perfil.Mrx

            elif elpw < perfil.esb_alma <= elrw:
                Fcr = (1.43 - 0.515 * perfil.esb_alma * perfil.raiz_fy_E) * perfil.mat.fy
                return Fcr * perfil.Wx

            else:
                Fcr = 1.52 * perfil.mat.E / perfil.esb_alma ** 2
                return Fcr * perfil.Wx

    @staticmethod
    def _Mn_Tubo(perfil):

        elp = 0.07 * perfil.mat.E / perfil.mat.fy
        elr = 0.31 * perfil.mat.E / perfil.mat.fy

        if perfil.esb <= elp:
            return perfil.Mplx
        elif elp < perfil.esb <= elr:
            return 0.021 * (perfil.mat.E / perfil.esb + perfil.mat.fy) * perfil.W
        else:
            Fcr = 0.33 * perfil.mat.E / perfil.esb
            return Fcr * perfil.W

    @staticmethod
    def Mrdx(perfil, Lb, Cb, theta_b=0.90):

        if perfil.tipo == 'U LAMINADO':
            return AISC360._Mnx_LTB_CW(perfil, Lb, Cb) * theta_b

        elif perfil.tipo == 'I LAMINADO':
            return min(AISC360._Mnx_LTB_CW(perfil, Lb, Cb), AISC360._Mnx_FLB(perfil)) * theta_b

        elif perfil.tipo == 'I SOLDADO':
            return min(AISC360._Mnx_LTB_NCW(perfil, Cb, Lb), AISC360._Mnx_FLB(perfil)) * theta_b

        elif perfil.tipo == 'T LAMINADO':
            return min(AISC360._Mnx_LTB_CW(perfil, Lb, Cb), AISC360._Mnx_FLB(perfil), AISC360._Mnx_WLB(perfil)) * \
                   theta_b

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Ix > perfil.Iy:
            return min(AISC360._Mnx_LTB_CW(perfil, Lb, Cb), AISC360._Mnx_FLB(perfil), AISC360._Mnx_WLB(perfil)) * \
                   theta_b

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Ix < perfil.Iy:
            return min(AISC360._Mnx_FLB(perfil), AISC360._Mnx_WLB(perfil)) * theta_b

        elif perfil.tipo == 'TUBO CIR':
            return AISC360._Mn_Tubo(perfil) * theta_b

        else:
            print('Método não implementado para perfis do tipo {}'.format(perfil.tipo))

    @staticmethod
    def _Mny_LTB(perfil, Cb, Lb):

        sqrt_JA = sqrt(perfil.J * perfil.A)

        # Comprimentos limites de plastificação(Lp) e início de escoamento(Lr)
        Lp = 0.13 * perfil.mat.E * perfil.rx * sqrt_JA / perfil.Mply
        Lr = 2 * perfil.mat.E * perfil.rx * sqrt_JA * (0.7 * perfil.Mry)

        Mcrx = 2 * perfil.E * Cb * sqrt_JA / (perfil.indice_esbeltez(Lb, Lb)[1])

        if Lb <= Lp:
            return perfil.Mply
        elif Lp > Lb >= Lr:
            Mny_flt = Cb * (perfil.Mply - (perfil.Mply - 0.7 * perfil.Mry) * (Lb - Lp) / (Lr - Lp))
            return Mny_flt if Mny_flt < perfil.Mply else perfil.Mply

    @staticmethod
    def _Mny_FLB(perfil):

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO'):

            elpf = 0.38 * perfil.raiz_E_fy
            elrf = perfil.raiz_E_fy if perfil.tipo in ('I LAMINADO', 'U LAMINADO') else 0.95 * perfil.raiz_E_fy * \
                sqrt(AISC360._kc(perfil.esb_alma))

            esb = perfil.esb_mesa

            if esb <= elpf:
                return min(perfil.Mply,  1.6 * perfil.Mry)
            elif elpf < esb <= elrf:
                return perfil.Mply - (perfil.Mply - 0.7 * perfil.Mry) * (esb - elpf) / (elrf - elpf)
            else:
                Fcr = 0.69 * perfil.mat.E / perfil.esb_mesa ** 2
                return Fcr * perfil.Wy

        elif perfil.tipo in ('TUBO RET', 'CAIXAO'):

            elpf = 1.12 * perfil.raiz_E_fy
            elrf = 1.4 * perfil.raiz_E_fy if perfil.tipo == 'TUBO RET' else 1.49 * perfil.raiz_E_fy

            c1 = 0.38 if perfil.tipo == 'TUBO RET' else 0.34
            Fy = perfil.mat.fy
            bef = AISC360._bef(perfil.bint, c1, elrf, perfil.esb_alma, Fy, Fy)

            esb = perfil.esb_alma

            if esb <= elpf:
                return perfil.Mply
            elif elpf < esb <= elrf:
                return min(perfil.Mply - (perfil.Mp - perfil.Mry) * (3.57 * esb * perfil.raiz_E_fy - 4), perfil.Mply)
            else:
                return perfil.mat.Fy * AISC360._Sey(perfil, bef)

    @staticmethod
    def _Mny_WLB(perfil):

        elpw = 2.42 * perfil.raiz_E_fy
        elrw = 1.4 * perfil.raiz_E_fy

        if perfil.esb_mesa < elpw:
            return perfil.Mply

        elif elpw > perfil.esb_mesa > elrw:
            Mn = perfil.Mply - (perfil.Mply - perfil.Mry) * (0.305 * perfil.esb_mesa * perfil.raiz_E_fy - 0.178)
            return min(Mn, perfil.Mplx)

        else:
            hc = 2 * (perfil.b / 2 - perfil.tw)
            aw = min(hc * perfil.tf / (perfil.h * perfil.tw), 10)
            Rpg = AISC360._Rpg(aw, perfil.esb_mesa, perfil.mat.E, perfil.mat.fy)
            Fcr = 0.9 * perfil.E * 4 / perfil.esb_alma ** 2
            return min(Rpg * perfil.mat.fy * perfil.Wy, Rpg * Fcr * perfil.Wy)

    @staticmethod
    def Mrdy(perfil, Lb, Cb, theta_b=0.90):

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO'):
            return AISC360._Mny_FLB(perfil) * theta_b

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Iy > perfil.Ix:
            return min( AISC360._Mny_FLB(perfil), AISC360._Mny_LTB(perfil, Cb, Lb), AISC360._Mny_WLB(perfil)) * theta_b

        elif perfil.tipo in ('TUBO RET', 'CAIXAO') and perfil.Iy < perfil.Ix:
            return min(AISC360._Mny_FLB(perfil), AISC360._Mny_WLB(perfil)) * theta_b

        else:
            print('Método não implementado para perfis do tipo {}'.format(perfil.tipo))

    @staticmethod
    def _Rpc(perfil, dados=False):

        rpc_dados = namedtuple('rpc_dados', 'elrw elpw Mp Myc')

        hc = 2 * (perfil.d - perfil.tfs - perfil.hcg)

        elrw = 5.7 * perfil.raiz_E_fy

        if perfil.bissimetrico:
            elpw = 3.76 * perfil.raiz_E_fy
        else:
            hp = 2 * (perfil.d - perfil.tfs - perfil.hpl)
            elpw = (hc / hp) * perfil.raiz_E_fy / (0.54 * perfil.Mpl / perfil.Mrx - 0.09) ** 2
            elpw = elpw if elpw >= elrw else elrw

        Mp = min(perfil.Mplx, 1.6 * perfil.Mrx)
        Myc = perfil.mat.fy * perfil.Ws

        if perfil.Iys / perfil.Iy > 0.23:

            Rpc = Mp / Myc
            esb = hc / perfil.tw

            if esb > elpw:
                Rpc2 = (Mp / Myc - (Mp / Myc - 1) * (esb - elpw) / (elrw - elpw))
                Rpc = Rpc2 if Rpc2 <= Rpc else Rpc
            return Rpc if not dados else (Rpc, rpc_dados(elrw, elpw, Mp, Mp, Myc))

        else:
            return 1 if not dados else (1, rpc_dados(elrw, elpw, Mp, Mp, Myc))

    @staticmethod
    def _Fl(Fy, Sxc, Sxi):

        if Sxi / Sxc >= 0.7:
            return 0.7 * Fy
        else:
            Fl = Fy * Sxi / Sxc
            return Fl if Fl >= 0.5 * Fy else 0.5 * Fy

    @staticmethod
    def _Rpg(aw, esb, E, fy):
        return 1 - aw / (1200 + 300 * aw) * (esb - 5.7 * sqrt(E / fy))

    @staticmethod
    def _Sex(perfil, bef):
        """ Módulo elástico efetivo, consideranco possível flambagem local"""

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

        Ix = Imsx + Aefm_sup * dmsy ** 2 + \
             Imix + Am_inf * dmiy ** 2 + \
             Iax + 2 * Aalma * da ** 2

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

        Iy = Iac + Aef_ac * dacx ** 2 + \
             Iat + Aalma * datx ** 2 + \
             Im + 2 * Amesa * dm ** 2

        Sef = Iy / (xcg)

        return Sef
