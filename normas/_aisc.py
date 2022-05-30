from perfis import *


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
    def _Aef(perfil, Fcr):

        Fy = perfil.material.fy

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

            kc = 4 / sqrt(perfil.esb_alma)
            kc = 0.35 if kc < 0.35 else kc
            kc = 0.76 if kc > 0.76 else kc
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
                return (0.038 * perfil.raiz_E_fy ** 2 / perfil.esb + 2/3) * perfil.A

    @staticmethod
    def Ncrd(perfil, klx, kly, klz, phi_c=0.9, data=False):

        Ncrd_dados = namedtuple('Ncrd_dados', 'Fe Fy_Fe Fcr Aef')

        Fe = perfil.par_estabilidade(klx, kly, klz).fe

        Fy_Fe = perfil.material.fy / Fe

        if Fy_Fe <= 2.25:
            Fcr = (0.658 ** Fy_Fe) * perfil.material.fy
        else:
            Fcr = 0.877 * Fe

        Aef = AISC360._Aef(perfil, Fcr)

        ncrd = Fcr * Aef * phi_c

        return ncrd if not data else (ncrd, Ncrd_dados(Fe, Fy_Fe, Fcr, Aef))

    # CORTANTE
    @staticmethod
    def Vrdx(perfil, phi_v=0.90, data = False):

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
            Cv2 = 1.51 * kv * perfil.raiz_E_fy ** 2 / (perfil.esb_mesa ** 2 * perfil.material.fy)

        Vrdx = perfil.Vplx * Cv2 * phi_v

        return Vrdx if not data else (Vrdx, Vrdx_dados(perfil.Vplx, kv, Cv2, elr, elp))

    @staticmethod
    def _Vn_perfil_IU(perfil, Cv2, a):

        Aw = perfil.Awy
        Afc = perfil.bfs * perfil.tfs if perfil.tipo == 'I SOLDADO' else perfil.bf * perfil.tf
        Aft = perfil.bfs * perfil.tfs if perfil.tipo == 'I SOLDADO' else perfil.bf * perfil.tf

        if 2 * Aw / (Afc + Aft) <= 2.5 and perfil.h/Afc <= 6 and perfil.h/Aft <= 6:
            return Aw * (Cv2 + (1 - Cv2)) / (1.15 * sqrt(1 + (a/perfil.h) ** 2))
        else:
            return Aw * (Cv2 + (1 - Cv2)) / (1.15 * (a/perfil.h + sqrt(1 + (a/perfil.h) ** 2)))

    @staticmethod
    def Vrdy(perfil, a=None,  phi_v=0.90):

        Vrdy_dados = namedtuple('Vrdy_dados', 'Vpl kv Cv2 elp elr')

        if perfil.tipo in ('I LAMINADO', 'I SOLDADO', 'U LAMINADO', 'U SOLDADO'):
            kv = 5.34 if a is None or a/perfil.h > 3 else 5 + 5 / (a / perfil.h) ** 2
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
            Cv2 = 1.51 * kv * perfil.raiz_E_fy ** 2 / (perfil.esb_alma ** 2 * perfil.material.fy)
            return perfil.Vply * Cv2 * phi_v


#MOMENTO FLETOR