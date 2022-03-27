from math import sqrt, pi
from secao import SecaoGenerica

class SecaoMetalica(SecaoGenerica):

    def __init__(self, a, ix, iy, j, wx, wy, zx, zy, xo, yo, cw, material, simetria):

        super().__init__(a, ix, iy, j, material, wx, wy, xo, yo, cw, simetria):

        self.zx = zx
        self.zy = zy

    # -------------------------------------------------------------------------------------
    # --------------------------Verificações de resistência--------------------------------
    # -------------------------------------------------------------------------------------

    # ----------------------------------NBR 8800-------------------------------------------

    # TRAÇÂO
    # --------

    def resistEscSecBruta_NBR8800(self, gama_a1=1.1):
        return self.a * self.material.fy / gama_a1

    # COMPRESSÃO
    # -----------
    def calcularIndEsbRed(self, klx, kly, klz, Q=1):
        ne = self.calcularNe(klx, kly, klz)
        ind_esb_red = sqrt(Q * self.a * self.material.fy) / ne
        return ind_esb_red

    def calcularFatorRedComp(self, indice_esb_reduzido):
        if indice_esb_reduzido <= 1.5:
            fator_red_comp = 0.658 ** (indice_esb_reduzido ** 2)
        else:
            fator_red_comp = 0.877 / (indice_esb_reduzido) ** 2
        return fator_red_comp

    def calcularFatorRedQs(self):
        pass

    def calcularFatorRedQa(self, fator_red_com):
        pass

    def calcularFatorRedQ(self, fator_red_com):
        qa = self.calcularFatorRedQa(fator_red_com)
        qs = self.calcularFatorRedQs()
        return qa * qs

    def ncrd_NBR8800(self, klx, kly, klz, gama_a1=1.1):

        indice_esb_red = self.calcularIndEsbRed(klx, kly, klz)
        fator_red_comp = self.calcularFatorRedComp(indice_esb_red)
        q = self.calcularFatorRedQ(fator_red_comp)

        indice_esb_red = self.calcularIndEsbRed(klx, kly, klz, q)
        fator_red_comp = self.calcularFatorRedComp(indice_esb_red)
        ncrd = fator_red_comp * q * self.a * self.material.fy / gama_a1
        return ncrd

    # CORTANTE
    # -----------
    @property
    def awx(self):
        pass

    @property
    def awy(self):
        pass

    def calcularParEsbLimSecCom(self, kv):
        return 1.1 * sqrt(kv * self.material.E / self.material.fy)

    def calcularParEsbLimSemiCom(self, kv):
        return 1.37 * sqrt(kv * self.material.E / self.material.fy)

    def calcularVpl(self, aw):
        vpl = 0.60 * aw * self.material.fy
        return vpl

    # EM X
    # -------------
    def calcularParEsbVrdx(self):
        pass

    def definirKvVrdx(self):
        pass

    def vrdx_NBR8800(self, a, gama_a1=1.1):

        par_esbeltez = self.calcularParEsbVrdx()
        kv = self.definirKvVrdx(a)
        par_esbeltez_lim_p = self.calcularParEsbLimSecCom(kv)
        par_esbeltez_lim_r = self.calcularParEsbLimSecSemiCom(kv)
        vplx = self.calcularVpl(self.awx)

        if par_esbeltez <= par_esbeltez_lim_p:
            vrdx = vplx / gama_a1
        elif par_esbeltez_lim_p < par_esbeltez <= par_esbeltez_lim_r:
            vrdx = (par_esbeltez_lim_p / par_esbeltez) * (vplx / gama_a1)
        else:
            vrdx = 1.24 * (par_esbeltez_lim_p / par_esbeltez) ** 2 * (vplx / gama_a1)

        return vrdx

    # CORTANTE EM Y
    # -------------
    def calcularParEsbVrdy(self):
        pass

    def definirKvVrdy(self, a):
        pass

    def vrdy_NBR8800(self, a, gama_a1=1.1):

        par_esbeltez = self.calcularParEsbVrdy()
        kv = self.definirKvVrdy(a)
        par_esbeltez_lim_p = self.calcularParEsbLimSecCom(kv)
        par_esbeltez_lim_r = self.calcularParEsbLimSecSemiCom(kv)
        vply = self.calcularVpl(self.awy)

        if par_esbeltez <= par_esbeltez_lim_p:
            vrdy = vply / gama_a1
        elif par_esbeltez_lim_p < par_esbeltez <= par_esbeltez_lim_r:
            vrdy = (par_esbeltez_lim_p / par_esbeltez) * (vply / gama_a1)
        else:
            vrdy = 1.24 * (par_esbeltez_lim_p / par_esbeltez) ** 2 * (vply / gama_a1)

        return vrdy

    # MOMENTO EM X
    # ------------

    @property
    def mplx(self):
        return self.zx * self.material.fy

    # Estado Limite FLT
    def calcularParEsbMrdxFLT(self, lb):
        pass

    def calcularParEsbLimPlasMrdxFLT(self):
        pass

    def calcularParEsbLimEscMrdxFLT(self):
        pass

    def mrxFLT(self):
        pass

    def mcrxFLT(self, cb, lb):
        pass

    def mrdxFLT(self, cb, lb):
        pass

    # Estado Limite FLM
    def calcularParEsbMrdxFLM(self):
        pass

    def calcularParEsbLimPlasMrdxFLM(self):
        pass

    def calcularParEsbLimEscMrdxFLM(self):
        pass

    def mrxFLM(self):
        pass

    def mcrxFLM(self):
        pass

    def mrdxFLM(self):
        pass

    # Estado Limite FLA
    def calcularParEsbMrdxFLA(self):
        pass

    def calcularParEsbLimPlasMrdxFLA(self):
        pass

    def calcularParEsbLimEscMrdxFLA(self):
        pass

    def mrxFLA(self):
        pass

    def mcrxFLA(self):
        pass

    def mrdxFLA(self):
        pass

    def mrdx_NBR8800(self, cb, lb):
        pass

    # MOMENTO EM Y
    # ------------

    @property
    def mply(self):
        return self.material.zy

    # Estado Limite FLT
    def calcularParEsbMrdyFLT(self):
        pass

    def calcularParEsbLimPlasMrdyFLT(self):
        pass

    def calcularParEsbLimEscMrdyFLT(self):
        pass

    def mryFLT(self):
        pass

    def mcryFLT(self, Cb):
        pass

    def mrdyFLT(self, Cb):
        pass

    # Estado Limite FLM
    def calcularParEsbMrdyFLM(self):
        pass

    def calcularParEsbLimPlasMrdyFLM(self):
        pass

    def calcularParEsbLimEscMrdyFLM(self):
        pass

    def mryFLM(self):
        pass

    def mcryFLM(self, par_esb):
        pass

    def mrdyFLM(self):
        pass

    # Estado Limite FLA
    def calcularParEsbMrdyFLA(self):
        pass

    def calcularParEsbLimPlasMrdyFLA(self):
        pass

    def calcularParEsbLimEscMrdyFLA(self):
        pass

    def mryFLA(self):
        pass

    def mcryFLA(self):
        pass

    def mrdyFLA(self):
        pass

    def mrdy_NBR8800(self, cb, lb):
        pass

class SecaoILaminada(SecaoMetalica):

    def __init__(self, nome, material):

        self.ht = ht
        self.bf = bf
        self.tf = tf
        self.tw = tw
        self.r = r

        simetria = [True, True]
        super().__init__(a, ix, iy, j, wx, wy, zx, zy, xo, yo, cw, material, simetria)


    @property
    def hw(self):
        return self.ht - 2 * self.tf
    @property
    def h(self):
        return self.ht - 2 * self.tf - 2 * self.r

    #-------------------------------------------------------------------------------------
    #--------------------------Verificações de resistência--------------------------------
    #-------------------------------------------------------------------------------------


    #----------------------------------NBR 8800-------------------------------------------
    #COMPRESSÃO
    #-----------
    def calcularFatorRedQs(self):
        b = self.bf/2
        t = self.tf
        qs = 1

        ind_esbeltez = b/t
        ind_esbeltez_lim_p = 0.56 * sqrt(self.material.e / self.material.fy)
        ind_esbeltez_lim_r = 1.03 * sqrt(self.material.e / self.material.fy)

        if ind_esbeltez_lim_p < ind_esbeltez <= ind_esbeltez_lim_r:
            qs = 1.415 - 0.74 * (b / t) * sqrt( self.material.fy / self.material.e)
        elif ind_esbeltez > ind_esbeltez_lim_r:
            qs = 0.69 * self.material.e / (self.material.fy * (b / t)**2)

        return qs

    def calcularFatorRedQa(self, fator_red_com):

        t = self.tw
        b = self.h
        tensao = self.fy * fator_red_com
        ca = 0.34

        bef = 1.92 * t * sqrt(self.material.e / tensao) * (1 - ca / (b/t) * sqrt(self.material.e / tensao))
        bef = bef if bef < b else b

        aef = self.a - (b - bef) * t

        return aef/self.a

    #CORTANTE
    #-----------
    @property
    def awx(self):
        return self.ht * self.tw
    @property
    def awy(self):
        return 2 * self.bf * self.tf

    #CORTANTE EM X
    def calcularParEsbVrdx(self):
        return self.h / self.tw

    def definirKvVrdx(self, a = None):
        if a is None or a/self.h > 3 or a/self.h > (260 / (self.h / self.tw))**2:
            kv = 5
        else:
            kv = 5 + 5 / (a / self.h)**2
        return kv

    #CORTANTE EM Y
    def calcularParEsbVrdy(self):
        return (self.bf/2) / self.tf

    def definirKvVrdy(self, a):
        return 1.2


    #MOMENTO EM X
    #------------
    #Estado Limite FLT
    def calcularParEsbMrdxFLT(self, lb):
        return lb / self.ry

    def calcularParEsbLimPlasMrdxFLT(self):
        return 1.76 * sqrt(self.material.e / self.material.fy)

    def calcularParEsbLimEscMrdxFLT(self):
        iy = self.iy
        j = self.j
        wx = self.wx
        e = self.material.e
        ry = self.ty
        cw = self.cw
        fy = self.fy

        beta_1 = (0.7 * fy * wx) / (e * j)
        par_esb_lim_esc = (1.38 * sqrt(iy*j) / (ry*j*beta_1)) * sqrt(1 + sqrt(1 + 27 * cw * beta_1**2 / iy))
        return par_esb_lim_esc

    def mrxFLT(self):
        return 0.7 * self.material.fy * self.wx

    def mcrxFLT(self, cb, lb):
        e = self.material.e
        iy = self.iy
        j = self.j
        cw = self.cw
        mcr = (cb * pi**2 * e * iy / lb**2) * sqrt(cw / iy * (1 + 0.039 * j * lb**2 / cw))
        return mcr

    def mrdxFLT(self, cb, lb, gamma_a1 = 1.1):

        par_esb = self.calcularParEsbMrdxFLT(lb)
        par_esb_lp = self.calcularParEsbLimPlasMrdxFLT()
        par_esb_le = self.calcularParEsbLimEscMrdxFLT()
        mr = self.mrxFLT()

        mrdx_flt = self.mplx / gamma_a1

        if par_esb_lp < par_esb < par_esb_le:
            mrdx_flt = cb*(self.mplx-(self.mplx - mr)*(par_esb - par_esb_lp)/(par_esb_le - par_esb_lp))/gamma_a1
        elif par_esb > par_esb_le:
            mrdx_flt = self.mcrxFLT() / gamma_a1
        return mrdx_flt

    #Estado Limite FLM
    def calcularParEsbMrdxFLM(self):
        return (self.bf/2) / self.tf
    def calcularParEsbLimPlasMrdxFLM(self):
        return 0.38 * sqrt(self.material.e / self.material.fy)
    def calcularParEsbLimEscMrdxFLM(self):
        return 0.83 * sqrt(self.material.e / (self.material.fy*0.7))
    def mrxFLM(self):
        return 0.7 * self.material.fy * self.wx
    def mcrxFLM(self, par_esb):
        return 0.69 * self.material.e * self.wx / par_esb ** 2

    def mrdxFLM(self, gamma_a1 = 1.1):

        par_esb = self.calcularParEsbMrdxFLM()
        par_esb_lp = self.calcularParEsbLimPlasMrdxFLM()
        par_esb_le = self.calcularParEsbLimEscMrdxFLM()
        mr = self.mrxFLM()

        mrdx_flm = self.mplx / gamma_a1

        if par_esb_lp < par_esb < par_esb_le:
            mrdx_flm = (self.mplx-(self.mplx - mr)*(par_esb - par_esb_lp)/(par_esb_le - par_esb_lp))/gamma_a1
        elif par_esb > par_esb_le:
            mrdx_flm = self.mcrxFLM(par_esb) / gamma_a1
        return mrdx_flm

    #Estado Limite FLA
    def calcularParEsbMrdxFLA(self):
        return self.h / self.tw
    def calcularParEsbLimPlasMrdxFLA(self):
        return 0.38 * sqrt(self.material.e / self.material.fy)
    def calcularParEsbLimEscMrdxFLA(self):
        return 0.5 * sqrt(self.material.e / self.material.fy)
    def mrxFLA(self):
        return self.material.fy * self.wx
    def mrdxFLA(self, gamma_a1 = 1.1):

        par_esb = self.calcularParEsbMrdxFLA()
        par_esb_lp = self.calcularParEsbLimPlasMrdxFLA()
        par_esb_le = self.calcularParEsbLimEscMrdxFLA()
        mr = self.mrxFLA()

        mrdx_fla = self.mpx / gamma_a1

        if par_esb_lp < par_esb < par_esb_le:
            mrdx_fla= (self.mplx-(self.mplx - mr)*(par_esb - par_esb_lp)/(par_esb_le - par_esb_lp))/gamma_a1
        return mrdx_fla

     #MOMENTO EM Y
    #------------
    #Estado Limite FLT
    def mrdyFLT(self, cb, lb, gamma_a1):
        return self.mply / gamma_a1

    #Estado Limite FLM
    def calcularParEsbMrdyFLM(self):
        return (self.bf/2) / self.tf
    def calcularParEsbLimPlasMrdyFLM(self):
        return 0.38 *  sqrt(self.material.e / self.material.fy)
    def calcularParEsbLimEscMrdyFLM(self):
        return 0.83 * sqrt(self.material.e / (0.7 * self.material.fy))
    def mryFLM(self):
        return 0.70 * self.fy * self.wy
    def mcryFLM(self, par_esb):
        return 0.69 * self.material.e / par_esb**2 * self.wy

    def mrdyFLM(self, gamma_a1 = 1.1):

        par_esb = self.calcularParEsbMrdyFLM()
        par_esb_lp = self.calcularParEsbLimPlasMrdyFLM()
        par_esb_le = self.calcularParEsbLimEscMrdyFLM()
        mr = self.mryFLM()

        mrdy_flt = self.mpy / gamma_a1

        if par_esb_lp < par_esb < par_esb_le:
            mrdy_flt = (self.mply-(self.mply - mr)*(par_esb - par_esb_lp)/(par_esb_le - par_esb_lp))/gamma_a1

        elif par_esb > par_esb_le:
            mrdy_flt = self.mcryFLM(par_esb) / gamma_a1

        return mrdy_flt

    #Estado Limite FLA
    def mrdyFLA(self, gamma_a1 = 1.1):
        return self.mply / gamma_a1