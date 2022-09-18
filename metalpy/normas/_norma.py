class Norma:

    """
    Esta é uma classe abastrata que deve servir como base para implementação
    de outras normativas de verificação
    """

    def Ntrd_brt():
        """Calcula a força de tração resistente do perfil de acordo com o ELU
            de escoamento da seção bruta"""
         #Dever ser implementada na classe concreta"
        NotImplementedError

    def Ncrd():
        """Calcula a força de compressão resistente do perfil"""
         #Dever ser implementada na classe concreta"
        NotImplementedError

    def Vrdx():
        """Calcula a força cortante resistente do perfil em relação ao eixo X"""
        #Dever ser implementada na classe concreta"
        NotImplementedError
    
    def Vrdy():
        """Calcula a força cortante resistente do perfil em relação ao eixo Y"""
        #Dever ser implementada na classe concreta"
        NotImplementedError
    
    def Mrdx():
        """Calcula o momento fletor resistente do perfil em relação ao eixo X"""
        #Dever ser implementada na classe concreta"
        NotImplementedError
    
    def Mrdy():
        """Calcula o momento fletor resistente do perfil em relação ao eixo X"""
        #Dever ser implementada na classe concreta"
        NotImplementedError
    
    def Trd():
        """Calcula o momento torsor resistente do perfil em relação ao eixo X"""
        #Dever ser implementada na classe concreta"
        NotImplementedError

    def verif_NM():
        """Verifica a combinação de esforços axiais e momentos fletores"""
        #Dever ser implementada na classe concreta"
        NotImplementedError
    
