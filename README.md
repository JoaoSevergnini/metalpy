MetalPy
=======

**MetalPy** (/metalpai/) é uma biblioteca em Python de código aberto para 
determinação da capacidade de resistênte de perfis estruturais de aço, esta biblioteca
conta com um módulo de perfis metálicos e um módulo de verificação no qual é fornecida
métodos de verifição de perfis de aço de acordo com a norma brasileira *NBR8800* e 
a norma americana *AISC360*.

Veja o exemplo abaixo:

~~~python
from metalpy.perfis import **
from metalpy.normas import Nbr8800

perfil = PerfilILam('W150X150', 'A350', norma = 'NBR8800')
perfil.Ntrd_ESB()
~~~~


