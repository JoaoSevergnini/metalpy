<h1 align="center">
<img src="logo/logo_MP.svg" width="300">
</h1><br>

**MetalPy** (/metalpai/) é uma biblioteca em Python de código aberto para 
determinação da capacidade de resistênte de perfis estruturais de aço. Esta biblioteca
conta com um módulo de perfis metálicos e um módulo de verificação no qual são fornecidos
métodos de verificação da capacidade resistente dos perfis.

O módulo de perfis metálicos conta com perfis dos tipos:

 * I laminado (W, H, HP) (Ver ...)
 * I soldados (Ver .. )
 * Seção caixão (Ver ...)
 * Tubos retangulares (Ver ...)
 * Tubos circulares (ver ... )

O módulo para verificação dos perfis de aço conta com a norma brasileira *NBR8800* e 
a norma americana *AISC360*, e possibilita a determinação da capacidade reistênte 
dos perfis de acordo com o métodos dos estados limites últimos. As funções que permitem
o determinação da capacidade resistente podem ser utlizadas de duas formas, de forma direta
como no exemplo a seguir:

~~~python
from metalpy.perfis import PerfilILam
from metalpy.normas import NBR8800
from metalpy.material import Aco

#Istância da classe Aco com propriedades em kN/cm²
AR350 = Aco(20000, 0.3, 35, 45)

#Criando uma instancia da classe PerfilILam que representa o perfil W150X150 com as propriedades em cm
W150x150 = PerfilILam('W150X150', mat = AR350, und = 'cm')

#Obtendo a resistência ao momento em relação ao eixo de maior inércia
Cb = 1
Lb = 450
NBR8800.Mrdx(W150x150, Cb, Lb)
~~~~

Ou de forma indireta, passando a norma como parâmtro de inicialização do perfil,
como demonstrado abaixo:

~~~python
from metalpy.perfis import PerfilILam
from metalpy.material import Aco

#Istância da classe Aco com propriedades em kN/cm²
AR350 = Aco(20000, 0.3, 35, 45)

#Criando uma instancia da classe PerfilILam que representa o perfil W150X150 com as propriedades em cm
W150X150 = PerfilILam('W150X150', mat = AR350, und = 'cm', norma = 'NBR8800')

Cb = 1
Lb = 450
W150X150.Mrdx(Cb, Lb)
~~~
