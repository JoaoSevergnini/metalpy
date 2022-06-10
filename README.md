MetalPy
=======

**MetalPy** (/metalpai/) é uma biblioteca em Python de código aberto para 
determinação da capacidade de resistênte de perfis estruturais de aço. Esta biblioteca
conta com um módulo de perfis metálicos e um módulo de verificação no qual é fornecida
métodos de verifição da capacidade resistente dos perfis.
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
from metalpy.perfis import **
from metalpy.normas import Nbr8800

#Criando um perfil W150X150 a partir da classe PerfilILam
W150x150 = PerfilILam('W150X150', 'A350')

#Obtendo a resistência ao momento em relação ao eixo de menor inércia
Cb = 1
Lb = 450
Nbr8800.Mrdx(W150x150, Cb, Lb)
~~~~

Ou de forma indireta, passando a norma como parâmtro de inicialização do perfil,
como demostrado abaixo:

~~~python
W150X150 = PerfilILam('W150X150', 'A350', norma = 'Nbr8800')

W150X150.Mrdx(Cb, Lb)
~~~
