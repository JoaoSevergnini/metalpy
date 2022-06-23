MetalPy
=======

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
from perfis import PerfilILam
from normas import NBR8800

#Criando um perfil W150X150 a partir da classe PerfilILam
W150x150 = PerfilILam('W150X150', 'A350')

#Obtendo a resistência ao momento em relação ao eixo de maior inércia
Cb = 1
Lb = 450
NBR8800.Mrdx(W150x150, Cb, Lb)
~~~~

Ou de forma indireta, passando a norma como parâmtro de inicialização do perfil,
como demonstrado abaixo:

~~~python
from perfis import PerfilILam

W150X150 = PerfilILam('W150X150', 'A350', norma = 'NBR8800')

Cb = 1
Lb = 450
W150X150.Mrdx(Cb, Lb)
~~~
