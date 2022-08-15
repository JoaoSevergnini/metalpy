"""
Este módulo apresenta os métodos de verificação da capacidade resistênte de
perfis metálicos, este módulo contempla a norma brasileira ABNT:NBR8800 -
"""

from metalpy.normas._aisc import AISC360
from metalpy.normas._nbr8800 import NBR8800

__all__ = (
    'AISC360',
    'NBR8800',
)