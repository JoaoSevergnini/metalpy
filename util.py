Und = {'mm': 1, 'cm': 10, 'm': 10e3}


class NumPositivo:

    def __init__(self, var):
        self.nome_atr = str('_'+var)

    def __set__(self, instance, valor):
        if valor >= 0:
            setattr(instance, self.nome_atr, valor)
        else:
            raise ValueError(self.nome_atr[1:] + ' deve ser >= 0')

    def __get__(self, instance, owner):
        return getattr(instance, self.nome_atr)


class PropGeo(NumPositivo):

    def __init__(self, var: str, dim: int):
        self.dim = dim
        self.un = None

        super().__init__(var)

    def __set__(self, instance, valor):

        if hasattr(instance, self.nome_atr):
            raise (self.nome_atr[1:] + ' não pode ser alterado')

        else:
            if hasattr(instance, 'und'):
                super().__set__(instance, valor / (Und[instance.und] ** self.dim))
            else:
                super().__set__(instance, valor)


