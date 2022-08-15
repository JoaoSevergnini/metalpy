from setuptools import setup, find_packages

setup(
    name='metalpy',
    version='0.0.1',
    description='Verificação de estruturas metálicas',
    author="João Severgnini",
    author_email="joao.a.severgnini@gmail.com",

    packages=find_packages(),

    url = 'https://github.com/JoaoSevergnini/metalpy',
    

    license='MIT',

    install_requires =['pandas']
)