from setuptools import setup

setup(name='dda_mie',
      version='0.1.0',
      description='First version of a DDA',
      url='https://github.com/nunodsousa/emDDA',
      author='Nuno de Sousa',
      author_email='nunodsousa@dipc.org',
      license='MIT',
      packages=['mole_dda'],
      install_requires=['termcolor', 'numpy', 'pandas', 'scipy>=1.3.1', 'mole_mie>=0.1.2', 'tensorflow'],
      zip_safe=False)