from distutils.core import setup, Extension

pyunigen_module = Extension(
    'pyunigen',
    sources=['src/pyunigen.cpp'],
    libraries=['unigen', 'approxmc', 'cryptominisat5'],
    language='C++', )

setup(
    name='pyunigen',
    ext_modules=[pyunigen_module]
)
