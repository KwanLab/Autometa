from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = { }
ext_modules = [ ]

if use_cython:
    ext_modules += [
        Extension("lca_functions", [ "lca_functions.pyx" ]),
    ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules += [
        Extension("lca_functions", [ "lca_functions.c" ]),
    ]

setup(
    name='lca_functions',
    cmdclass = cmdclass,
    ext_modules=ext_modules,
)
