from setuptools import setup, Extension
from Cython.Build import cythonize

ext_modules = cythonize(
    [
        Extension(
            "cyrcantile._base",          # имя модуля
            ["cyrcantile/_base.pyx"],    # путь к исходнику
        )
    ],
    language_level=3,
)

setup(name="cyrcantile", ext_modules=ext_modules)
