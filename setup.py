from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        "cyrcantile._base",              # имя подмодуля
        ["cyrcantile/_base.pyx"],
    )
]

setup(
    name="cyrcantile",
    ext_modules=cythonize(extensions, language_level=3),
)
