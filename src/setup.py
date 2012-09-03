from distutils.core import setup, Extension

mucal = Extension('mucal',
                  sources = ['mucalmodule.c', 'mucal.c'],
                  include_dirs = ['./'],
                  )

setup(name = "Mucal",
      version = "0.1",
      description = "",
      ext_modules = [mucal])
