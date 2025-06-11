from setuptools import setup, Extension
import pybind11
import subprocess

def get_root_flags():
    cflags = subprocess.check_output(['root-config', '--cflags']).decode('utf-8').strip()
    glibs = subprocess.check_output(['root-config', '--glibs']).decode('utf-8').strip()
    return cflags.split(), glibs.split()

cflags, glibs = get_root_flags()


ext_modules = [
    Extension(
        'Interface',
        ['Interface/Interface.cc', 'Interface/CNBDSIMClass.cc'],
        include_dirs=[pybind11.get_include(),
                      "./include/",
                      "/usr/local/include/bdsim/",
                      "/usr/local/include/bdsim/analysis/",
                      "/usr/local/include/bdsim/parser/",
                      "/usr/local/include/Geant4/"],
        library_dirs=['/opt/bdsim/bdsim-build/'],
        libraries=['bdsim', 'bdsimRootEvent'],
        extra_compile_args=cflags + ['-DPYBIND'],
        extra_link_args=glibs,
        language='c++',
    ),
]

setup(
    name='Interface',
    version='0.1',
    ext_modules=ext_modules,
)
