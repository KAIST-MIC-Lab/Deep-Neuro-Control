from setuptools import Extension
from setuptools import setup
import logging
import pathlib
import platform
import setuptools
import setuptools.command.build_ext
import setuptools.command.install
import shutil
import subprocess
import sys,os,re

class InstallationError(Exception): pass

major,minor,_,_,_ = sys.version_info
setupdir = pathlib.Path(__file__).resolve().parent

python_versions = [(3, 9), (3, 10), (3, 11), (3, 12), (3, 13)]
if (major,minor) not in python_versions: raise InstallationError("Unsupported python version")


class build_ext(setuptools.command.build_ext.build_ext):
    """
    Extend the default `build_ext` command replacing the extension
    building functionality with one that simply copies a pre-built
    extension module.
    """
    def build_extension(self,ext):
        tgtdir = pathlib.Path(self.build_lib).joinpath(*ext.name.split('.')[:-1])

        try: os.makedirs(tgtdir)
        except OSError: pass
        for s in ext.sources:
            logging.info("copying %s -> %s" % (s,tgtdir))
            shutil.copy(s,tgtdir)

class install(setuptools.command.install.install):
    """
    Extend the default install command, adding an additional operation
    that installs the dynamic MOSEK libraries.
    """
    libdir   = ['../../bin']
    instlibs = [('libtbb.12.8.dylib', 'libtbb.12.dylib'), ('libmosek64.11.0.dylib', 'libmosek64.11.0.dylib')]
    
    def findlib(self,lib):
        for p in self.libdir:
            f = pathlib.Path(p).joinpath(lib)
            if f.exists():
                return f
        raise InstallationError(f"Library not found: {lib}")
    
    def install_libs(self):
        mskdir = pathlib.Path(self.install_lib).joinpath('mosek')
        for lib,tgtname in [ (self.findlib(lib),tgtname) for (lib,tgtname) in self.instlibs ]:
            logging.info(f"copying {lib} -> {mskdir}")
            shutil.copy(lib,mskdir.joinpath(tgtname))
                
        for dirpath,dirnames,filenames in os.walk(mskdir):
            for fname in filenames:
                if pathlib.Path(fname).suffix == ".so":
                    self.osxfixlibs(os.path.join(dirpath,fname), self.instlibs, "@rpath/")
    
    def osxfixlibs(self,libfile, libs, prefix=''):
        """
        Replace path in the dynamic library reference in the DYLIB
        `libfile` for all libraries listed in `libs` (disregarding their
        current paths in the `libfile`).
    
        To see current library references, use `otool -L LIBNAME`.
    
        Example: If `prefix` is `@rpath/` the DYLIB `libfile` contains a
        reference to `/Users/john/libs/libtest.dylib`, and `libtest.dylib`
        is listed in libs, the reference will be changed to
        `@rpath/libtest.dylib`.
        """
        L = [ l.strip().split(' ')[0] for l in subprocess.check_output(['otool','-L',libfile]).decode('utf-8').split('\n')[1:] ]
        d = { os.path.basename(f) : f for f in L }
    
        args = []
        for l in libs:
            tgtlib = l[1]
            if tgtlib in d:
                args.extend([ '-change', d[tgtlib], prefix+tgtlib])
        if len(args) > 0:
            cmd = [ 'install_name_tool' ]
            cmd.extend(args)
            cmd.append(libfile)
            subprocess.call(cmd)
    
    def run(self):
        super().run()
        self.execute(self.install_libs, (), msg="Installing native libraries")

os.chdir(setupdir)
setup(name =        'Mosek',
      version =     '11.0.27',
      packages =    ['mosek', 'mosek.fusion', 'mosek.fusion.impl'],
      ext_modules = [ Extension("mosek._msk",sources=[r"mosek/_msk.cpython-39-darwin.so",r"mosek/_msk.cpython-310-darwin.so",r"mosek/_msk.cpython-311-darwin.so",r"mosek/_msk.cpython-312-darwin.so",r"mosek/_msk.cpython-313-darwin.so"]),Extension("mosek.fusion.impl.fragments",sources=[r"mosek/fusion/impl/fragments.cpython-39-darwin.so",r"mosek/fusion/impl/fragments.cpython-310-darwin.so",r"mosek/fusion/impl/fragments.cpython-311-darwin.so",r"mosek/fusion/impl/fragments.cpython-312-darwin.so",r"mosek/fusion/impl/fragments.cpython-313-darwin.so"]) ],
      cmdclass =    { "build_ext" : build_ext, "install" : install })
