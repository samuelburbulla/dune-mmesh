try:
    from dune.packagemetadata import metaData
except ImportError:
    from packagemetadata import metaData
from skbuild import setup

mmeshVersion = '1.3.0.dev20211006'
duneVersion  = '2.8.0rc1'

metadata = metaData(duneVersion)[1]
metadata['version'] = mmeshVersion

# auto-generate pyproject.toml with duneVersion when building sdist
from skbuild.command.sdist import sdist
class mysdist(sdist):
    def run(self):
        requires = ['setuptools', 'wheel', 'scikit-build', 'cmake', 'ninja', 'requests']
        requires += metadata['install_requires']
        with open('pyproject.toml', 'w') as f:
            f.write("[build-system]\n")
            f.write("requires = ['"+"', '".join(requires)+"']\n")
            f.write("build-backend = 'setuptools.build_meta'\n")
        sdist.run(self)
metadata['cmdclass'] = {'sdist': mysdist}

setup(**metadata)
