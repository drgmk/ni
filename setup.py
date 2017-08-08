from setuptools import setup

setup(name='ni',
      version='0.1',
      description='zodi levels for nulling interferometry',
      url='http://github.com/drgmk/ni',
      author='Grant M. Kennedy',
      author_email='gkennedy@ast.cam.ac.uk',
      license='MIT',
      packages=['ni'],
      classifiers=['Programming Language :: Python :: 3'],
      install_requires=['numpy'],
      zip_safe=False)