from setuptools import setup, find_packages

setup(name="afra",
      version="1.2.0",
      description="AliCPT Foreground Removal Analysis",
      author="Jiaxin Wang, Jian Yao, Le Zhang", 
      license="GPLv3",
      url="https://github.com/afra-dev/afra",
      packages=find_packages(),
      dependency_links=[],
      python_requires='>=3.5',
      install_requires=['numpy', 'scipy', 'healpy', 'pymaster', 'dynesty', 'emcee'],
      zip_safe=False,
      classifiers=["Development Status :: 4 - Beta",
                   "Topic :: Utilities",
                   "License :: OSI Approved :: GNU General Public License v3 "
                   "or later (GPLv3+)"],)
