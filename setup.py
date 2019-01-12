from distutils.core import setup


# I followed this tutorial to have both the git repository matched with the pip
# repository: https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56
setup(
        name='class4gl',
        version='0.1.2',
        license='gpl-3.0',        # https://help.github.com/articles/licensing-a-repository
        description = 'a framework to investigate the dynamics of the atmospheric boundary layer weather balloons worldwide', # Give a short description
        author = 'Hendrik Wouters',                        # Type in your name
        author_email = 'hendrik.wouters@ugent.be',         # Type in your E-Mail
        url = 'https://github.com/hendrikwout/class4gl',   # Provide either the link to your github or to your website
        download_url ='https://github.com/hendrikwout/class4gl/archive/v0.1.tar.gz',
        # I explain this later on
        keywords = ['atmospheric boundary layer', 'weather balloons',
                    'land--atmosphere interactions'],   # Keywords
        packages=['class4gl'],
        # packages=find_packages(),
        install_requires=['beautifulsoup4','pyyaml','pysolar','basemap','xarray'],
        # long_description=open('README.md').read(),
        classifiers=[
                'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4
                #'Intended Audience :: Atmospheric scientists',
                #'Topic :: modelling of the atmospheric boundary layer',
                # 'License :: gpl-3.0',   
                'Programming Language :: Python :: 3',      
                # 'Programming Language :: Python :: 3.4',
                # 'Programming Language :: Python :: 3.5',
                # 'Programming Language :: Python :: 3.6',
                # 'Programming Language :: Python :: 3.7',
              ],
)
