from setuptools import setup

setup(
    name='colorstamps',
    version='0.1.0',
    description='2D colormaps for every occation',
    url='https://github.com/trygvrad/colorstamps',
    author='Trygve Ræder',
    author_email='tmara@dtu.dk',
    license='MIT',
    packages=['colorstamps'],
    install_requires=['colorspacious',
                      'numpy',
                      'matplotlib',
                      'scipy',
                      ],

    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
    ],
)
