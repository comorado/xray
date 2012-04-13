from distutils.core import setup

config = {
    'description': 'A new project',
    'author': 'Brian Mattern',
    'url': 'http://github.com/rephorm/xray',
    'author_email': 'rephorm@rephorm.com',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['xray'],
    'scripts': [],
    'name': 'xray'
    }

setup(**config)
