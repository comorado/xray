from distutils.core import setup

config = {
    'description': 'A new project',
    'author': 'Brian Mattern',
    'url': 'https://github.com/comorado/xray',
    'author': 'Egor Klevak',
    'author_email': 'comorado@gmail.com',
    'version': '0.2',
    'install_requires': ['nose'],
    'packages': ['xray'],
    'scripts': [],
    'name': 'xray'
    }

setup(**config)
