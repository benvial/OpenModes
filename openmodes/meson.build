
subdir('src')

python_sources = ['integration.py', '__init__.py', 
'sources.py', 'ipython.py', 'model.py', 
'impedance.py', 'helpers.py', 'material.py', 
'array.py', 'modes.py', 'simulation.py', 
'parts.py', 'visualise.py', 'eig.py', 
'constants.py', 'basis.py', 'version.py', 
'multipole.py']

py.install_sources(
  python_sources,
  pure: false,
  subdir: 'openmodes'
)

subdir('operator')
subdir('geometry')
subdir('mesh')
subdir('external')
subdir('templates')
subdir('static')