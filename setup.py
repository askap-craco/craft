import os
import re
import glob

from setuptools import find_packages, setup

regexp = re.compile(r'.*__version__ = [\'\"](.*?)[\'\"]', re.S)

base_package = 'craft'
base_path = os.path.dirname(__file__)

init_file = os.path.join(base_path, 'src', 'craft', '__init__.py')
with open(init_file, 'r') as f:
    module_content = f.read()

    match = regexp.match(module_content)
    if match:
        version = match.group(1)
    else:
        raise RuntimeError(
            'Cannot find __version__ in {}'.format(init_file))

with open('README.rst', 'r') as f:
    readme = f.read()

with open('CHANGELOG.rst', 'r') as f:
    changes = f.read()

def parse_requirements(filename):
    ''' Load requirements from a pip requirements file '''
    with open(filename, 'r') as fd:
        lines = []
        for line in fd:
            line.strip()
            if line and not line.startswith("#"):
                lines.append(line)
    return lines

requirements = parse_requirements('requirements.txt')


if __name__ == '__main__':
    setup(
        name='craft',
        description='Testing utilities for CRACO',
        long_description='\n\n'.join([readme, changes]),
        license='MIT license',
        url='https://github.com/strocode/craft',
        version=version,
        author='Keith Bannister',
        author_email='keith.bannister@csiro.au',
        maintainer='Keith Bannister',
        maintainer_email='keith.bannister@csiro.au',
        install_requires=requirements,
        keywords=['craft'],
        package_dir={'': 'src'},
        packages=find_packages('src'),
        package_data={'craft':['data/antenna_locations/*']},
        scripts=glob.glob('scripts/*'),
        entry_points = {
            'console_scripts':['uvfrbsim=craft.uvfrbsim:_main',
                               'craco_pipeline=craft.craco_pipeline:_main',
                               'craco_plan=craft.craco_plan:_main',
                               'craco_img_krnl=craft.craco_img_krnl:_main',
                               'craftcor=craft.craftcor:_main',
                               'fredfof_ics=craft.fredfof_ics:_main',
                               'fredfof=craft.fredfof:_main',
                               'plot_allbeams=craft.plot_allbeams:_main',
                               'quickcorr=craft.quickcorr:_main',
                               'vcraft2fil=craft.vcraft2fil:_main',
                               'uvfits2fil=craft.uvfits2fil:_main',
                               'antenna_locations=craft.antenna_locations:_main',
                               'filadd=craft.filadd:_main'
            ]
            },
        zip_safe=False,
        classifiers=['Development Status :: 3 - Alpha',
                     'Intended Audience :: Developers',
                     'Programming Language :: Python :: 3.6',
                     'Programming Language :: Python :: 3.7',
                     'Programming Language :: Python :: 3.8']
    )
