import re
from setuptools import setup, find_packages

def get_reqs():
    reqs = []
    for line in open("requirements.txt", "r").readlines():
        reqs.append(line)
    return reqs

def get_property(prop, project):
    result = re.search(
        r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
        open(project + "/__init__.py").read(),
    )
    return result.group(1)

setup(
    name='starsign', 
    version=get_property("__version__", "starsign"),
    description="StarSign: more fun than normal horoscopes!",
    url="https://github.com/richardson-t/starsign",
    packages=find_packages(),
    keywords="Zodiac Star Fun Horoscope",
    install_requires=get_reqs())
