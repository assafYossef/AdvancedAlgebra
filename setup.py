from setuptools import setup, find_packages

# Read the requirements from requirements.txt file
with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='galwa',
    version='1.0.0',
    packages=find_packages(),
    install_requires=required,  # Use the requirements
    author='Assaf Yosef and Ido Nahum'
)