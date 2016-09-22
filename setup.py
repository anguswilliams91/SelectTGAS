from setuptools import setup

setup(
    include_package_data=True,
    name="tgas_sf",
    version="0.1",
    author="Angus Williams",
    author_email="anguswilliams91@gmail.com",
    packages=['tgas_sf'],
    package_dir={'tgas_sf':'tgas_sf'},
    package_data={
        'tgas_sf':['maps/*.dict'],
    },
    install_requires=['numpy','matplotlib','healpy']
    )