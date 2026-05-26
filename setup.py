from setuptools import setup, find_packages

setup(
    name="simprepper",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        # List your dependencies here, e.g., 'requests', 'pandas'
    ],
    entry_points={
        'console_scripts': [
            # This creates a command called 'run-my-tool'
            # format: 'command_name = package.module:function'
            'simprepper = simprepper.sim_prep_amber_gromacs:main',
            'simprepper-forcefields = simprepper.utils:find_forcefields',
        ],
    },
)