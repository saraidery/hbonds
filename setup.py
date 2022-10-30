from setuptools import setup

extra_requirements = {
    "tests": ["pytest", "coverage", "pytest-cov"],
}

setup(
    name="H-bonds",
    author="",
    description="Produce H-bond statistics for periodic cell",
    install_requires=["numpy", "scipy"],
    extras_require=extra_requirements,
)
