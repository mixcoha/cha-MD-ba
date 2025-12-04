from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read().splitlines()

setup(
    name="cha-md-ba",
    version="0.1.0",
    author="Edgar Mixcoha",
    author_email="your.email@example.com",
    description="Molecular Dynamics Simulation Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/cha-md-ba",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=4.0",
            "black>=21.0",
            "isort>=5.0",
            "mypy>=0.9",
            "flake8>=3.9",
            "pre-commit>=2.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "cha-md-ba=cha_md_ba.cli.main:cli",
        ],
    },
)