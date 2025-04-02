from setuptools import setup, find_packages

setup(
    name="cha-md-ba",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "biopython>=1.79",
        "mdanalysis>=2.0.0",
        "pytest>=6.2.5",
        "python-dotenv>=0.19.0",
        "click>=8.0.1",
        "rich>=10.12.0",
    ],
    entry_points={
        "console_scripts": [
            "cha-md=cha_md_ba.cli:main",
            "cha-md-prepare=cha_md_ba.prepare:main",
            "cha-md-run=cha_md_ba.run:main",
            "cha-md-analyze=cha_md_ba.analyze:main",
        ],
    },
    author="Edgar Mixcoha",
    author_email="edgarmixcoha@gmail.com",
    description="Herramientas para automatizar simulaciones de dinámica molecular",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/mixcoha/cha-MD-ba",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
) 