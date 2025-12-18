"""
setup script for ar signaling simulation package
"""

from setuptools import setup, find_packages

with open("readme.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="ar-signaling-sim",
    version="0.1.0",
    author="biosim",
    description="Comprehensive computational model of androgen receptor signaling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/androgen-receptor-signaling-sim",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "black>=22.0",
            "flake8>=5.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "ar-sim=scripts.run_simulation:main",
            "ar-dashboard=src.visualization.interactive:launch_dashboard",
        ],
    },
)
