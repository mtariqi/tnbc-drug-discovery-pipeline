from setuptools import setup, find_packages

setup(
    name="tnbc-genomics-agent",
    version="1.0.0",
    author="Your Name",
    author_email="your@email.com",
    description="AI-powered RTK/nRTK redundancy analysis pipeline for Triple-Negative Breast Cancer genomics",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/tnbc-genomics-agent",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "anthropic>=0.28.0",
        "pandas>=2.0.0",
        "numpy>=1.26.0",
        "biopython>=1.83",
    ],
    extras_require={
        "dev": ["pytest>=8.0.0", "pytest-cov>=5.0.0"],
        "notebook": ["jupyter>=1.0.0", "ipykernel>=6.29.0"],
    },
    entry_points={
        "console_scripts": [
            "tnbc-pipeline=pipeline:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)
