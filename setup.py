"pyweather setup module"

desc = ("to explore the opportunities of creating large-scale "
        "weather/climate simulation in Python")

def main():

    from setuptools import setup, find_packages

    console_scripts = ["slabs=pyslabs.command:main"]
    keywords = ["Parallel I/O", "pyslabs"]

    setup(
        name="pyweather",
        version="0.1.0",
        description=desc,
        long_description=desc,
        author=author,
        author_email="youngsung.kim.act2@gmail.com",
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
        ],
        keywords=["weather", "climate", "GPU"],
        packages=find_packages(exclude=["tests"]),
        include_package_data=True,
        project_urls={
            "Bug Reports": "https://github.com/grnydawn/pyweather/issues",
            "Source": "https://github.com/grnydawn/pyweather",
        }
    )

if __name__ == '__main__':
    import multiprocessing
    multiprocessing.freeze_support()
    main()
