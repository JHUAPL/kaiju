# Project Kaiju
Simulation tools for space plasmas using fluid and particle models.

To build:
Create a build folder.
From inside the build folder, call cmake on the base kaiju directory.
From inside the build folder, call "make" to build the application.

## Using Python Tools

To install the Python tools, run:
```bash
pip install kaipy
```

## Documentation

For help and more information, visit the [Kaiju Wiki](https://kaiju-docs.readthedocs.io).

## Compilation Requirements

**Note:** This version requires `clawpack` to compile.  
Copy `claw.F` to `kaiju/src/rcm` before building.
