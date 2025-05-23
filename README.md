## Software requirements

To compile this project, you need to have the following software installed:
- CMake (version 3.10 or higher)
- C++ compiler that supports C++23 standard

## Build instructions

After cloning the repository, navigate to the project directory and run the following commands:

1. Create a build directory:
```bash
mkdir build && cd build
```
2. Run CMake to configure the project:
```bash
cmake ..
```
3. Build the project:
```bash
make
```

This will create the `streaming-masked-superstring` executable in the
`build/src` directory.

## Usage

Compute the masked superstring with the approximate algorithm by using the
`compute` subcommand:
```bash
streaming-masked-superstring compute <input-fasta> <output-fasta>
```

Compute the exact masked superstring by using the `exact` subcommand:
```bash
streaming-masked-superstring exact <input-fasta> <output-fasta>
```

Compare and report accuracy of the approximate algorithm using the `compare`
subcomand:
```bash
streaming-masked-superstring compare <approximate-output> <exact-output>
```
