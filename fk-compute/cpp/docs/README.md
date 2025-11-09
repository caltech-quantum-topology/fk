# FK Computation C++ Library Documentation

This directory contains the complete Sphinx-based documentation for the FK Computation C++ Library.

## Documentation Structure

The documentation is organized as follows:

- **index.rst** - Main documentation entry point
- **introduction.rst** - Mathematical background and library overview
- **getting_started.rst** - Installation and setup instructions
- **user_guide.rst** - Comprehensive usage guide with examples
- **examples.rst** - Practical code examples and use cases
- **build_instructions.rst** - Detailed compilation and build instructions

## Building the Documentation

### Prerequisites

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Required system packages:
- Python 3.7+
- Sphinx 4.0+
- Doxygen (for C++ API documentation)

### Build Commands

```bash
# Build HTML documentation
make html

# Clean previous builds
make clean

# Build and serve locally
make html && python3 -m http.server 8000 -d _build/html
```

The generated documentation will be available in `_build/html/index.html`.

## Documentation Features

### Current Features

- **Comprehensive Guides**: Complete installation, user guide, and build instructions
- **Rich Examples**: Practical code examples for all major functionality
- **Mathematical Background**: Detailed explanation of FK invariants and algorithms
- **Multiple Output Formats**: HTML (primary), with support for PDF and EPUB
- **Clean Theme**: Uses sphinx-rtd-theme for professional appearance

### Planned Features (Future)

- **Automated API Documentation**: Using Doxygen + Breathe + Exhale for automatic C++ API docs
- **Interactive Examples**: Jupyter notebook integration for live examples
- **Search Integration**: Enhanced search capabilities
- **Multi-language Support**: Documentation translations

## Quick Start

To build and view the documentation:

```bash
cd docs
pip install -r requirements.txt
make html
python3 -m http.server 8000 -d _build/html
# Open http://localhost:8000 in your browser
```

## Content Overview

The documentation provides complete coverage of:

- Mathematical foundations and algorithmic approaches
- Installation and build procedures for all platforms
- Comprehensive API usage with practical examples
- Performance optimization and best practices
- Integration patterns and advanced topics

## Support

The documentation system is designed to be comprehensive, maintainable, and user-friendly. It provides complete coverage of the FK Computation C++ Library while being easily extensible for future development.