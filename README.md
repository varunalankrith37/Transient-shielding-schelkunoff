# Transient-shielding-schelkunoff
MATLAB code for time-domain electromagnetic shielding analysis


# Time-Domain Electromagnetic Shielding Analysis

MATLAB implementation of Schelkunoff's electromagnetic shielding theory transformed to the time domain using FFT-based convolution methods.

## Description

This repository contains the MATLAB codes used in the Master's thesis:

**"Analysis of Transient Shielding Efficiency According to Schelkunoff"**

Author: Varun Alankrith Heman  
Institution: Otto von Guericke University Magdeburg  
Supervisor: Prof. Dr.-Ing. R. Vick, Dr.-Ing. M. Magdowski  
Year: 2026

## Files

- `gaussiancoppernewtest.m` - Gaussian pulse through copper shield with validation
- `DECu.m` - Double exponential pulse simulations
- `DEConcrete.m` - Double exponential pulse simulations
- `DS_Cu.m` - Damped sinusoidal pulse simulations
- `fourier.m` - Custom Fourier transform function 
- `invfourier.m` - Custom inverse Fourier transform function 

## Requirements

- MATLAB R2020a or later
- No additional toolboxes required

## Usage

Run individual .m files in MATLAB. Each script is self-contained and generates publication-quality PDF figures.

## License

MIT License - Free to use with attribution

## Citation

If you use this code, please cite the thesis:
```
V. A. Heman, "Analysis of Transient Shielding Efficiency According to Schelkunoff,"
Master's Thesis, Otto von Guericke University Magdeburg, 2026.
