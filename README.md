# LPN and Regular LPN Estimator

## Overview

This repository contains implementations of algorithms to estimate the computational costs for solving LPN (Learning Parity with Noise) and Regular LPN problems using the RP (Reduce and Prange's) algorithm and its faster variant, Quick-RP.

## Files

- **estimator (Reduce and Prange).sage**: Script for analyzing the costs associated with solving LPN and Regular LPN using the RP algorithm.
- **Quick-RP.sage**: Script for faster estimations using the Quick-RP method.

## Installation

You need to have SageMath installed to run these scripts. It can be downloaded from the official [SageMath website](https://www.sagemath.org/download.html).

## Usage

### Estimating Costs with RP Algorithm

For LPN: RP(m, n, t)
For Regular LPN: regular-RP(m, n, t)

### Estimating Costs with Quick-RP Algorithm

For LPN: Esti_RP_fixed_level3(m, t, n)
For Regular LPN: Esti_regularRP_fixed_level3(m, t, n)


### Parameters

- `m`: Number of samples
- `n`: Dimension of the secret
- `t`: Number of nonzero entries

