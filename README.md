# LPN and Regular LPN Estimator

## Overview

This repository contains implementations of algorithms to estimate the computational costs for solving LPN (Learning Parity with Noise) and Regular LPN problems using the RP (Reduce and Prange's) algorithm and its faster variant, Quick-RP.

## Files

- **estimator (Reduce and Prange).sage**: Script for analyzing the costs associated with solving LPN and Regular LPN using the RP algorithm and RSD for algebraic_RP algorithm..
- **Quick-RP.sage**: Script for faster estimations using the Quick-RP method.

## Installation

You need to have SageMath installed to run these scripts. It can be downloaded from the official [SageMath website](https://www.sagemath.org/download.html).

## Usage

### Estimating Costs with RP Algorithm

For LPN: Esti(m, n, t), 
For Regular LPN: Esti_for_regular_LPN(m, n, t)

### Estimating Costs with Quick-RP Algorithm

For LPN: Quick_RP(m, n, t), 
For Regular LPN: Quick_regular_RP(m, n, t)

### Estimating Costs with Algebraic-RP Algorithm

For RSD over binary field: RSD_for_high_degree_F2(m, n, t), 
For RSD : RSD_for_high_degree(m, n, t)



### Parameters

- `m`: Number of samples
- `n`: Dimension of the secret
- `t`: Number of nonzero entries

