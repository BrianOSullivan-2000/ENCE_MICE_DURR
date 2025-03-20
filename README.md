# ENCE_MICE_DURR

This repository contains the code used for the paper:

[Infilling of high-dimensional rainfall networks through multiple imputation by chained equations](https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/joc.8513)

The paper introduces two imputation methods applied to a high-dimensional rainfall network in the Republic of Ireland. These are Elastic-Net Chained Equations (ENCE) and Mutiple Imputation by Chained Equations using Direct Use of Regularisation (MICE DURR). Both methods use a series of regularised regression models to predict missing values at each rainfall station, where MICE DURR also applies Multiple Imputation to account for bias/uncertainty that may be caused by missingness present in the data. These methods are designed with high-dimensional monthly rainfall data in mind, but they are flexible, regression-based approaches applicable to a variety of datasets.

# Citation

Please cite the following paper if you use this code:
O'Sullivan, B., & Kelly, G. (2024). 
Infilling of high-dimensional rainfall networks through multiple imputation by chained equations. 
International Journal of Climatology, 44(9), 3075–3091. 
https://doi.org/10.1002/joc.8513
