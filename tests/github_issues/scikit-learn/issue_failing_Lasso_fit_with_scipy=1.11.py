"""
This file is for making a minimum failing example for the issue that the sklearn Lasso fit is not compatible with scipy 1.11 when using sparse X matrices.
"""

import numpy as np
from scipy.sparse import csc_array
from sklearn.linear_model import Lasso

if __name__ == '__main__':
    # Create sparse matrix
    X = csc_array(np.array([[1, 2, 0], [0, 0, 3], [4, 0, 5]]))  # also fails with dok_array
    y = np.array([1, 2, 3])

    lasso = Lasso(fit_intercept=True)  # fit_intercept=True gets a more intentionally but still unexpected TypeError, with False it looks like a bug.
    lasso.fit(X, y)     # Fit Lasso to show that it fails. LinearRegression and Ridge work fine.
    print('Done!')
