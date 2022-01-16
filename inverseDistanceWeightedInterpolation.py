import pandas as pd
import matplotlib.pyplot as plt

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2016 Paul Brodersen <paulbrodersen+idw@gmail.com>

# Author: Paul Brodersen <paulbrodersen+idw@gmail.com>

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
Inverse distance weighting (IDW)
--------------------------------

Compute the score of query points based on the scores of their k-nearest neighbours,
weighted by the inverse of their distances.

Example:
--------

# import idw

# 'train'
idw_tree = idw.tree(X1, z1)

# 'test'
spacing = np.linspace(-5., 5., 100)
X2 = np.meshgrid(spacing, spacing)
grid_shape = X2[0].shape
X2 = np.reshape(X2, (2, -1)).T
z2 = idw_tree(X2)


For a more complete example see demo().

"""

import numpy as np
from scipy.spatial import cKDTree

class tree(object):
    """
    Compute the score of query points based on the scores of their k-nearest neighbours,
    weighted by the inverse of their distances.

    @reference:
    https://en.wikipedia.org/wiki/Inverse_distance_weighting

    Arguments:
    ----------
        X: (N, d) ndarray
            Coordinates of N sample points in a d-dimensional space.
        z: (N,) ndarray
            Corresponding scores.
        leafsize: int (default 10)
            Leafsize of KD-tree data structure;
            should be less than 20.

    Returns:
    --------
        tree instance: object

    Example:
    --------

    # 'train'
    idw_tree = tree(X1, z1)

    # 'test'
    spacing = np.linspace(-5., 5., 100)
    X2 = np.meshgrid(spacing, spacing)
    X2 = np.reshape(X2, (2, -1)).T
    z2 = idw_tree(X2)

    See also:
    ---------
    demo()

    """
    def __init__(self, X=None, z=None, leafsize=10):
        if not X is None:
            self.tree = cKDTree(X, leafsize=leafsize )
        if not z is None:
            self.z = np.array(z)

    def fit(self, X=None, z=None, leafsize=10):
        """
        Instantiate KDtree for fast query of k-nearest neighbour distances.

        Arguments:
        ----------
            X: (N, d) ndarray
                Coordinates of N sample points in a d-dimensional space.
            z: (N,) ndarray
                Corresponding scores.
            leafsize: int (default 10)
                Leafsize of KD-tree data structure;
                should be less than 20.

        Returns:
        --------
            idw_tree instance: object

        Notes:
        -------
        Wrapper around __init__().

        """
        return self.__init__(X, z, leafsize)

    def __call__(self, X, k=3, eps=1e-6, p=2, regularize_by=1e-9):
        """
        Compute the score of query points based on the scores of their k-nearest neighbours,
        weighted by the inverse of their distances.

        Arguments:
        ----------
            X: (N, d) ndarray
                Coordinates of N query points in a d-dimensional space.

            k: int (default 6)
                Number of nearest neighbours to use.

            p: int or inf
                Which Minkowski p-norm to use.
                1 is the sum-of-absolute-values "Manhattan" distance
                2 is the usual Euclidean distance
                infinity is the maximum-coordinate-difference distance

            eps: float (default 1e-6)
                Return approximate nearest neighbors; the k-th returned value
                is guaranteed to be no further than (1+eps) times the
                distance to the real k-th nearest neighbor.

            regularise_by: float (default 1e-9)
                Regularise distances to prevent division by zero
                for sample points with the same location as query points.

        Returns:
        --------
            z: (N,) ndarray
                Corresponding scores.
        """
        self.distances, self.idx = self.tree.query(X, k, eps=eps, p=p)
        self.distances += regularize_by
        weights = self.z[self.idx.ravel()].reshape(self.idx.shape)
        mw = np.sum(weights/self.distances, axis=1) / np.sum(1./self.distances, axis=1)
        return mw

    def transform(self, X, k=6, p=2, eps=1e-6, regularize_by=1e-9):
        """
        Compute the score of query points based on the scores of their k-nearest neighbours,
        weighted by the inverse of their distances.

        Arguments:
        ----------
            X: (N, d) ndarray
                Coordinates of N query points in a d-dimensional space.

            k: int (default 6)
                Number of nearest neighbours to use.

            p: int or inf
                Which Minkowski p-norm to use.
                1 is the sum-of-absolute-values "Manhattan" distance
                2 is the usual Euclidean distance
                infinity is the maximum-coordinate-difference distance

            eps: float (default 1e-6)
                Return approximate nearest neighbors; the k-th returned value
                is guaranteed to be no further than (1+eps) times the
                distance to the real k-th nearest neighbor.

            regularise_by: float (default 1e-9)
                Regularise distances to prevent division by zero
                for sample points with the same location as query points.

        Returns:
        --------
            z: (N,) ndarray
                Corresponding scores.

        Notes:
        ------

        Wrapper around __call__().
        """
        return self.__call__(X, k, eps, p, regularize_by)

def example():
    import matplotlib.pyplot as plt

    # create sample points with structured scores
    X1 = 10 * np.random.rand(4, 2) -5

    def func(x, y):
        return np.sin(x**2 + y**2) / (x**2 + y**2)

    z1 = func(X1[:,0], X1[:,1])

    # 'train'
    idw_tree = tree(X1, z1)

    # 'test'
    spacing = np.linspace(-5., 5., 100)
    X2 = np.meshgrid(spacing, spacing)
    grid_shape = X2[0].shape

    X2 = np.reshape(X2, (2, -1)).T
    z2 = idw_tree(X2)

    # plot
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True, figsize=(10,3))
    ax1.contourf(spacing, spacing, func(*np.meshgrid(spacing, spacing)))
    ax1.set_title('Ground truth')
    ax2.scatter(X1[:,0], X1[:,1], c=z1, linewidths=0)
    ax2.set_title('Samples')
    ax3.contourf(spacing, spacing, z2.reshape(grid_shape))
    ax3.set_title('Reconstruction')
    plt.show()
    return

def example1():
    import matplotlib.pyplot as plt

    # create sample points with structured scores
    X1 = np.array([[-5,-5],[-1,-1],[1,2],[2,3],[3,4]])
    print(X1)
    #X1 = 10 * np.random.rand(4, 2) - 5
    #print(X1)
    def func(x, y):
        return np.sin(x**2 + y**2) / (x**2 + y**2)

    z1 = func(X1[:,0], X1[:,1])
    #z1 = np.array([1,2,3,4,np.nan])
    print(z1)

    # 'train'
    idw_tree = tree(X1, z1)

    # 'test'
    spacing = np.linspace(-5., 5., 10)
    X2 = np.meshgrid(spacing, spacing)
    grid_shape = X2[0].shape

    X2 = np.reshape(X2, (2, -1)).T
    z2 = idw_tree(X2)

    # plot
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True, figsize=(10,3))
    ax1.contourf(spacing, spacing, func(*np.meshgrid(spacing, spacing)))
    ax1.set_title('Ground truth')
    ax2.scatter(X1[:,0], X1[:,1], c=z1, linewidths=0)
    ax2.set_title('Samples')
    ax3.contourf(spacing, spacing, z2.reshape(grid_shape))
    ax3.set_title('Reconstruction')
    plt.show()
    return

#data = pd.read_excel('D:/python_data/climate.xlsx',index_col='time')
#print(data.iloc[0,5:8])


def songjunwu():
    # create sample points with structured scores

    X1 = np.array([[439752.337,5784565.37],[469129.083,5821941.479],[440490.526,5818941.098]])
    z1 = np.array([1,2,3])



    # 'train'
    idw_tree = tree(X1, z1)



    # 'test'
    spacingx = np.linspace(439500, 469500, 61)
    spacingy = np.linspace(5784000, 5822500, 78)
    #print(spacingx, spacingy)
    X2 = np.meshgrid(spacingx, spacingy)
    grid_shape = X2[0].shape
    X2 = np.reshape(X2, (2, -1)).T
    z2 = idw_tree(X2)

    spacingx1 = np.linspace(442500, 453000, 22)
    spacingy1 = np.linspace(5798500, 5812500, 29)
    #print(spacingx1,spacingy1)
    X21 = np.meshgrid(spacingx1, spacingy1)
    grid_shape1 = X21[0].shape
    X21 = np.reshape(X21, (2, -1)).T
    z21 = idw_tree(X21)
    #print('aaaaaaaaaaaaaaaaaaaaaaa',z21.reshape(grid_shape1))



    # plot
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True, figsize=(10,3))
    #ax1.scatter(X1[:,0], X1[:,1], c=z1, linewidths=0)
    #ax1.set_title('Samples')
    l1 = ax2.contourf(spacingx1, spacingy1, z21.reshape(grid_shape1),cmap=plt.cm.hot)
    ax2.set_title('Smaller Reconstruction')
    l2 = ax3.contourf(spacingx, spacingy, z2.reshape(grid_shape))
    ax3.set_title('Reconstruction')
    plt.colorbar(l1)
    plt.colorbar(l2)

    plt.show()

    return
def run():
    import matplotlib.pyplot as plt
    # create sample points with structured scores

    X1 = np.array([[439752.337,5784565.37],[469129.083,5821941.479],[440490.526,5818941.098]])
    z1 = [100,2,3]

    # 'train'
    idw_tree = tree(X1, z1)

    # 'test'
    spacingx = np.linspace(442500, 453000, 22)
    spacingy = np.linspace(5798500, 5812500, 29)
    X2 = np.meshgrid(spacingx, spacingy)
    grid_shape = X2[0].shape
    X2 = np.reshape(X2, (2, -1)).T
    z2 = idw_tree(X2)
    z3 = z2.reshape(grid_shape)

    # plot
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True, figsize=(10,3))
    #ax1.scatter(X1[:,0], X1[:,1], c=z1, linewidths=0)
    #ax1.set_title('Samples')
    #l1 = ax2.contourf(spacingx1, spacingy1, z21.reshape(grid_shape1),cmap=plt.cm.hot)
    #ax2.set_title('Smaller Reconstruction')
    l2 = ax3.contourf(spacingx, spacingy, z2.reshape(grid_shape))
    ax3.set_title('Reconstruction')
    #plt.colorbar(l1)
    plt.colorbar(l2)

    plt.show()

    return(z3)




#example1()

