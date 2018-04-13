diffGrad 
=======
**diffGrad** is a class which can be used for computing the derivatives and hessians of a function using finite difference. Many schemes have been implemented. 

Notice that the hessian computation remains incomplete and needs to be improved.

Features
------
**diffGrad** is able to 

* Compute derivatives with the following schemes 
	- forward and backward finite differences of order 1 to 5 (BDx and FDx with x={1,...,5})
	- central finite differences of order 2 to 8 (CDx with x={2,4,6,8})
* Minimize the number of computations (especially the responses at the central points is done only one time)
* Use a specific stepsize in every direction
* Generate the set of sample points which can be used externally for computing responses. These responses can be loaded by the class in order to compute the gradients.

First start
------
* unConstrained problems: `Example_unConstrained.m`


Download
------

The toolbox can be downloaded [here](https://bitbucket.org/luclaurent/diffGrad/downloads).

If you use `git`, you can clone the repository using the following command

    git clone git@bitbucket.org:luclaurent/diffGrad.git diffGrad






License ![GNU GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
----

    diffGrad - A toolbox to compute derivatives and hessians using finite differences
    Copyright (C) 2018  Luc LAURENT <luc.laurent@lecnam.net>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.