#!/usr/bin/env bash
# -*- coding:  utf-8 -*-
# Author: Gaurav Sablok
# date: 2023-08-05
# MIT License
# a numpy simulator writted in bash for the quick recall
# main ✗ $ sh numpy_simulator.sh
# ask for the numpy command:options
# The following options are there
# the numpy options in this are:
#     1. arange: provides the range implementation of numpy
#     2. add provides how to apply the add acorss numpy arrays
#     3. subtract provides how to apply across the numpy arrays
#     4. multiply provides how to apply multiply across the numpy arrays
#     5.divide provides how to apply divide across the numpy arrays
#     6.reciprocal provides how to apply reciprocal across the numpy arrays
#     7.power provides how to apply power across the numpy arrays
#     8.mod provides how to apply mod or remainder across the numpy arrays
#     9.sum provides how to apply sum across the numpy arrays
#     10.diagnonalsum provides how to apply sum across diagonals across the numpy arrays
#     11. average provides how to apply average across the numpy arrays
#     12. prod provides how to apply product across the numpy arrays
# (base) gauravsablok@Gauravs-MacBook-Pro:~/Desktop/GitHub/chpc_cluster
# main ✗ $ sh numpy_simulator.sh
# ask for the numpy command:
# please use the options command to list the options

read -p "ask for the numpy command:" numpy
if [[ $numpy == "" ]]; then
    echo "please use the options command to list the options"
fi
if [[ $numpy == options ]]; then
    echo "The following options are there"
    echo "the numpy options in this are:
    1. arange: provides the range implementation of numpy
    2. add provides how to apply the add acorss numpy arrays
    3. subtract provides how to apply across the numpy arrays
    4. multiply provides how to apply multiply across the numpy arrays
    5.divide provides how to apply divide across the numpy arrays
    6.reciprocal provides how to apply reciprocal across the numpy arrays
    7.power provides how to apply power across the numpy arrays
    8.mod provides how to apply mod or remainder across the numpy arrays
    9.sum provides how to apply sum across the numpy arrays
    10.diagnonalsum provides how to apply sum across diagonals across the numpy arrays
    11. average provides how to apply average across the numpy arrays
    12. prod provides how to apply product across the numpy arrays"
fi
if [[ $numpy == arange ]]; then
    echo "(np.arange(integer, dtype="i"). 
                                 reshape(rows,columns)"
    echo "np.arange(start,stop, dtype = "i"). 
                              reshape(rows,columns)"
fi
if [[ $numpy == add ]]; then
    echo "np.add(array_1, array_2)"
fi
if [[ $numpy == subtract ]]; then
    echo "np.subtract(array_1, array_2)"
fi
if [[ $numpy == multiply ]]; then
    echo "np.multiply(array_1, array_2)"
fi
if [[ $numpy == divide ]]; then
    echo "np.divide(array_1, array_2)"
fi
if [[ $numpy == reciprocal ]]; then
    echo "this method returns the 
                        argument inverse"
    echo "np.reciprocal(np.array
                (np.arange(1,10)).reshape(3,3))"
    echo "if the value is greater 
                than 1 then it will return 0"
fi
if [[ $numpy == power ]]; then
    echo "np.power(np.arange(1,10).
            reshape(3,3),2)"
    echo "it will return the power 
            of the each element of the array"
fi
if [[ $numpy == mod ]]; then
    echo "np.mod and np.remainder are the same
        and they return the remainder of the modulus
        np.mod(np.arange(1,10), np.arange(1,10))
        np.remainder(np.mod(np.arange(1,10), np.arange(1,10)))"
fi
if [[ $numpy == sum ]]; then
    echo "this returns the sum of the array"
    echo "np.arange(1,10).sum()"
fi
if [[ $numpy == diagonalsum ]]; then
    echo "to do the diagonal sum along
        the axis you have to do the diagonal sum
        np.np.arange(1,10).reshape(3,3).sum(axis = 1)
        for sum across the horizontal axis and
        np.arange(1,10).reshape(3,3).sum(axis = 0)
        for sum across the vertical axis"
fi
if [[ $numpy == average ]]; then
    echo "np.average(np.arange(1,10).reshape(3,3))"
    echo "this functions will calculate the
    average of the array"
fi
if [[ $numpy == prod ]]; then
    echo "np.prod(np.arange(1,10).reshape(3,3))
        this will give the product of all the array elements"
    echo "np.prod(np.arange(1,10).reshape(3,3), axis = 0)
        this will give the product of the elements across the 
        horizontal axis and to gain across the vertical axis
        np.prod(np.arange(1,10).reshape(3,3), axis = 1)"
fi
