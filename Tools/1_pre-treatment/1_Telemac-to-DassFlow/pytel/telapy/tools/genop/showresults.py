# -*- coding: utf-8 -*-
"""
Print convergence messages

Author(s) : Fabrice Zaoui

Copyright EDF 2016-2017

"""


def printinfo(i, nitermax, fvalopt, indopt, nsimul):
    """
    Show basic info if verbose mode is on

    :param 'i': iteration number
    :param 'nitermax': maximum number of iterations
    :param 'fvalopt': best cost function
    :param 'indopt': best individual
    :param 'nsimul': Number od simulations

    :return: none
    """
    if i == 0:
        print('-----  genop optimization -----')
        print('-------------------------------')
        print('|     iter     |     Cost     |')
        print('-------------------------------')
    if nitermax != -1:
        print('      %3d       %e' % ((i+1), fvalopt))
    if (i+1 == nitermax) | (nitermax == -1):
        print('-------------------------------------\n')
        print('--> Maximum number of iterations: %d' % (i+1))
        print('-------------------------------------\n')
        print('--> Number of simulations: %d' % nsimul)
        print('-------------------------------------\n')
        print('--> Best solution found:')
        print('\t * cost function: %e' % fvalopt)
        print('\t * solution:\t')
        print('\t'.join('\t\t{}: {}'.format(*k) for k in enumerate(indopt)))
        print('-------------------------------------\n')
    return
