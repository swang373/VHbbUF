"""
These functions emulate the overloaded TCut
operators which are not supported by PyROOT.
"""

def add(*cuts):

    n_summands = len(cuts)

    result = ''
    
    if (n_summands < 1):
        return result

    elif (n_summands == 1):
        result += '{!s}'.format(cuts[0])
        return result

    else:
        for cut in cuts[:-1]:
            result += '({!s})&&'.format(cut)
        result += '({!s})'.format(cuts[-1])
        return result
    
def mult(*cuts):

    assert (len(cuts) > 1), 'Not enough factors!'
    
    result = ''
    for cut in cuts[:-1]:
        result += '({!s})*'.format(cut)
    result += '({!s})'.format(cuts[-1])
    return result

