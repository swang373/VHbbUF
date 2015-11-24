"""
These functions substitute the overloaded TCut
operators which are not supported within PyROOT.
"""

def add(*cuts):

    assert (len(cuts) > 1), 'Not enough summands!'

    result = ''
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

