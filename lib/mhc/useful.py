"""
useful
"""

from itertools import izip


def remove_blank_fields(fields, limit):
    """
    Removes blank elements from a list of strings.  Limit refers to the number of non-blank elements in the new list; once this number is achieved, no more blank elements are removed.
    
    @param fields: the list of strings to remove blank elements from
    @param limit: the number of non-blank elements after which no more blank elements will be removed
    
    @return: a new list with blank elements removed
    """
    new_fields = []
    count = 0
    for field in fields:
        if len(new_fields)>=limit: break
        count += 1
        if len(field)==0: continue
        new_fields.append(field)
    return new_fields + fields[count:]


def argmax(*args):
    if type(args[0]) in [list, tuple]:
        array = args[0]
    else:
        array = args
    return max(izip(array, xrange(len(array))))
