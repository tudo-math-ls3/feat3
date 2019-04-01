# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.

# This code was taken from https://stackoverflow.com/questions/3812849/how-to-check-whether-a-directory-is-a-sub-directory-of-another-directory/17624617#17624617
import sys
import os

def os_path_split_asunder(path):
    """
    http://stackoverflow.com/a/4580931/171094
    """
    parts = []
    while True:
        newpath, tail = os.path.split(path)
        if newpath == path:
            assert not tail
            if path: parts.append(path)
            break
        parts.append(tail)
        path = newpath
    parts.reverse()
    return parts

def is_subdirectory(potential_subdirectory, expected_parent_directory):
    """
    Is the first argument a sub-directory of the second argument?

    :param potential_subdirectory:
    :param expected_parent_directory:
    :return: True if the potential_subdirectory is a child of the expected parent directory

    >>> is_subdirectory('/var/test2', '/var/test')
    False
    >>> is_subdirectory('/var/test', '/var/test2')
    False
    >>> is_subdirectory('var/test2', 'var/test')
    False
    >>> is_subdirectory('var/test', 'var/test2')
    False
    >>> is_subdirectory('/var/test/sub', '/var/test')
    True
    >>> is_subdirectory('/var/test', '/var/test/sub')
    False
    >>> is_subdirectory('var/test/sub', 'var/test')
    True
    >>> is_subdirectory('var/test', 'var/test')
    True
    >>> is_subdirectory('var/test', 'var/test/fake_sub/..')
    True
    >>> is_subdirectory('var/test/sub/sub2/sub3/../..', 'var/test')
    True
    >>> is_subdirectory('var/test/sub', 'var/test/fake_sub/..')
    True
    >>> is_subdirectory('var/test', 'var/test/sub')
    False
    """
    def _get_normalized_parts(path):
        return os_path_split_asunder(os.path.realpath(os.path.abspath(os.path.normpath(path))))

    # Make absolute and handle symbolic links, split into components
    sub_parts = _get_normalized_parts(potential_subdirectory)
    parent_parts = _get_normalized_parts(expected_parent_directory)

    if len(parent_parts) > len(sub_parts):
        # A parent directory never has more path segments than its child
        return False

    # We expect the zip to end with the short path, which we know to be the parent
    return all(part1==part2 for part1, part2 in zip(sub_parts, parent_parts))
