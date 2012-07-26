#!/usr/bin/env python
# by TR

import unittest
modules = 'util stations events trace stream rf xcorr data noisexcorr' #imaging

def suite():
    test_list = []
    for module in modules.split():
        exec ('from sito.tests import test_' + module)
        test_list.append(eval('test_%s.suite()' % module))
    return unittest.TestSuite(test_list)

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
