#!/usr/bin/env python

"""Running all the PBS scripts"""

import os
import re

# Make
print('Now making')
os.system('make')
# Find all files
files = [f for f in os.listdir('.') if os.path.isfile(f)]
# Find files matching a specific pattern
string_to_find = 'PBSDriver-\d'
files = [f for f in files if re.search(string_to_find, f) != None]
# Call the file
for f in files:
    print('\n'*3 + '='*60)
    print('Running the {} script'.format(f))
    print('='*60 + '\n'*3)
    os.system('python ' + f)
