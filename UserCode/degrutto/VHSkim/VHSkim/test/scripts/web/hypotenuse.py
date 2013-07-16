#!/usr/bin/env python
# file: hypotenuse.py

import sys, math

if len(sys.argv) != 3:  # the program name and the two arguments
  # stop the program and print an error message
  sys.exit("Must provide two positive numbers")

# Convert the two arguments from strings into numbers
x = float(sys.argv[1])
y = float(sys.argv[2])

print "Hypotenuse =", math.sqrt(x**2+y**2)
