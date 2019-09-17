#!/usr/bin/env python
# Usage: python blosum.py blosum62.txt
#        Then, enter input in "row col" format -- e..g, "s f".
import sys

class InvalidPairException(Exception):
  pass

class Matrix:

  def __init__(self, matrix_filename):
    self._load_matrix(matrix_filename)


  # load matrix from txt file
  def _load_matrix(self, matrix_filename):
    with open(matrix_filename) as matrix_file:
      matrix = matrix_file.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0) # get header from first line
    columns = header.split()
    matrix = {}

    for row in lines:
      entries = row.split()
      row_name = entries.pop(0) # ignore column header
      matrix[row_name] = {}

      if len(entries) != len(columns):
        raise Exception('Improper entry number in row')
      for column_name in columns:
        matrix[row_name][column_name] = entries.pop(0)

    self._matrix = matrix



  def lookup_score(self, a, b):
    a = a.upper()
    b = b.upper()

#    print("a: " + a + " ...b: " + b) # DEBUG

    if a not in self._matrix or b not in self._matrix[a]:
      raise InvalidPairException('[%s, %s]' % (a, b))
    return self._matrix[a][b]
