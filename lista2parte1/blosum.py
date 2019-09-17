#!/usr/bin/env python
# Lookup_table class to load, store, and use the BLOSUM lookup table

import sys

class InvalidPairException(Exception):
	pass

class Lookup_table:

	def __init__(self, table_filename):
		self._load_table(table_filename)


	# load table from txt file
	def _load_table(self, table_filename):
		with open(table_filename) as table_file:
			table = table_file.read()
			lines = table.strip().split('\n')

			header = lines.pop(0) # get header from first line
			columns = header.split()
			table = {}

			for row in lines:
				entries = row.split()
				row_name = entries.pop(0) # ignore column header
				table[row_name] = {}

				if len(entries) != len(columns):
					raise Exception('Improper entry number in row')
				for column_name in columns:
					table[row_name][column_name] = entries.pop(0)

			self._table = table


	# returns the score on the position (a, b) on the table
	def lookup_score(self, a, b):
		a = a.upper()
		b = b.upper()
	
		# check if a and b are valid
		if a not in self._table or b not in self._table[a]:
			raise InvalidPairException('[%s, %s]' % (a, b))
		return self._table[a][b]
