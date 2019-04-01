#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.

import sys
# This is a very primitive script used for importing Maple generated C code into a hpp file.

# First argument is the name of this file, so strip that
args = sys.argv[1:]

if( ( (len(args) < 2) or len(args) > 3) or ("help" in " ".join(args)) or ("?" in " ".join(args)) ):
  print("Usage: import_lines [source] [destination] [output file]")
  print("If [output_file] is not present, this script will write to output.txt")
  sys.exit(1)

# source is a file that contains (few!) lines starting with a keyword and destination is a file containing some lines
# starting with the same keywords. These lines are then replaced with the lines from source with the corresponding
# keywords

source_file = open(args[0],'r')
dest_file = open(args[1],'r')

# Read output file or set default
if(len(args) == 3):
  output_file = open(args[2],'w')
else:
  output_file = open("output.txt",'w')

# Iterate over the destination file
for num, dest_line in enumerate(dest_file, 1):
  output_line = dest_line

  data = dest_line.split()
  # Make sure we do not try this on lines with zero or one items
  if(len(data) > 1):
    dest_first, dest_rest = data[0], data[1:]
    # Search for matching keyword in source_file
    for line in source_file:
      replace_data = line.split()
      # Make sure there is enough data for [keyword, value] split
      if(len(replace_data)>1):
        first, rest = replace_data[0], replace_data[2:]
	# If the keyword is right, replace set the output line to the one from the source file
        if dest_first == first:
          print "Replaced line ",  num,  " starting with",  dest_first
          output_line = str(data[0])+' '+str(data[1])+' '+str(''.join(replace_data[2:])+"\n")
    # Rewind the source file
    source_file.seek(0)

  output_file.write(output_line)

source_file.close()
dest_file.close()
output_file.close()

print("Imported lines from " + args[0] + " to content of " +  args[1] + " and wrote output to " + str(output_file))
sys.exit(0)
