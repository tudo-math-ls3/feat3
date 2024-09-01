import os
import sys

# patch_lst entry: [line-no, original, patched]
def patch_file(directory, filename, patch_lst):
  filename = os.path.join(directory, filename)
  # create backup file if it doesn't exist
  filename_b = filename + ".backup"
  if not os.path.isfile(filename_b):
    os.rename(filename, filename_b)
  print("Patching '%s'..." % filename)
  # open backup file for reading
  fi = open(filename_b, "rt")
  # open file for writing
  fo = open(filename, "wt")
  # loop over all input file lines
  lno = 0
  for line in fi:
    lno = lno + 1
    if (len(patch_lst) > 0) and (patch_lst[0][0] == lno):
      # this line is to be patched
      if line.strip() != patch_lst[0][1]:
        print("ERROR: when processing file '%s': in line %i" % (filename, lno))
        print("expected : '%s'" % patch_lst[0][1])
        print("but found: '%s'" % line.strip())
        print("Patch aborted!")
        sys.exit(1)
      # okay replace line
      fo.write(patch_lst[0][2] + "\n")
      # remove patch line
      patch_lst = patch_lst[1:]
    else:
      fo.write(line)
  # ensure that all patches were applied
  if len(patch_lst) > 0:
    print("ERROR: when processing file '%s': end of file found, but there are still patches left")
    print("Patch aborted!")
    sys.exit(1)
  # okay, that's it
  fo.close()
  fi.close()
