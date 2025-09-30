import sys
import os



def calculate_bounding_box(surface_file: str) -> dict:
  x_lim = [sys.float_info.max,sys.float_info.min]
  y_lim = [sys.float_info.max,sys.float_info.min]
  z_lim = [sys.float_info.max,sys.float_info.min]
  ##read in file
  first_line = True
  with open(surface_file, 'r') as fob:
    for line in fob:
      tmp = line.split()
      if not line:
        continue
      if("OFF" in line):
        continue
      if("#" == line[0]):
        continue
      if(first_line):
        first_line = False
        continue
      if(len(tmp) > 3):
        break
      x_lim[0] = min(float(tmp[0]), x_lim[0])
      x_lim[1] = max(float(tmp[0]), x_lim[1])
      y_lim[0] = min(float(tmp[1]), y_lim[0])
      y_lim[1] = max(float(tmp[1]), y_lim[1])
      z_lim[0] = min(float(tmp[2]), z_lim[0])
      z_lim[1] = max(float(tmp[2]), z_lim[1])

  return {"x_min":x_lim[0], "x_max":x_lim[1], "y_min":y_lim[0], "y_max":y_lim[1], "z_min":z_lim[0], "z_max":z_lim[1]}

def main():

  if(len(sys.argv) != 2):
    print(f"Error: Expected exactly one filename as argument, got {len(sys.argv)-1}.")
    sys.exit()
  print("Hello")
  filename = os.fsdecode(sys.argv[1])

  dic = calculate_bounding_box(filename)
  for key, val in dic.items():
    print(f"{key}:{val}")



if __name__=="__main__":
    main()