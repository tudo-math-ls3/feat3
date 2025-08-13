import argparse
import os
import sys
import re
import platform
from functools import reduce

def read_csv_entry(line : str) -> list[str]:
  tmp = []
  for l in line.split(','):
    if l:
      tmp.append(l)
  return tmp

def parse_likwid_csv(likwid_output : str) -> dict:
  return_dict = {}
  for m in re.finditer(r'(STRUCT|TABLE),(.+?),(\d+),*\n?(.*?)(?=$|STRUCT|TABLE|NUMA domains)', likwid_output, flags=re.DOTALL):
    if(m.group(1) == "STRUCT"):
      loc_dict = {}
      loc_name = m.group(2)
      size = int(m.group(3))
      for sub_m in re.finditer(r'^(.+?):((?:,[^,]+)+),*$', m.group(4), flags=re.MULTILINE):
        loc_dict[sub_m.group(1)] = sub_m.group(2).lstrip(',')
      if(len(loc_dict) != size):
        print(f'Something went wrong in {loc_name}')
        print(m.group(4))
      loc_dict['type'] = "STRUCT"
      return_dict[loc_name] = loc_dict
    elif(m.group(1) == "TABLE"):
      loc_dict = {}
      loc_name = m.group(2)
      size = int(m.group(3))
      table = m.group(4).splitlines()
      list_names = table[0].rstrip(',').split(',')
      list_of_values = []
      for i in range(len(list_names)):
        list_of_values.append(['def' for _ in range(size)])
      for k,l in enumerate(table[1:]):
        tmp = l.rstrip(',').split(',')
        for i, v in enumerate(tmp):
          list_of_values[i][k] = v
      for i, name in enumerate(list_names):
        loc_dict[name] = list_of_values[i]
      loc_dict['type'] = "TABLE"
      return_dict[loc_name] = loc_dict
    else:
      print(f'No known case {m.group(1)}')
  ## sadly, there appears to be an error in the csv output of likwid, so we grap for numa domains explicit
  m = re.search(r'^NUMA domains:,(\d+?),*$', likwid_output, flags=re.MULTILINE)
  return_dict['NUMA domains'] = m.group(1)

  return return_dict

def construct_mappings(likwid_output : str) -> dict:
  likwid_dict = parse_likwid_csv(likwid_output)
  ## parse general information
  cores_per_socket = int(likwid_dict['Hardware Thread Topology']['Cores per socket'])
  sockets = int(likwid_dict['Hardware Thread Topology']['Sockets'][0])
  threads_per_core = int(likwid_dict['Hardware Thread Topology']['Threads per core'])
  print(f'Likwid Hardware Info: Sockets {sockets}, Cores Per Socket {cores_per_socket}, Threads per Core {threads_per_core}')
  ## construct our mapping hwthread -> cpu
  hw_to_cpu = [int(cpu_id) for cpu_id in likwid_dict['Topology']['Core']]
  ## sanity check
  if(len(hw_to_cpu) != cores_per_socket * sockets * threads_per_core):
    print("Number of hwthreads does not fit hardware information")
  ## also, pairing cpu -> hwthreads could be practical
  cpu_to_hw = []
  for i in range(cores_per_socket*sockets):
    cpu_to_hw.append([int(k) for k, x in enumerate(hw_to_cpu) if x == i])
  ## sanity check
  threads_per = reduce(lambda x,y: max(x, len(y)) ,cpu_to_hw, 0)
  if(threads_per != threads_per_core):
    print("Error between hardware info and parsed mapping")
  ## now construct list of numa -> cores, by first constructing numa -> hwthread, mapping to cores and then removing duplicate values of core ids
  num_numas = int(likwid_dict['NUMA domains'])
  numa_to_cores = []
  for i in range(num_numas):
    keyname = f"NUMA Topology {i}"
    loc_map = [hw_to_cpu[int(hw_thread)] for hw_thread in likwid_dict[keyname]['Processors'].split(',')]
    numa_to_cores.append(list(set(loc_map)))

  loc_dict = {'numas'    : numa_to_cores,
              'cpu_to_hw': cpu_to_hw,
              'cores_per_socket': cores_per_socket,
              'sockets'  : sockets,
              'threads_per_core' : threads_per_core,
              'num_numas'        : num_numas
              }

  return loc_dict

def construct_mapping_hex(numa_dict : dict, target_ranks_per_node : int, cull_hw_threads : bool, mpi_exec : str) -> str:
  if target_ranks_per_node < 0:
    target_ranks_per_node = numa_dict['num_numas']
  ## sanity check, is ranks per node a multiple of number of numa domains?
  if(target_ranks_per_node%numa_dict['num_numas'] != 0):
    print(f"Target ranks per node {target_ranks_per_node} is not a multiple of number of numa domains {numa_dict['num_numas']}")
    sys.exit(-1)
  if(target_ranks_per_node>len(numa_dict['cpu_to_hw'])):
    print(f"More target ranks per node {target_ranks_per_node} than pyhsical cores {len(numa_dict['cpu_to_hw'])}")
    sys.exit(-1)

  ## now construction is straight forward, for each numa domain we simply construct target_per/numa_domains mappings
  ranks_per_numa = target_ranks_per_node//numa_dict['num_numas']
  cores_per_rank = len(numa_dict['numas'][0])//ranks_per_numa
  print(f"Ranks per numa: {ranks_per_numa}, Cores per rank {cores_per_rank}")
  hex_string = ""
  cpu_to_hw = numa_dict['cpu_to_hw']
  for i in range(numa_dict['num_numas']):
    cur_numa_to_core = numa_dict['numas'][i]
    for k in range(0, len(cur_numa_to_core), cores_per_rank):
      ## get our target cores
      cur_slice = cur_numa_to_core[k:(k+cores_per_rank)]
      hex_val = 0
      for val in cur_slice:
        if(cull_hw_threads):
          # print(f" val {val} , cpuid {cpu_to_hw[val][0]}")
          hex_val |= 2**(cpu_to_hw[val][0])
        else:
          for values in cpu_to_hw[val]:
            hex_val |= 2**(values)
      # output_string += f",{hex(hex_val)}"
      hex_string += f",{hex(hex_val)}"

  hex_string = hex_string.lstrip(',')
  output_string = ""
  num_thread_per_rank = cores_per_rank if cull_hw_threads else cores_per_rank*numa_dict['threads_per_core']
  if(mpi_exec == "srun"):
    output_string = f" Plain Hexstring: {hex_string}\n Number of processes per node: {target_ranks_per_node}\n Number of threads per rank {num_thread_per_rank}\n" \
                    + f" export OMP_NUM_THREADS={cores_per_rank if cull_hw_threads else cores_per_rank*numa_dict['threads_per_core']}\n" \
                    + f" export OMP_PROC_BIND=true\n" \
                    + f" srun --tasks-per-node {target_ranks_per_node} --cpu-bind=verbose,mask_cpu:{hex_string}"
  elif(mpi_exec == "mpirun"):
    print("mpirun does not support a relative hexmap for processor ids.\n"
          "We instead have to write a rank file to do explicit multithreaded binding.\n"
          "ToDo:, write rank file for this case")
    sys.exit(-1)
  else:
    print(f"Not implemented {mpi_exec} yet.")
    sys.exit(-1)

  return output_string



if __name__ == '__main__':
  description = (
    "This tool is intended to be used in conjunction with likwid-topology to extract numa domain information and\n"
    "hyperthread information to explicitly construct a pinning mask to be used with srun or mpirun.\n"
    "Additional provides the number of threads per rank you can start the program with.\n"
    "For now, you are required to provide a file in which the output of likwid-topology -O was written.\n"
    "WARNING: You HAVE to use the -O option, else the parsing fails."
  )
  parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("--likwid-file", required=True, help="Path to the file the output of likwid-topology -O was written.")
  parser.add_argument("--ranks-per-node", default="as-numa", help="The number of mpi processes per node you want to use.\n"
                       "For now MUST be a multiple of the number of numa domains reported by likwid-topology and\n"
                       "at maximum the number of physical cores on the machine.\n"
                       "Optionally, if not provided, value as-numa or negative, the number of numa domains equalls the number of ranks.")
  parser.add_argument("--mpi-exec", default="srun", help="The intended mpi executable to be used.\n"
                       "For now only mpirun and srun supported.")
  parser.add_argument("--allow-hyperthreading", action='store_true', help="The default behavior is to ignore hyperthreading, effectively cutting\n"
                      "the number of available threads to the number of pyhsical cores.\n"
                      "If you provide this option, the hyperthreads will be used in the output flag accordingly.")

  args = parser.parse_args()
  likwid_input = ""
  if(args.likwid_file):
    if not os.path.exists(args.likwid_file):
      print(f"File {args.likwid_file} does not exist.")
      sys.exit(-1)
    with open(args.likwid_file) as f:
      likwid_input = f.read()
  ranks_per_node = -1
  if(args.ranks_per_node != "as-numa"):
    ranks_per_node = int(args.ranks_per_node)
  mpi_exec = args.mpi_exec
  cull_hyperthreading = not args.allow_hyperthreading

  print("FEAT numa binding tool")
  print(f"Started with {args}")
  numa_dict = construct_mappings(likwid_input)
  map1 = construct_mapping_hex(numa_dict, ranks_per_node, cull_hyperthreading, mpi_exec)

  print(map1)