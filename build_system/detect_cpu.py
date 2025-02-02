# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
import platform
from build_system.feat_util import get_output
from build_system.feat_util import is_found

def detect_cpu():
  cputype = "unknown"

  # detect arm architecture
  if "arm" in platform.machine() or "aarch64" in platform.machine():
    if platform.system() == "Linux":
      d = {}
      f = open("/proc/cpuinfo")
      input = f.read()
      lines = input.split("\n")

      for line in lines:
        values = line.split("\t: ")
        if len(values) == 2:
          d[values[0].strip()] = values[1].strip()
        else:
          values = line.split(":")
          if len(values) == 2:
            d[values[0].strip()] = values[1].strip()

      cpu_architecture = d["CPU architecture"]
      cpu_implementer = d["CPU implementer"]
      cpu_part = d["CPU part"]

      # ARM
      if cpu_implementer == "0x41":
        if cpu_part == "0xc07":   # Raspberry Pi 2
          cputype = "cortexa7"
        elif cpu_part == "0xc0f":
          cputype = "cortexa15"
        elif cpu_part == "0xd03": # Raspberry Pi 3
          cputype = "cortexa53"
        elif cpu_part == "0xd0b": # Raspberry Pi 5
          cputype = "cortexa76"
      if cpu_implementer == "0x43": # Cavium/Marvell
        if cpu_part == "0x0af":
          cputype = "thunderx2"
      if cpu_implementer == "0x46": # Fujitsu
        if cpu_part == "0x001":
          cputype = "a64fx"
      if cpu_implementer == "0x4e":
        if cpu_part == "0x004":
          cputype = "armv8"

    else:
      print ("detect_cpu: operating system not supported")
      return cputype

    return cputype

  # detect power architecture
  if "ppc64" in platform.machine():
    if platform.system() == "Linux":
      d = {}
      f = open("/proc/cpuinfo")
      input = f.read()
      lines = input.split("\n")

      for line in lines:
        values = line.split("\t: ")
        if len(values) == 2:
          d[values[0].strip()] = values[1].strip()
        else:
          values = line.split(":")
          if len(values) == 2:
            d[values[0].strip()] = values[1].strip()

      cpu_name = d["cpu"]

      # PowerPC
      if "POWER7" in cpu_name:
        cputype = "power7"

    else:
      print ("detect_cpu: operating system not supported")
      return cputype

    return cputype

  #no special arch detected, assume x86 from now
  # read int cpu information
  if platform.system() == "Linux":
    d = {}
    f = open("/proc/cpuinfo")
    input = f.read()
    lines = input.split("\n")

    for line in lines:
      items = line.split()
      if items[0] == "vendor_id":
        vendor_id = items[2]
      if items[0] == "cpu" and items[1] == "family":
        cpu_family = int(items[3])
      if items[0] == "model" and items[1] != "name":
        model = int(items[2])
      if items[0] == "stepping":
        stepping = int(items[2])
        #assume that all three fields have been found in this order and abort, cause we now all we need
        break

  elif platform.system() == "Darwin":
    # vendor_id
    vendor_id = get_output("sysctl -A machdep.cpu.vendor")
    vendor_id = vendor_id[0]
    vendor_id = vendor_id.split(":")
    vendor_id = vendor_id[1].strip()

    # cpu_family
    cpu_family = get_output("sysctl -A machdep.cpu.family")
    cpu_family = cpu_family[0]
    cpu_family = cpu_family.split(":")
    cpu_family = int(cpu_family[1].strip())

    # model
    model = get_output("sysctl -A machdep.cpu.model")
    model = model[0]
    model = model.split(":")
    model = int(model[1].strip())

  elif platform.system() == "Windows":
    line = platform.processor()
    values = line.split(" ")
    vendor_id = values[7]
    cpu_family = int(values[2])
    model = int(values[4])

  elif platform.system() == "FreeBSD":
    if not is_found("cpuid"):
      print ("detect_cpu: you need to install the package 'cpuid' to enable feat for cpu detection.")
      return cputype

    lines = get_output("cpuid")

    cpu_family = ""
    model = ""

    for line in lines:
      if "Vendor ID:" in line:
        vendor_id = line.split('"')[1]
      if "Family" in line and cpu_family == "":
        cpu_family = int(line.split(" ")[1])
      if "Model" in line and model == "":
        model = int(line.split(" ")[1])

  else:
    print ("detect_cpu: operating system not supported")
    return cputype

  #map cpu information to cpu string
  # INTEL
  if vendor_id == "GenuineIntel":
    if cpu_family == 4:
      cputype = "i486"
    elif cpu_family == 5:
      cputype = "pentium"
    elif cpu_family == 6:
      if model < 2:
        cputype = "pentiumpro"
      elif model < 7:
        cputype = "pentium2"
      elif model < 9:
        cputype = "pentium3"
      elif model == 9:
        cputype = "pentiumm"
      elif model < 12:
        cputype = "pentium3"
      elif model == 13:
        cputype = "pentiumm"
      elif model == 14:
        cputype = "coresolo"
      elif model == 15:
        cputype = "coreduo"
      elif model == 23:
        cputype = "penryn"
      elif model == 26:
        cputype = "nehalem"
      elif model == 37:
        cputype = "westmere"
      elif model == 42:
        cputype = "sandybridge"
      elif model == 44:
        cputype = "westmere"
      elif model == 45:
        cputype = "sandybridge"
      elif model == 58:
        cputype = "ivybridge"
      elif model == 60:
        cputype = "haswell"
      elif model == 62:
        cputype = "ivybridge"
      elif model == 63:
        cputype = "haswell"
      elif model == 69:
        cputype = "haswell"
      elif model == 79:
        cputype = "broadwell"
      elif model == 85 and stepping < 5:
        cputype = "skylake-sp"
      elif model == 85 and stepping >= 5:
        cputype = "cascadelake"
      elif model == 87:
        cputype = "knightslanding"
      elif model == 94:
        cputype = "skylake"
      elif model == 106:
        cputype = "ice-lake"
      elif model == 142:
        cputype = "kaby-lake"
      elif model == 143:
        cputype = "sapphirerapids"
      elif model == 151:
        cputype = "alder-lake"
      elif model == 158:
        cputype = "coffee-lake"
    elif cpu_family == 7:
      cputype ="itanium"
    elif cpu_family == 15:
      if model == 2:
        cputype = "pentium4m"
      elif model < 3:
        cputype = "pentium4"
      elif model < 5:
        cputype = "nocona"
    elif cpu_family == 32:
      cputype = "itanium2"

  #AMD
  if vendor_id == "AuthenticAMD":
    if cpu_family == 4:
      cputype = "amd486"
    elif cpu_family == 5:
      if model < 6:
        cputype = "k5"
      else :
        cputype = "k6"
    elif cpu_family == 6:
      if model < 6:
        cputype = "athlon"
      else:
        cputype = "athlonxp"
    elif cpu_family == 15:
      if model in [5, 33, 37]:
        cputype = "opteron"
      elif model in [35,43,75,107]:
        cputype = "athlon64x2"
      elif model in [65, 67]:
        cputype = "opteronx2"
      elif model in [72,  104]:
        cputype = "turionx2"
      else:
        cputype = "athlon64"
    elif cpu_family == 16:
      if model == 2:
        cputype = "barcelona"
      elif model == 4:
        cputype = "shanghai"
      elif model in [5, 6]:
        cputype = "athlonII"
      elif model == 8:
        cputype = "istanbul"
      elif model == 9:
        cputype = "magnycours"
      elif model == 10:
        cputype = "phenomII"
    elif cpu_family == 17:
      if model == 3:
        cputype = "athlon64x2"
    elif cpu_family == 23:
      if model == 1:
        cputype = "zen"
      elif model == 49:
        cputype = "zen2"
    elif cpu_family == 25:
      if model == 1:
        cputype = "zen3"
      if model == 17:
        cputype = "zen4"
      if model == 33:
        cputype = "zen3"
      if model == 80:
        cputype = "zen3"

  # TODO insert sparc support here once it has been implemented properly

  if "unknown" in cputype:
    print ("Warning: cputype unknown - vendor_id: " + vendor_id + ", cpu_family: " + str(cpu_family) + ", model: " + str(model))

  return cputype
