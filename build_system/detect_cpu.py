import platform
from build_system.feast_util import get_output

def detect_cpu():
  cputype = "unknown"

# detect arm architecture
  if "arm" in platform.machine():
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
        if cpu_part == "0xc07":
          cputype = "cortexa7"
        elif cpu_part == "0xc0f":
          cputype = "cortexa15"

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
      values = line.split("\t: ")
      if len(values) == 2:
        d[values[0].strip()] = values[1].strip()

    vendor_id = d["vendor_id"]
    cpu_family = int(d["cpu family"])
    model = int(d["model"])

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

  # TODO insert sparc support here once it has been implemented properly

  print ("Warning: cputype unknown - vendor_id: " + vendor_id + ", cpu_family: " + str(cpu_family) + ", model: " + str(model))
  return cputype
