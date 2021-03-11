// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

///////////////////////////////////////////////////////////////////////////////
//
// Checkpoint Info Tool
// =======================
//
// This is a tool which reads in the number of processes used to create the checkpointing file "filename.cp" created by the CheckpointControl class.
//
// To use this tool, just provide the checkpoint filename with its extension as first parameter.
//
// Usage: ./checkpoint_info path/to/checkpoint_file.cp
//
// \author Maximilian Esser
//
///////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<fstream>
#include<string>
#include<kernel/util/dist_file_io.hpp>

int main(int argc, char ** argv)
{
  std::size_t relevant_bytes = 32u;
  if(argc < 2)
  {
    std::cerr << "No checkpoint-file given!" << std::endl;
    std::cerr << "Usage: ./checkpoint_info path/to/checkpoint_file.cp!" << std::endl;
    return 1;
  }
  FEAT::String filename = argv[1];
  std::size_t pos = filename.rfind('.');
  FEAT::String extension = filename.substr(pos + 1);
  FEAT::String name = filename.substr(0, pos);
  if(extension != "cp")
  {
    std::cerr << "No valid checkpoint-file given, please reassure you added .cp at the end!" << std::endl;
    return 2;
  }
  std::ios_base::openmode bin = std::ifstream::in | std::ifstream::binary;
  std::ifstream file(filename.c_str(), bin);
  if (! file.is_open())
  {
     std::cerr << "Unable to open checkpoint-file " << filename << std::endl;
     return 2;
  }
  file.seekg (0, file.end);
  std::streamoff length = file.tellg();
  if(length < 0)
  {
    std::cerr << "Something went wrong while reading in checkpoint-file!" << std::endl;
    return 2;
  }
  file.seekg (0, file.beg);
  if(std::size_t(length) < relevant_bytes)
  {
    std::cerr << "File to small for a checkpoint-file!" << std::endl;
    return 3;
  }
  std::vector<char> buffer(relevant_bytes);
  file.read(buffer.data(), std::streamoff(relevant_bytes));
  FEAT::DistFileIO placeholder = FEAT::DistFileIO();
  std::uint64_t magic = placeholder.magic_combined;
  std::uint64_t * uiarray(reinterpret_cast<std::uint64_t *>(buffer.data()));
  if(uiarray[0] != magic)
  {
    std::cerr << "Encoding is not compatible!" << std::endl;
    return 4;
  }
  std::uint64_t rank_number = uiarray[2];
  std::cout << "The checkpoint-file requires " << rank_number << " processes." << std::endl;
  return 0;
}
