#include <Windows.h>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <direct.h>

using namespace std;

// path check file
#define KERNEL_PROJECT_PATH "./kernel/kernel.vc11.vcxproj"
#define BUILD_MODES_PATH "./visual_studio/vc11-build-modes.xml"

class VC11Gen
{
public:
  // build ids
  vector<string> build_ids;
  // kernel GUID string
  string kernel_guid;
  // project GUID string
  string project_guid;
  // project name
  string project_name;
  // project path relative to root dir
  string project_path;
  // root path relative to project dir
  string root_path;

  // header file list
  set<string> hpp_list;
  // source file list
  set<string> cpp_list;

public:
  VC11Gen() :
    kernel_guid(read_kernel_guid()),
    project_guid(gen_random_guid())
  {
    parse_build_modes();
  }

  void parse_build_modes()
  {
    cout << "Parsing build modes..." << endl;
    ifstream ifs(BUILD_MODES_PATH, std::ios_base::in);
    if(!ifs.is_open() || !ifs.good())
      throw string("Failed to parse '" + string(BUILD_MODES_PATH) + "'");

    string line;
    line.reserve(512);
    size_t pos0, pos1;
    while(!ifs.eof() && ifs.good())
    {
      getline(ifs, line);
      if((pos0 = line.find("'$(Configuration)'=='")) == line.npos)
        continue;
      if((pos1 = line.find_last_of("'")) == line.npos)
        break;
      build_ids.push_back(line.substr(pos0 + 21, pos1 - pos0 - 21));
    }
  }

  static string read_kernel_guid()
  {
    cout << "Checking Kernel Project File '" << KERNEL_PROJECT_PATH << "'..." << endl;
    // try to open the stream
    ifstream ifs(KERNEL_PROJECT_PATH, std::ios_base::in);
    if(ifs.is_open() && ifs.good())
    {
      string line;
      line.reserve(512);
      size_t pos0, pos1;
      while(!ifs.eof() && ifs.good())
      {
        getline(ifs, line);
        if((pos0 = line.find("<ProjectGuid>")) == line.npos)
          continue;
        if((pos1 = line.find("</ProjectGuid>")) == line.npos)
          break;
        return line.substr(pos0 + 13, pos1 - pos0 - 13);
      }
    }
    throw string("Failed to parse '" + string(KERNEL_PROJECT_PATH) + "'");
  }

  static int gen_random_int()
  {
    static int s((int)time(NULL));
    static int x(362436069);
    static int y(521288629);
    static int z(88675123);
    int t = s ^ (s << 11);
    s = x;
    x = y;
    y = z;
    return z = z ^ (z >> 19) ^ t ^ (t >> 8);
  }

  // generates a random GUID
  static string gen_random_guid()
  {
    // generate 128 random bits
    int b0(gen_random_int());
    int b1(gen_random_int());
    int b2(gen_random_int());
    int b3(gen_random_int());

    // modify b1 and b2
    b1 = (b1 & ~0xF000) | 0x4000;
    b2 = (b2 & ~0x4000) | 0x8000;

    // print to stream
    ostringstream oss;
    oss << setfill('0') << hex << uppercase;
    oss << "{";
    oss << setw(8) << b0 << "-";
    oss << setw(4) << ((b1 >> 16) & 0xFFFF) << "-";
    oss << setw(4) << (b1 & 0xFFFF) << "-";
    oss << setw(4) << (b2 & 0xFFFF) << "-";
    oss << setw(4) << ((b2 >> 16) & 0xFFFF);
    oss << setw(8) << b3;
    oss << "}";

    return oss.str();
  }

  bool gen_project_dirs(bool& existed, int& depth)
  {
    vector<string> dirs;
    depth = 0u;

    // separate paths
    size_t n0(0);
    while(n0 != project_path.npos)
    {
      // find separators
      size_t n1 = project_path.find_first_of('\\', n0);
      size_t n2 = project_path.find_first_of('/', n0);
      size_t n3(0);
      if((n1 == project_path.npos) && (n2 == project_path.npos))
      {
        dirs.push_back(project_path.substr(n0));
        break;
      }
      else if(n1 == project_path.npos)
        n3 = n2;
      else if(n2 == project_path.npos)
        n3 = n1;
      else
        n3 = min(n1, n2);

      // separate strings
      dirs.push_back(project_path.substr(n0, n3 - n0));
      n0 = n3 + 1u;
    }

    // eliminate invalid paths
    for(size_t i(0); i < dirs.size(); )
    {
      if(dirs[i].empty())
        dirs.erase(dirs.begin() + i);
      else
        ++i;
    }

    depth = int(dirs.size());
    existed = true;

    // create directories
    string path(".");
    for(size_t i(0); i < dirs.size(); ++i)
    {
      path += "/" + dirs[i];
      // try to create the directory
      if(_mkdir(path.c_str()) != 0)
      {
        if(errno != EEXIST)
        {
          cout << "ERROR: failed to create '" << path << "' !" << endl;
          return false;
        }
      }
      else
      {
        // directory created
        existed = false;
      }
    }

    // okay
    return true;
  }

  bool find_headers()
  {
    WIN32_FIND_DATA find_data;
    HANDLE find_handle;

    // find hpp files
    memset(&find_data, 0, sizeof(WIN32_FIND_DATA));
    find_handle = FindFirstFile(string("./" + project_path + "/*.hpp").c_str(), &find_data);
    if(find_handle !=  INVALID_HANDLE_VALUE)
    {
      bool found(true);
      while(found)
      {
        hpp_list.insert(find_data.cFileName);
        found = FindNextFile(find_handle, &find_data) != FALSE;
      }
      FindClose(find_handle);
    }

    if(hpp_list.empty())
      return false;

    cout << endl << "The following header files have been found:" << endl;
    set<string>::iterator it(hpp_list.begin()), jt(hpp_list.end());
    for(; it != jt; ++it)
      cout << "- " << *it << endl;

    string cmd;
    cout << endl << "Do you want these files to be included in the project file?" << endl;
    cout << "Type 'y' for yes or 'n' for no" << endl;
    cout << ">";
    cin >> cmd;
    cout << endl;

    return (cmd == "y") || (cmd == "yes");
  }

  bool find_sources()
  {
    WIN32_FIND_DATA find_data;
    HANDLE find_handle;

    // find cpp files
    memset(&find_data, 0, sizeof(WIN32_FIND_DATA));
    find_handle = FindFirstFile(string("./" + project_path + "/*.cpp").c_str(), &find_data);
    if(find_handle !=  INVALID_HANDLE_VALUE)
    {
      bool found(true);
      while(found)
      {
        cpp_list.insert(find_data.cFileName);
        found = FindNextFile(find_handle, &find_data) != FALSE;
      }
      FindClose(find_handle);
    }

    if(cpp_list.empty())
      return false;

    cout << endl << "The following source files have been found:" << endl;
    set<string>::iterator it(cpp_list.begin()), jt(cpp_list.end());
    for(; it != jt; ++it)
      cout << "- " << *it << endl;

    string cmd;
    cout << endl << "Do you want these files to be included in the project file?" << endl;
    cout << "Type 'y' for yes or 'n' for no" << endl;
    cout << ">";
    cin >> cmd;
    cout << endl;

    return (cmd == "y") || (cmd == "yes");
  }

  bool write_solution()
  {
    // build solution filename
    string sln_path(project_path + "\\" + project_name + ".vc11.sln");
    cout << "Writing solution file        '" << sln_path << "'..." << endl;
    ofstream ofs(sln_path, ios_base::out|ios_base::trunc);
    if(!ofs.is_open() || !ofs.good())
    {
      cout << "ERROR: Failed to create '" + sln_path + "'" << endl;
      return false;
    }

    // write utf-8 bom
    ofs << char(-17) << char(-69) << char(-65) << endl;

    // write header
    ofs << "Microsoft Visual Studio Solution File, Format Version 11.00" << endl;
    ofs << "# Visual Studio 2012" << endl;

    // include app project
    ofs << "Project(\"{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}\") = \"" << project_name << ".vc11\", \"";
    ofs << project_name << ".vc11.vcxproj\", \"" << project_guid << "\"" << endl;
    ofs << "EndProject" << endl;

    // include kernel project
    ofs << "Project(\"{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}\") = \"kernel.vc11\", \"";
    ofs << root_path << "\\kernel\\kernel.vc11.vcxproj\", \"" << kernel_guid << "\"" << endl;
    ofs << "EndProject" << endl;

    // write solution configurations
    ofs << "Global" << endl;
    ofs << "\tGlobalSection(SolutionConfigurationPlatforms) = preSolution" << endl;
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it )
    {
      ofs << "\t\t" << *it << "|Win32 = " << *it << "|Win32" << endl;
      ofs << "\t\t" << *it << "|x64 = " << *it << "|x64" << endl;
    }
    ofs << "\tEndGlobalSection" << endl;

    // write project configurations
    ofs << "\tGlobalSection(ProjectConfigurationPlatforms) = postSolution" << endl;
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it )
    {
      ofs << "\t\t" << project_guid << "." << *it << "|Win32.ActiveCfg = " << *it << "|Win32" << endl;
      ofs << "\t\t" << project_guid << "." << *it << "|Win32.Build.0 = " << *it << "|Win32" << endl;
      ofs << "\t\t" << project_guid << "." << *it << "|x64.ActiveCfg = " << *it << "|x64" << endl;
      ofs << "\t\t" << project_guid << "." << *it << "|x64.Build.0 = " << *it << "|x64" << endl;
    }
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it )
    {
      ofs << "\t\t" << kernel_guid << "." << *it << "|Win32.ActiveCfg = " << *it << "|Win32" << endl;
      ofs << "\t\t" << kernel_guid << "." << *it << "|Win32.Build.0 = " << *it << "|Win32" << endl;
      ofs << "\t\t" << kernel_guid << "." << *it << "|x64.ActiveCfg = " << *it << "|x64" << endl;
      ofs << "\t\t" << kernel_guid << "." << *it << "|x64.Build.0 = " << *it << "|x64" << endl;
    }
    ofs << "\tEndGlobalSection" << endl;

    // write solution properties
    ofs << "\tGlobalSection(SolutionProperties) = preSolution" << endl;
    ofs << "\t\tHideSolutionNode = FALSE" << endl;
    ofs << "\tEndGlobalSection" << endl;

    // end-of-file
    ofs << "EndGlobal" << endl;
    ofs.close();

    return true;
  }

  bool write_project()
  {
    // build project filename
    string prj_path(project_path + "\\" + project_name + ".vc11.vcxproj");
    cout << "Writing project file        '" << prj_path << "'..." << endl;
    ofstream ofs(prj_path, ios_base::out|ios_base::trunc);
    if(!ofs.is_open() || !ofs.good())
    {
      cout << "ERROR: Failed to create '" + prj_path + "'" << endl;
      return false;
    }

    // write utf-8 bom
    ofs << char(-17) << char(-69) << char(-65);

    // write file header
    ofs << "<\?xml version=\"1.0\" encoding=\"utf-8\"\?>" << endl;
    ofs << "<Project DefaultTargets=\"Build\" ToolsVersion=\"4.0\" "
      "xmlns=\"http://schemas.microsoft.com/developer/msbuild/2003\">" << endl;

    // write global project properties
    ofs << "  <!-- global project properties -->" << endl;
    ofs << "  <PropertyGroup Label=\"Globals\">" << endl;
    ofs << "    <ProjectGuid>" << project_guid << "</ProjectGuid>" << endl;
    ofs << "    <FeastAppName>" << project_name << "</FeastAppName>" << endl;
    ofs << "    <FeastRootPath>$(ProjectDir)" << root_path << "</FeastRootPath>" << endl;
    ofs << "  </PropertyGroup>" << endl;

    // write common config import
    ofs << "  <!-- import common config -->" << endl;
    ofs << "  <Import Project=\"$(FeastRootPath)\\visual_studio\\vc11-common-config.xml\" />" << endl;

    // write header inclusion list
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <!-- Header File List -->" << endl;
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <ItemGroup Label=\"Header-Files\">" << endl;
    for(set<string>::iterator it(hpp_list.begin()); it != hpp_list.end(); ++it)
      ofs << "    <ClInclude Include=\"" << *it << "\" />" << endl;
    ofs << "  </ItemGroup>" << endl;

    // write source inclusion list
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <!-- Source File List -->" << endl;
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <ItemGroup Label=\"Source-Files\">" << endl;
    for(set<string>::iterator it(cpp_list.begin()); it != cpp_list.end(); ++it)
      ofs << "    <ClCompile Include=\"" << *it << "\" />" << endl;
    ofs << "  </ItemGroup>" << endl;

    // write project configurations
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <!-- Project Configurations -->" << endl;
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <ItemGroup Label=\"ProjectConfigurations\">" << endl;
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it)
    {
      ofs << "    <ProjectConfiguration Include=\"" << *it << "|Win32\">" << endl;
      ofs << "      <Configuration>" << *it << "</Configuration>" << endl;
      ofs << "      <Platform>Win32</Platform>" << endl;
      ofs << "    </ProjectConfiguration>" << endl;
      ofs << "    <ProjectConfiguration Include=\"" << *it << "|x64\">" << endl;
      ofs << "      <Configuration>" << *it << "</Configuration>" << endl;
      ofs << "      <Platform>x64</Platform>" << endl;
      ofs << "    </ProjectConfiguration>" << endl;
    }
    ofs << "  </ItemGroup>" << endl;

    // write kernel project reference
    ofs << "  <!-- kernel project reference -->" << endl;
    ofs << "  <ItemGroup>" << endl;
    ofs << "    <ProjectReference Include=\"$(FeastRootPath)\\kernel\\kernel.vc11.vcxproj\">" << endl;
    ofs << "      <Project>" << kernel_guid << "</Project>" << endl;
    ofs << "    </ProjectReference>" << endl;
    ofs << "  </ItemGroup>" << endl;

    // write app-target import
    ofs << "  <!-- import lib target -->" << endl;
    ofs << "  <Import Project=\"$(FeastRootPath)\\visual_studio\\vc11-target-app.xml\" />" << endl;

    // end-of-file
    ofs << "</Project>" << endl;
    ofs.close();

    return true;
  }

  void main()
  {
    // write kernel and project GUID
    cout << "Kernel  GUID: " << kernel_guid << endl;
    cout << "Project GUID: " << project_guid << endl;

    // read project name
    cout << "Please enter the project name:" << endl << ">";
    cin >> project_name;
    if(!cin.good())
      return;
    cout << endl;

    // read project path
    cout << endl << "Please enter the path to the project directory:" << endl;
    cout << ">";
    cin >> project_path;
    if(!cin.good())
      return;
    cout << endl;

    // generate project directories
    int path_depth(0);
    bool path_existed(false);
    if(!gen_project_dirs(path_existed, path_depth))
      return;

    // build relative kernel project path
    root_path = "..";
    for(int i(1); i < path_depth; ++i)
      root_path += "\\..";

    // build header and source lists
    if(!find_headers())
      hpp_list.clear();
    if(!find_sources())
      cpp_list.clear();

    // write solution file
    if(!write_solution())
      return;

    // write project file
    if(!write_project())
      return;

    cout << endl << "Project files for '" << project_name << "' have been written successfully" << endl;
  }
};


// application main entrypoint
int main(int /*argc*/, char** /*argv*/)
{
  // print header
  cout << endl << "FEAST Visual Studio 2012 Project File Generator" << endl << endl;

  // print warning
  cout << "WARNING: This tool is still under construction, so use with care!" << endl << endl;

  try
  {
    VC11Gen generator;
    generator.main();
  }
  catch(string& s)
  {
    cout << "ERROR: " << s << endl;
  }
  return 0;
}
