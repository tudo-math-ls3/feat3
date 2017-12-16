/**
 * \brief FEAT Visual Studio 2015 Project File Generator
 *
 * This tool can be used to generate FEAT application/tool project files for
 * Microsoft (R) Visual Studio 2015.
 *
 * \author Peter Zajac
 */
#if !defined(_MSC_VER) || (_MSC_VER < 1900)
#error This tool can only be compiled using Visual Studio 2015 or higher.
#endif

#include <Windows.h>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <deque>
#include <direct.h>

using namespace std;

#define KERNEL_PROJECT_PATH "./build_system/vc14/kernel.vc14.vcxproj"
#define BUILD_MODES_PATH "./build_system/vc14/build-modes.xml"

class ArgList :
  public deque<string>
{
public:
  explicit ArgList(int argc, char** argv)
  {
    for(int i(1); i < argc; ++i)
      push_back(argv[i]);
  }

  string next()
  {
    string arg(front());
    pop_front();
    return arg;
  }
};

class Generator
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

  // specifies whether this is a test project
  bool testing;

  // header file list
  set<string> hpp_list;
  // source file list
  set<string> cpp_list;
  // CUDA source file list
  set<string> cu_list;

public:
  Generator() :
    kernel_guid(read_kernel_guid()),
    project_guid(gen_random_guid()),
    testing(false)
  {
    parse_build_modes();
  }

  /// parses all build mode strings from 'vc14-build-modes.xml'
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

  /// reads out the kernel project GUID
  static string read_kernel_guid()
  {
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

  /// generates a pseudo-random 32-bit int
  static int gen_random_int()
  {
    static int s((int)time(nullptr));
    static int x(362436069);
    static int y(521288629);
    static int z(88675123);
    int t = s ^ (s << 11);
    s = x;
    x = y;
    y = z;
    return z = z ^ (z >> 19) ^ t ^ (t >> 8);
  }

  /// generates a random GUID
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

  static vector<string> split_path(string path, bool nodots = true)
  {
    vector<string> dirs;
    size_t n0(0);
    while(n0 != path.npos)
    {
      // find separators
      size_t n1 = path.find_first_of('\\', n0);
      size_t n2 = path.find_first_of('/', n0);
      size_t n3(0);
      if((n1 == path.npos) && (n2 == path.npos))
        n3 = path.npos;
      else if(n1 == path.npos)
        n3 = n2;
      else if(n2 == path.npos)
        n3 = n1;
      else
        n3 = min(n1, n2);

      // extract dirname
      string dirname(n3 == path.npos ? path.substr(n0) : path.substr(n0, n3 - n0));
      n0 = (n3 == path.npos ? n3 : n3 + 1u);

      // skip empty paths
      if(dirname.empty())
        continue;

      // let's see if we have to skip dots
      if(nodots)
      {
        // skip current-dir marker '.'
        if(dirname.compare(".") == 0)
          continue;

        // if we have a parent-dir marker '..', check if the last parsed dir name (if any)
        // is not also a parent-dir marker...
        if((dirname.compare("..") == 0) && (!dirs.empty()) && (dirs.back().compare("..") != 0))
        {
          dirs.pop_back();
          continue;
        }
      }

      // push dir name
      dirs.push_back(dirname);
    }

    // return directory vector
    return std::move(dirs);
  }

  static string remove_ext(string name)
  {
    size_t n = name.find_last_of('.');
    if(n == name.npos)
      return name;
    else
      return name.substr(0, n);
  }

  static bool gen_dirs(vector<string>& dirs)
  {
    bool existed = true;
    // create directories
    string path(".");
    for(size_t i(0); i < dirs.size(); ++i)
    {
      path += "/" + dirs[i];
      // try to create the directory
      if(_mkdir(path.c_str()) != 0)
      {
        if(errno != EEXIST)
          throw string("Failed to create '" + path + "' !");
      }
      else
      {
        // directory created
        existed = false;
      }
    }

    // okay
    return existed;
  }

  /// generates the project directory tree (if it does not already exist)
  bool gen_project_dirs(int& depth)
  {
    vector<string> dirs(split_path(project_path));
    depth = int(dirs.size());
    return gen_dirs(dirs);
  }

  /// searches for *.<ext> files in the application directory
  void find_files(set<string>& lst, string ext)
  {
    WIN32_FIND_DATA find_data;
    HANDLE find_handle;

    // find hpp files
    memset(&find_data, 0, sizeof(WIN32_FIND_DATA));
    find_handle = FindFirstFile(string("./" + project_path + "/*." + ext).c_str(), &find_data);
    if(find_handle !=  INVALID_HANDLE_VALUE)
    {
      bool found(true);
      while(found)
      {
        lst.insert(find_data.cFileName);
        found = FindNextFile(find_handle, &find_data) != FALSE;
      }
      FindClose(find_handle);
    }
  }


  void find_files()
  {
    find_files(hpp_list, "hpp");
    find_files(cpp_list, "cpp");
    find_files(cu_list, "cu");
  }

  /// writes the application solution file
  void write_solution()
  {
    // build solution filename
    string sln_path(project_path + "\\" + project_name + ".vc14.sln");
    cout << "Writing solution file: '" << sln_path << "'..." << endl;
    ofstream ofs(sln_path, ios_base::out|ios_base::trunc);
    if(!ofs.is_open() || !ofs.good())
      throw string("Failed to create '" + sln_path + "'");

    // write utf-8 bom
    ofs << char(-17) << char(-69) << char(-65) << endl;

    // write header
    ofs << "Microsoft Visual Studio Solution File, Format Version 12.00" << endl;
    ofs << "# Visual Studio 14" << endl;
    ofs << "VisualStudioVersion = 14.0.23107.0" << endl;

    // include app project
    ofs << "Project(\"{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}\") = \"" << project_name << ".vc14\", \"";
    ofs << project_name << ".vc14.vcxproj\", \"" << project_guid << "\"" << endl;
    ofs << "EndProject" << endl;

    // include kernel project
    ofs << "Project(\"{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}\") = \"kernel.vc14\", \"";
    ofs << root_path << "\\build_system\\vc14\\kernel.vc14.vcxproj\", \"" << kernel_guid << "\"" << endl;
    ofs << "EndProject" << endl;

    // write solution configurations
    ofs << "Global" << endl;
    ofs << "\tGlobalSection(SolutionConfigurationPlatforms) = preSolution" << endl;
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it )
    {
      ofs << "\t\t" << *it << "|x64 = " << *it << "|x64" << endl;
      ofs << "\t\t" << *it << "|x86 = " << *it << "|x86" << endl;
    }
    ofs << "\tEndGlobalSection" << endl;

    // write project configurations
    ofs << "\tGlobalSection(ProjectConfigurationPlatforms) = postSolution" << endl;
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it )
    {
      ofs << "\t\t" << project_guid << "." << *it << "|x64.ActiveCfg = " << *it << "|x64" << endl;
      ofs << "\t\t" << project_guid << "." << *it << "|x64.Build.0 = " << *it << "|x64" << endl;
      ofs << "\t\t" << project_guid << "." << *it << "|x86.ActiveCfg = " << *it << "|Win32" << endl;
      ofs << "\t\t" << project_guid << "." << *it << "|x86.Build.0 = " << *it << "|Win32" << endl;
    }
    for(vector<string>::iterator it(build_ids.begin()); it != build_ids.end(); ++it )
    {
      ofs << "\t\t" << kernel_guid << "." << *it << "|x64.ActiveCfg = " << *it << "|x64" << endl;
      ofs << "\t\t" << kernel_guid << "." << *it << "|x64.Build.0 = " << *it << "|x64" << endl;
      ofs << "\t\t" << kernel_guid << "." << *it << "|x86.ActiveCfg = " << *it << "|Win32" << endl;
      ofs << "\t\t" << kernel_guid << "." << *it << "|x86.Build.0 = " << *it << "|Win32" << endl;
    }
    ofs << "\tEndGlobalSection" << endl;

    // write solution properties
    ofs << "\tGlobalSection(SolutionProperties) = preSolution" << endl;
    ofs << "\t\tHideSolutionNode = FALSE" << endl;
    ofs << "\tEndGlobalSection" << endl;

    // end-of-file
    ofs << "EndGlobal" << endl;
    ofs.close();
  }

  /// writes the application project file
  void write_project()
  {
    // build project filename
    string prj_path(project_path + "\\" + project_name + ".vc14.vcxproj");
    cout << "Writing project file.: '" << prj_path << "'..." << endl;
    ofstream ofs(prj_path, ios_base::out|ios_base::trunc);
    if(!ofs.is_open() || !ofs.good())
      throw string("Failed to create '" + prj_path + "'");

    // write utf-8 bom
    ofs << char(-17) << char(-69) << char(-65);

    // write file header
    ofs << "<\?xml version=\"1.0\" encoding=\"utf-8\"\?>" << endl;
    ofs << "<Project DefaultTargets=\"Build\" ToolsVersion=\"12.0\" "
      "xmlns=\"http://schemas.microsoft.com/developer/msbuild/2003\">" << endl;

    // write global project properties
    ofs << "  <!-- global project properties -->" << endl;
    ofs << "  <PropertyGroup Label=\"Globals\">" << endl;
    ofs << "    <ProjectGuid>" << project_guid << "</ProjectGuid>" << endl;
    ofs << "    <FeatAppName>" << project_name << "</FeatAppName>" << endl;
    ofs << "    <FeatRootPath>$(ProjectDir)" << root_path << "</FeatRootPath>" << endl;
    ofs << "  </PropertyGroup>" << endl;

    // write common config import
    ofs << "  <!-- import common config -->" << endl;
    ofs << "  <Import Project=\"$(FeatRootPath)\\build_system\\vc14\\common-config.xml\" />" << endl;

    // write header inclusion list
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <!-- Header File List -->" << endl;
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <ItemGroup Label=\"Header-Files\">" << endl;
    if(testing)
      ofs << "    <ClInclude Include=\"" << root_path << "\\test_system\\*.hpp\" />" << endl;
    for(set<string>::iterator it(hpp_list.begin()); it != hpp_list.end(); ++it)
      ofs << "    <ClInclude Include=\"" << *it << "\" />" << endl;
    ofs << "  </ItemGroup>" << endl;

    // write source inclusion list
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <!-- Source File List -->" << endl;
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <ItemGroup Label=\"Source-Files\">" << endl;
    if(testing)
      ofs << "    <ClCompile Include=\"" << root_path << "\\test_system\\test_system.cpp\" />" << endl;
    for(set<string>::iterator it(cpp_list.begin()); it != cpp_list.end(); ++it)
      ofs << "    <ClCompile Include=\"" << *it << "\" />" << endl;
    ofs << "  </ItemGroup>" << endl;

    // write CUDA inclusion list
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <!-- CUDA File List -->" << endl;
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <ItemGroup Label=\"CUDA-Files\">" << endl;
    for(set<string>::iterator it(cu_list.begin()); it != cu_list.end(); ++it)
      ofs << "    <CudaCompile Include=\"" << *it << "\" />" << endl;
    ofs << "  </ItemGroup>" << endl;

    // write kernel project reference
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <!-- Final Imports -->" << endl;
    ofs << "  <!-- ********************************************************************* -->" << endl;
    ofs << "  <ItemGroup>" << endl;
    ofs << "    <ProjectReference Include=\"$(FeatRootPath)\\build_system\\vc14\\kernel.vc14.vcxproj\">" << endl;
    ofs << "      <Project>" << kernel_guid << "</Project>" << endl;
    ofs << "    </ProjectReference>" << endl;
    ofs << "  </ItemGroup>" << endl;

    // write app-target import
    ofs << "  <Import Project=\"$(FeatRootPath)\\build_system\\vc14\\target-app.xml\" />" << endl;

    // end-of-file
    ofs << "</Project>" << endl;
    ofs.close();
  }

  bool generate_dirs()
  {
    // generate project directories
    int path_depth(0);
    bool path_existed(gen_project_dirs(path_depth));

    // build relative kernel project path
    root_path = "..";
    for(int i(1); i < path_depth; ++i)
      root_path += "\\..";

    return path_existed;
  }

  void generate_files()
  {
    // write kernel and project GUID
    cout << endl;
    cout << "Kernel  GUID: " << kernel_guid << endl;
    cout << "Project GUID: " << project_guid << endl;
    cout << endl;

    if (!hpp_list.empty())
    {
      cout << "The following header files have been included in the project:" << endl;
      set<string>::iterator it(hpp_list.begin()), jt(hpp_list.end());
      for (; it != jt; ++it)
        cout << "- " << *it << endl;
      cout << std::endl;
    }

    if (!cpp_list.empty())
    {
      cout << "The following source files have been included in the project:" << endl;
      set<string>::iterator it(cpp_list.begin()), jt(cpp_list.end());
      for (; it != jt; ++it)
        cout << "- " << *it << endl;
      cout << std::endl;
    }

    if (!cu_list.empty())
    {
      cout << "The following CUDA files have been included in the project:" << endl;
      set<string>::iterator it(cu_list.begin()), jt(cu_list.end());
      for (; it != jt; ++it)
        cout << "- " << *it << endl;
      cout << std::endl;
    }

    // write solution file
    write_solution();

    // write project file
    write_project();

    // okay
    cout << endl << "Project files for '" << project_name << "' have been written successfully" << endl;
  }

  void main()
  {
    // read project path
    cout << endl << "Please enter the path to the project directory:" << endl;
    cout << ">";
    cin >> project_path;
    if(!cin.good())
      return;
    cout << endl;

    // read project name
    cout << "Please enter the project name:" << endl << ">";
    cin >> project_name;
    if(!cin.good())
      return;
    cout << endl;

    // generate directories
    if(generate_dirs())
    {
      // search for files if the dir tree existed
      find_files();

      if(!hpp_list.empty())
      {
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

      }

      if(!cpp_list.empty())
      {
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

        if((cmd != "y") && (cmd != "yes"))
          cpp_list.clear();
      }

      if(!cu_list.empty())
      {
        cout << endl << "The following CUDA files have been found:" << endl;
        set<string>::iterator it(cu_list.begin()), jt(cu_list.end());
        for(; it != jt; ++it)
          cout << "- " << *it << endl;

        string cmd;
        cout << endl << "Do you want these files to be included in the project file?" << endl;
        cout << "Type 'y' for yes or 'n' for no" << endl;
        cout << ">";
        cin >> cmd;
        cout << endl;

        if((cmd != "y") && (cmd != "yes"))
          cu_list.clear();
      }
    }

    // generate
    generate_files();
  }

  void parse_list(set<string>& list, string str)
  {
    std::size_t pos(0), pos2;
    for(;;)
    {
      pos2 = str.find(';', pos);
      if(pos2 == str.npos)
      {
        list.insert(str.substr(pos));
        return;
      }
      else
        list.insert(str.substr(pos, pos2 - pos - 1));
    }
  }

  void main(ArgList args)
  {
    // check for test
    if(args.front().compare("-test") == 0)
    {
      main_test(args);
      return;
    }

    // fetch project path
    project_path = args.next();

    // fetch project name
    if(!args.empty() && (args.front().front() != '-'))
      project_name = args.next();
    else
    {
      // try to extract project name from project path
      size_t pos1 = project_path.find_last_of('\\');
      size_t pos2 = project_path.find_last_of('/');
      size_t pos = project_path.npos;
      if((pos1 != pos) && (pos2 != pos))
        pos = (pos1 < pos ? pos1 : pos2);
      else if(pos1 != pos)
        pos = pos1;
      else if(pos2 != pos)
        pos = pos2;
      else
        throw string("Could not extract project name from '" + project_path + "' !");

      // extract project name from project path
      project_name = project_path.substr(pos+1);
    }

    // check remaining arguments
    while(!args.empty())
    {
      // fetch next arguments
      string arg(args.next());

      if(arg.substr(0,5).compare("-hpp:") == 0)
      {
        string str(arg.substr(5));
        if(str == "*")
          find_files(hpp_list, "hpp");
        else
          parse_list(hpp_list, str);
      }
      else if(arg.substr(0,5).compare("-cpp:") == 0)
      {
        string str(arg.substr(5));
        if(str == "*")
          find_files(cpp_list, "cpp");
        else
          parse_list(cpp_list, str);
      }
      else if(arg.substr(0,4).compare("-cu:") == 0)
      {
        string str(arg.substr(4));
        if(str == "*")
          find_files(cu_list, "cu");
        else
          parse_list(cu_list, str);
      }
      else // unknown argument
        throw string("Unknown argument: '"  + arg + "'");
    }

    // generate directories
    generate_dirs();

    // generate files
    generate_files();
  }

  void main_test(ArgList args)
  {
    // remove the '-test'
    args.pop_front();
    testing = true;

    // fetch test file path
    string test_path = args.next();

    // split the path
    vector<string> test_dirs(split_path(test_path));

    // assume that the test is in the kernel directory
    if(test_dirs.front().compare("kernel") != 0)
    {
      cout << "ERROR: '" << test_path << "' does not refer to a kernel test" << endl;
      return;
    }

    // build the project path
    project_path = "testing";
    for(size_t i(1); i+1 < test_dirs.size(); ++i)
      project_path += "\\" + test_dirs.at(i);

    // set project name
    project_name = remove_ext(test_dirs.back());

    // generate dirs
    generate_dirs();

    // add test path
    cpp_list.insert(root_path + "\\" + test_path);

    // generate files
    generate_files();
  }
};

// application main entrypoint
int main(int argc, char** argv)
{
  // print header
  cout << endl << "FEAT Visual Studio 2015 Project File Generator" << endl << endl;

  // check for usage/help
  if(argc > 1)
  {
    string arg(argv[1]);
    if((arg.compare("-?") == 0) || (arg.compare("/?") == 0) || (arg.compare("-help") == 0))
    {
      // print help message
      //       123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
      cout << "This tool can be used to generate FEAT application project files for" << endl;
      cout << "Microsoft (R) Visual Studio 2015." << endl;
      cout << endl;
      cout << "Important: This tool must be executed from the FEAT root directory, as it" << endl;
      cout << "           will not be able to resolve references otherwise." << endl;
      cout << endl;
      cout << "USAGE: vc14gen.exe [-test] <path> [<name>] [options]" << endl;
      cout << endl;
      cout << "If the '-test' option is given, then <path> specifies the path to a kernel" << endl;
      cout << "test source file for which a test project is to be generated. A corresponding" << endl;
      cout << "test project file will be generated in the 'testing' directory of the FEAT" << endl;
      cout << "root directory." << endl;
      cout << endl;
      cout << "If the '-test' option is not given, then this tool will generate an application" << endl;
      cout << "project file, with the following supported parameters:" << endl;
      cout << "<path>             The path to the application directory relative to the FEAT" << endl;
      cout << "                   root directory." << endl;
      cout << endl;
      cout << "<name>             The name of the project; will be used as the filename for" << endl;
      cout << "                   the project files. If not given, the name of the last" << endl;
      cout << "                   directory of the path is taken." << endl;
      cout << endl;
      cout << "Further options:" << endl;
      cout << "-hpp:<list>        Include all C++ header files of the semicolon-separated" << endl;
      cout << "                   list. Specify '-hpp:*' to include all C++ header files" << endl;
      cout << "                   in the application directory automatically." << endl;
      cout << endl;
      cout << "-cpp:<list>        Include all C++ source files in the semicolon-separated" << endl;
      cout << "                   list. Speficy '-cpp:*' to include all C++ source files" << endl;
      cout << "                   in the application directory automatically." << endl;
      cout << endl;
      cout << "-cu:<list>         Include all CUDA source files in the semicolon-separated" << endl;
      cout << "                   list. Speficy '-cu:*' to include all CUDA source files" << endl;
      cout << "                   in the application directory automatically." << endl;
      cout << endl;
      return 0;
    }
  }

  try
  {
    Generator generator;
    if(argc < 2)
      generator.main();
    else
    {
      ArgList args(argc, argv);
      generator.main(args);
    }
  }
  catch(string& s)
  {
    cout << endl << "ERROR: " << s << endl;
  }
  return 0;
}
