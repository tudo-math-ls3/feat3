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
#define KERNEL_PROJECT_PATH "./kernel/vc10/kernel.vcxproj"

// path to template directory
#define TMPL_PATH "./tools/vc10_gen/templates/"

// template markers
#define TMPL_KRN_GUID "[KRN_GUID]"
#define TMPL_KRN_PATH "[KRN_PATH]"
#define TMPL_PRJ_GUID "[PRJ_GUID]"
#define TMPL_PRJ_NAME "[PRJ_NAME]"
#define TMPL_REL_ROOT "[REL_ROOT]"

// reads the kernel project file GUID
bool read_kernel_guid(string& guid)
{
  cout << "Checking Kernel Project File '" << KERNEL_PROJECT_PATH << "'..." << endl;
  // try to open the stream
  ifstream ifs(KERNEL_PROJECT_PATH, std::ios_base::in);
  if(!ifs.is_open() || !ifs.good())
    return false;

  string line;
  line.reserve(512);
  size_t pos0, pos1;
  while(!ifs.eof() && ifs.good())
  {
    getline(ifs, line);
    if((pos0 = line.find("<ProjectGuid>")) == line.npos)
      continue;
    if((pos1 = line.find("</ProjectGuid>")) == line.npos)
      return false;  // missing marker
    guid = line.substr(pos0 + 13, pos1 - pos0 - 13);
    return true;
  }

  return false;
}

int gen_random_int()
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
string gen_random_guid()
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

void replace_marker(string& str, string marker, string rep_str)
{
  size_t mlen(marker.size());
  size_t rlen(rep_str.size());
  size_t pos;
  while((pos = str.find(marker)) != str.npos)
  {
    str.replace(pos, mlen, rep_str, 0, rlen);
  }
}

bool gen_dirs(string project_root, bool& existed, int& depth)
{
  vector<string> dirs;
  depth = 0u;

  // separate paths
  size_t n0(0);
  while(n0 != project_root.npos)
  {
    // find separators
    size_t n1 = project_root.find_first_of('\\', n0);
    size_t n2 = project_root.find_first_of('/', n0);
    size_t n3(0);
    if((n1 == project_root.npos) && (n2 == project_root.npos))
    {
      dirs.push_back(project_root.substr(n0));
      break;
    }
    else if(n1 == project_root.npos)
      n3 = n2;
    else if(n2 == project_root.npos)
      n3 = n1;
    else
      n3 = min(n1, n2);

    // separate strings
    dirs.push_back(project_root.substr(n0, n3));
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

  // also create the vc10 directory
  depth = int(dirs.size());
  dirs.push_back("vc10");
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
    else if(i + 1u < dirs.size())
    {
      // directory created
      existed = false;
    }
  }

  // okay
  return true;
}

bool parse_lists(string project_root, set<string>& hpp_list, set<string>& cpp_list)
{
  WIN32_FIND_DATA find_data;
  HANDLE find_handle;

  // find hpp files
  memset(&find_data, 0, sizeof(WIN32_FIND_DATA));
  find_handle = FindFirstFile(string("./" + project_root + "/*.hpp").c_str(), &find_data);
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

  // find cpp files
  memset(&find_data, 0, sizeof(WIN32_FIND_DATA));
  find_handle = FindFirstFile(string("./" + project_root + "/*.cpp").c_str(), &find_data);
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

  if(hpp_list.empty() && cpp_list.empty())
    return false;

  if(!hpp_list.empty())
  {
    cout << endl << "The following header files have been found:" << endl;
    set<string>::iterator it(hpp_list.begin()), jt(hpp_list.end());
    for(; it != jt; ++it)
      cout << "- " << *it << endl;
  }

  if(!cpp_list.empty())
  {
    cout << endl << "The following source files have been found:" << endl;
    set<string>::iterator it(cpp_list.begin()), jt(cpp_list.end());
    for(; it != jt; ++it)
      cout << "- " << *it << endl;
  }

  string cmd;
  cout << endl << "Do you want these files to be included in the project file?" << endl;
  cout << "Type 'y' for yes or 'n' for no" << endl;
  cout << ">";
  cin >> cmd;
  cout << endl;

  return (cmd == "y") || (cmd == "yes");
}

bool write_solution(string sln_path, string krn_guid, string prj_guid, string prj_name, string krn_path)
{
  cout << "Writing solution file        '" << sln_path << "'..." << endl;
  ifstream ifs(TMPL_PATH "solution.xml", ios_base::in);
  if(!ifs.is_open() || !ifs.good())
  {
    cout << "ERROR: Failed to open template file '" TMPL_PATH "solution.xml'" << endl;
    return false;
  }
  ofstream ofs(sln_path, ios_base::out|ios_base::trunc);
  if(!ofs.is_open() || !ofs.good())
  {
    cout << "ERROR: Failed to create '" + sln_path + "'" << endl;
    return false;
  }

  // write utf-8 bom
  ofs << char(-17) << char(-69) << char(-65);

  string line;
  line.reserve(512);

  while(!ifs.eof())
  {
    // read line
    getline(ifs, line);

    // replace markers
    replace_marker(line, TMPL_KRN_GUID, krn_guid);
    replace_marker(line, TMPL_PRJ_GUID, prj_guid);
    replace_marker(line, TMPL_PRJ_NAME, prj_name);
    replace_marker(line, TMPL_KRN_PATH, krn_path);

    // write line
    ofs << line << endl;
  }

  ofs.close();
  ifs.close();
  return true;
}

bool write_project(string prj_path, string krn_guid, string prj_guid, string prj_name, string krn_path,
  string rel_root, const set<string>& hpp_set, const set<string>& cpp_set)
{
  cout << "Writing project file         '" << prj_path << "'..." << endl;
  ofstream ofs(prj_path, ios_base::out|ios_base::trunc);
  if(!ofs.is_open() || !ofs.good())
  {
    cout << "ERROR: Failed to create '" + prj_path + "'" << endl;
    return false;
  }

  // write utf-8 bom
  ofs << char(-17) << char(-69) << char(-65);

  string line;
  line.reserve(512);

  // write project head
  ifstream ifs;
  ifs.open(TMPL_PATH "project_head.xml", ios_base::in);
  if(!ifs.is_open() || !ifs.good())
  {
    cout << "ERROR: Failed to open template file '" TMPL_PATH "project_head.xml'" << endl;
    return false;
  }
  while(!ifs.eof())
  {
    // read line
    getline(ifs, line);

    // replace markers
    replace_marker(line, TMPL_KRN_GUID, krn_guid);
    replace_marker(line, TMPL_PRJ_GUID, prj_guid);
    replace_marker(line, TMPL_PRJ_NAME, prj_name);
    replace_marker(line, TMPL_KRN_PATH, krn_path);
    replace_marker(line, TMPL_REL_ROOT, rel_root);

    // write line
    ofs << line << endl;
  }
  ifs.close();

  // write header list
  if(!hpp_set.empty())
  {
    ofs << "  <ItemGroup>" << endl;
    set<string>::const_iterator it(hpp_set.begin()), jt(hpp_set.end());
    for(; it != jt; ++it)
      ofs << "    <ClInclude Include=\"..\\" << *it << "\" />" << endl;
    ofs << "  </ItemGroup>" << endl;
  }

  // write source list
  if(!cpp_set.empty())
  {
    ofs << "  <ItemGroup>" << endl;
    set<string>::const_iterator it(cpp_set.begin()), jt(cpp_set.end());
    for(; it != jt; ++it)
      ofs << "    <ClCompile Include=\"..\\" << *it << "\" />" << endl;
    ofs << "  </ItemGroup>" << endl;
  }

  // write project tail
  ifs.open(TMPL_PATH "project_tail.xml", ios_base::in);
  if(!ifs.is_open() || !ifs.good())
  {
    cout << "ERROR: Failed to open template file '" TMPL_PATH "project_tail.xml'" << endl;
    return false;
  }
  while(!ifs.eof())
  {
    // read line
    getline(ifs, line);

    // write line
    ofs << line << endl;
  }
  ifs.close();

  ofs.close();
  return true;
}

bool write_filters(string flt_path, const set<string>& hpp_set, const set<string>& cpp_set)
{
  cout << "Writing project filters file '" << flt_path << "'..." << endl;
  ofstream ofs(flt_path, ios_base::out|ios_base::trunc);
  if(!ofs.is_open() || !ofs.good())
  {
    cout << "ERROR: Failed to create '" + flt_path + "'" << endl;
    return false;
  }

  // write utf-8 bom
  ofs << char(-17) << char(-69) << char(-65);

  ofs << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
  ofs << "<Project ToolsVersion=\"4.0\" xmlns=\"http://schemas.microsoft.com/developer/msbuild/2003\">" << endl;

  // generate filter guids
  string hpp_guid(gen_random_guid());
  string cpp_guid(gen_random_guid());

  // write filter headers
  ofs << "  <ItemGroup>" << endl;
  ofs << "    <Filter Include=\"Header Files\">" << endl;
  ofs << "      <UniqueIdentifier>" << hpp_guid << "</UniqueIdentifier>" << endl;
  ofs << "    </Filter>" << endl;
  ofs << "    <Filter Include=\"Source Files\">" << endl;
  ofs << "      <UniqueIdentifier>" << cpp_guid << "</UniqueIdentifier>" << endl;
  ofs << "    </Filter>" << endl;
  ofs << "  </ItemGroup>" << endl;

  // write header list
  if(!hpp_set.empty())
  {
    ofs << "  <ItemGroup>" << endl;
    set<string>::const_iterator it(hpp_set.begin()), jt(hpp_set.end());
    for(; it != jt; ++it)
    {
      ofs << "    <ClInclude Include=\"..\\" << *it << "\">" << endl;
      ofs << "      <Filter>Header Files</Filter>" << endl;
      ofs << "    </ClInclude>" << endl;
    }
    ofs << "  </ItemGroup>" << endl;
  }

  // write source list
  if(!cpp_set.empty())
  {
    ofs << "  <ItemGroup>" << endl;
    set<string>::const_iterator it(cpp_set.begin()), jt(cpp_set.end());
    for(; it != jt; ++it)
    {
      ofs << "    <ClCompile Include=\"..\\" << *it << "\">" << endl;
      ofs << "      <Filter>Source Files</Filter>" << endl;
      ofs << "    </ClCompile>" << endl;
    }
    ofs << "  </ItemGroup>" << endl;
  }

  ofs << "</Project>" << endl;
  ofs.close();

  return true;
}

bool write_user(string usr_path, string rel_root)
{
  cout << "Writing project user file    '" << usr_path << "'..." << endl;
  ofstream ofs(usr_path, ios_base::out|ios_base::trunc);
  if(!ofs.is_open() || !ofs.good())
  {
    cout << "ERROR: Failed to create '" + usr_path + "'" << endl;
    return false;
  }

  // write utf-8 bom
  ofs << char(-17) << char(-69) << char(-65);

  string line;
  line.reserve(512);

  // write project head
  ifstream ifs;
  ifs.open(TMPL_PATH "project_user.xml", ios_base::in);
  if(!ifs.is_open() || !ifs.good())
  {
    cout << "ERROR: Failed to open template file '" TMPL_PATH "project_user.xml'" << endl;
    return false;
  }
  while(!ifs.eof())
  {
    // read line
    getline(ifs, line);

    // replace markers
    replace_marker(line, TMPL_REL_ROOT, rel_root);

    // write line
    ofs << line << endl;
  }
  ifs.close();

  // okay
  return true;
}

void gen_app_project(string kernel_guid)
{
  // read project root
  string project_root;

  // read project path
  cout << endl << "Please enter the path to the project directory:" << endl;
  cout << ">";
  cin >> project_root;
  cout << endl;

  int path_depth(0);
  bool path_existed(false);
  if(!gen_dirs(project_root, path_existed, path_depth))
    return;

  // read project name
  string project_name;
  cout << "Please enter the project name:" << endl << ">";
  cin >> project_name;
  cout << endl;

  // generate project and solution guids
  string project_guid(gen_random_guid());
  cout << "Project GUID: " << project_guid << endl;

  // build paths
  string sln_path = "./" + project_root /*+ "/" + project_name*/ + "/vc10/" + project_name + ".sln";
  string prj_path = "./" + project_root /*+ "/" + project_name*/ + "/vc10/" + project_name + ".vcxproj";
  string flt_path = "./" + project_root /*+ "/" + project_name*/ + "/vc10/" + project_name + ".vcxproj.filters";
  string usr_path = "./" + project_root /*+ "/" + project_name*/ + "/vc10/" + project_name + ".vcxproj.user";

  // build relative kernel project path
  string krn_path;
  for(int i(0); i < path_depth; ++i)
    krn_path += "..\\";
  string rel_root("$(ProjectDir)" + krn_path + "..");
  krn_path += "..\\kernel\\vc10\\kernel.vcxproj";

  // build header and source lists
  set<string> hpp_set, cpp_set;
  if(!parse_lists(project_root, hpp_set, cpp_set))
  {
    hpp_set.clear();
    cpp_set.clear();
  }

  // write solution file
  if(!write_solution(sln_path, kernel_guid, project_guid, project_name, krn_path))
    return;

  // write project file
  if(!write_project(prj_path, kernel_guid, project_guid, project_name, krn_path, rel_root, hpp_set, cpp_set))
    return;

  // write filters
  if(!write_filters(flt_path, hpp_set, cpp_set))
    return;

  // write user file
  if(!write_user(usr_path, rel_root))
    return;

  cout << endl << "Project files for '" << project_name << "' have been written successfully" << endl;
}

// application main entrypoint
int main(int /*argc*/, char** /*argv*/)
{
  // print header
  cout << endl << "FEAST Visual Studio 2010 Project File Generator" << endl << endl;

  // print warning
  cout << "WARNING: This tool is still under construction, so use with care!" << endl << endl;

  // read kernel project guid
  string kernel_guid;
  if(!read_kernel_guid(kernel_guid))
  {
    cout << "ERROR: you must execute 'vc10_gen' from the FEAST root directory" << endl << endl;
    return 1;
  }
  cout << "Kernel GUID: " << kernel_guid << endl;

  // generate application project
  gen_app_project(kernel_guid);

  return 0;
}
