// ************************************************************************************************
// This file is used to generate the file 'vc11-build-modes.xml' which is used by
// the Visual Studio 2012 build system.
//
// WARNING:
// DO NOT MODIFY OR COMPILE THIS FILE UNLESS YOU KNOW *EXACTLY* WHAT YOU ARE DOING !
//
// Note:
// Ignoring the warning above might screw up the VS2012 build system.
//
// \author Peter Zajac
// ************************************************************************************************

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <fstream>

using namespace std;

struct tag
{
  string tag_on;
  string tag_off;
  string xml_on;
  string xml_off;

  tag() {}

  tag(string _tag_on, string _xml_tag, string _tag_off = "")
    : tag_on(_tag_on), tag_off(_tag_off)
  {
    xml_on  = "<" + _xml_tag + ">true</"  + _xml_tag + ">";
    xml_off = "<" + _xml_tag + ">false</" + _xml_tag + ">";
  }

  tag(string _tag_on, string _tag_off, string _xml_on, string _xml_off)
    : tag_on(_tag_on), tag_off(_tag_off), xml_on(_xml_on), xml_off(_xml_off)
  {
  }

  string t(bool b)
  {
    return !b ? tag_on : tag_off;
  }

  string x(bool b)
  {
    return !b ? xml_on : xml_off;
  }
};

int main(int argc, char* argv[])
{
  // create platform set
  vector<string> platforms;
  platforms.push_back("Win32");
  platforms.push_back("x64");
  //platforms.push_back("ARM"); // also offered by VS2012, but currently not used by FEAST

  // create tag set
  vector<tag> tags;
  tags.push_back(tag("dbg","DebugMode", "opt"));
  tags.push_back(tag("cuda", "BackendCUDA"));
  tags.push_back(tag("mkl", "BackendMKL"));
  tags.push_back(tag("mpi", "", "<SerialMode>false</SerialMode>", "<SerialMode>true</SerialMode>"));
  tags.push_back(tag("omp", "EnableOMP"));

  // create the build-mode set
  set<string> build_modes;
  size_t m = tags.size();
  size_t n = size_t(1) << m;
  for(size_t i(0); i < n; ++i)
  {
    string mode = tags[0].t((i & 0x1) != 0);
    for(size_t j(1); j < m; ++j)
    {
      string t(tags[j].t(((i >> j) & 0x1) != 0));
      if(!t.empty())
        mode += "-" + t;
    }
    build_modes.insert(mode);
  }

  std::ofstream ofs("vc11-build-modes.xml", std::ios_base::out|std::ios_base::trunc);
  if(!ofs.is_open() || !ofs.good())
  {
    cout << "ERROR: failed to open 'vc11-build-modes.xml' for writing" << endl;
    return 1;
  }

  ofs << "<Project ToolsVersion=\"4.0\" xmlns=\"http://schemas.microsoft.com/developer/msbuild/2003\">" << endl;

  // write project configurations
  ofs << "  <!-- ********************************************************************* -->" << endl;
  ofs << "  <!-- Project Configurations -->" << endl;
  ofs << "  <!-- ********************************************************************* -->" << endl;
  ofs << "  <ItemGroup Label=\"ProjectConfigurations\">" << endl;
  for(set<string>::iterator it(build_modes.begin()); it != build_modes.end(); ++it)
  {
    for(vector<string>::iterator jt(platforms.begin()); jt != platforms.end(); ++jt)
    {
      ofs << "    <ProjectConfiguration Include=\"" << *it << "|" << *jt << "\">" << endl;
      ofs << "      <Configuration>" << *it << "</Configuration>" << endl;
      ofs << "      <Platform>" << *jt << "</Platform>" << endl;
      ofs << "    </ProjectConfiguration>" << endl;
    }
  }
  ofs << "  </ItemGroup>" << endl;

  // write build-mode mapping
  ofs << "  <!-- ********************************************************************* -->" << endl;
  ofs << "  <!-- Configuration to Build-Mode mapping -->" << endl;
  ofs << "  <!-- ********************************************************************* -->" << endl;
  for(set<string>::iterator it(build_modes.begin()); it != build_modes.end(); ++it)
  {
    ofs << "  <PropertyGroup Label=\"BuildMode\" Condition=\"'$(Configuration)'=='" << *it << "'\">" << endl;
    for(vector<tag>::iterator jt(tags.begin()); jt != tags.end(); ++jt)
    {
      ofs << "    " << jt->x(it->find(jt->t(false)) == it->npos) << endl;
    }
    ofs << "  </PropertyGroup>" << endl;
  }

  ofs << "</Project>";

  ofs.close();
  return 0;
}
