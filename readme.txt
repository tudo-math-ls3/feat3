FEAST README

How to build feast:

In the top-lvl directory run:
cmake .
make


If you want to use out-of-source builds:
Delete (if it exists) the CMakeCache.txt file in $(PATH_TO_FEAST_SRC_TOP_LVL_DIR).
Switch to the target directory (e.g. nobackup/feastbuild/) and run
cmake $(PATH_TO_FEAST_SRC_TOP_LVL_DIR)
make
with $(PATH_TO_FEAST_SRC_TOP_LVL_DIR) beeing the directory of the feast sources.


Available make targets in the top lvl build directory are:
help - show all available make targets
tests - build all basic tests
check - build and run all basic tests
doc - build doxygen documentation (found in the folder doc)


Adding new tests or applications can be done via the local folder's CMakeLists.txt file.
Every folders CMakeLists.txt contains a list of tests named test_list.
All you have to do is add your files name to this list.
To add for example the file new-test.cpp to the list represented by
SET ( test_list base_header-test )
it has to be changed to
SET ( test_list base_header-test new-test)
containing the new file name (without the file ending .cpp).
The created test will also be named new-test.

To add an application myapp.cpp to the local application folder's CMakeLists.txt file, all that has to be
done is adding
ADD_EXECUTABLE(myappbin myapp.cpp)
to the CMakeLists.txt file.
Here myappbin is the executable's name and myapp is the sourcefile's name.
