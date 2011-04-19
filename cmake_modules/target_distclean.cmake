# add custom target distclean
# cleans and removes cmake generated files etc.
# Jan Woetzel 04/2003
#

IF (UNIX)
    ADD_CUSTOM_TARGET (distclean @echo cleaning for source distribution)
    ADD_DEPENDENCIES( "distclean" clean )
    SET(DISTCLEANED
        cmake.depends
        cmake.check_depends
        CMakeCache.txt
        feast_config.hpp
        cmake.check_cache
        Makefile
        core core.*
        gmon.out
        CMakeFiles
        *~
        )

    ADD_CUSTOM_COMMAND(
        DEPENDS clean
        COMMENT "distribution clean"
        COMMAND rm
        ARGS    -Rf CMakeTmp ${DISTCLEANED}
        TARGET  distclean
        )
ENDIF(UNIX)
