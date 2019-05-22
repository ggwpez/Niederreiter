TEMPLATE = app
CONFIG -= console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += debug_and_release

SOURCES += \
	src/main.cpp \
	src/bgc.cpp \
	src/printer.cpp \
	src/helper.cpp \
	src/tests.cpp \
	src/rand_helper.cpp \
	src/perm_gf2.cpp \
	src/hasher.cpp \
	src/ncs.cpp \
	src/binom.cpp \
	src/serializer.cpp \
	src/arg_parse.cpp

HEADERS += \
	include/printer.hpp \
	include/helper.hpp \
	include/tests.hpp \
	include/rand_helper.hpp \
	include/perm_gf2.hpp \
	include/hasher.hpp \
	include/bgc.hpp \
	include/ncs.hpp \
	include/binom.hpp \
	include/serializable.hpp \
	include/serializer.hpp \
	include/arg_parse.hpp

INCLUDEPATH += include/

#QMAKE_CXXFLAGS += -std=c++11
QMAKE_LIBS += -lntl -lpthread -lgmp
QMAKE_CXXFLAGS_DEBUG += -DDEBUG -DNTL_RANGE_CHECK -g3
QMAKE_CXXFLAGS_RELEASE += -Wall -march=native -msse4.2 -mpopcnt
#QMAKE_LFLAGS += -fuse-ld=gold
#-Werror
#QMAKE_CXXFLAGS -= -fPIC
GIT_OUTPUT=$$system(git rev-parse --short HEAD)
DEFINES += GIT_HASH='0x$$GIT_OUTPUT'
