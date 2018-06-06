TEMPLATE = app
CONFIG -= console
CONFIG -= app_bundle
CONFIG -= qt c++11

SOURCES += \
	src/main.cpp \
	src/bgc.cpp \
	src/printer.cpp \
	src/helper.cpp \
	src/tests.cpp \
    src/rand.cpp

HEADERS += \
	include/bgc.hpp \
	include/printer.hpp \
	include/helper.hpp \
	include/tests.hpp \
    include/rand.hpp

QMAKE_CXXFLAGS += -std=c++17
QMAKE_LIBS += -lntl -lgmp -lm -lpthread
QMAKE_CXXFLAGS_DEBUG += -DDEBUG
INCLUDEPATH += include/
