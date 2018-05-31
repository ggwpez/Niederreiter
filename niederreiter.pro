TEMPLATE = app
CONFIG -= console
CONFIG -= app_bundle
CONFIG -= qt c++11

SOURCES += \
	src/main.cpp \
	src/bgc.cpp \
	src/printer.cpp \
	src/helper.cpp

HEADERS += \
	include/bgc.hpp \
	include/printer.hpp \
	include/helper.hpp

QMAKE_CXXFLAGS += -std=c++17
QMAKE_LIBS += -lntl -lgmp -lm -lpthread
QMAKE_CXXFLAGS_DEBUG += -DDEBUG
INCLUDEPATH += include/
