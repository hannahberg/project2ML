TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib

SOURCES += main.cpp
SOURCES += solver.cpp
HEADERS += solver.h
SOURCES += bruteforce.cpp
HEADERS += bruteforce.h
SOURCES += impsamp.cpp
HEADERS += impsamp.h
SOURCES += interact.cpp
HEADERS += interact.h
#HEADERS += catch.hpp

LIBS += -larmadillo -llapack -lblas
