TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle

SOURCES += main.cpp \
    demreader.cpp

LIBS += \
    -losg -losgDB -losgViewer -losgGA -losgUtil -losgQt \
    -lOpenThreads -losgAnimation -losgText -losgSim -losgTerrain -lgdal
LIBS += -L$$PWD/../../../../../usr/local/lib/ -lopencv_core \
-lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs
HEADERS += \
    util.h \
    demreader.h

