
unix:!macx: LIBS += -L$$PWD/../../../../../usr/local/lib/ -lntl

INCLUDEPATH += $$PWD/../../../../../usr/local/ntl-10.3.0/include
DEPENDPATH += $$PWD/../../../../../usr/local/ntl-10.3.0/include

unix:!macx: PRE_TARGETDEPS += $$PWD/../../../../../usr/local/lib/libntl.a

unix:!macx: LIBS += -L$$PWD/../../../../../usr/local/gmp-6.1.2/.libs/ -lgmp

INCLUDEPATH += $$PWD/../../../../../usr/local/gmp-6.1.2
DEPENDPATH += $$PWD/../../../../../usr/local/gmp-6.1.2

unix:!macx: PRE_TARGETDEPS += $$PWD/../../../../../usr/local/gmp-6.1.2/.libs/libgmp.a

SOURCES += \
    main1.cxx \
    test.cpp \
    reduce.cpp \
    gen.cpp \
    checkroot.cpp

OTHER_FILES += \
    cop_mat.txt \
    reduced_mat.txt \
    param.txt \
    size_reduced_mat.txt \
    rounded_mat.txt

HEADERS += \
    stdinc.h \
    inch.h \
    test.h \
    reduce.h \
    gen.h \
    checkroot.h

CONFIG += c++11

