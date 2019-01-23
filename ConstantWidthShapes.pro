#-------------------------------------------------
#
# Project created by QtCreator 2018-12-15T16:24:02
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ConstantWidthShapes
TEMPLATE = app
INCLUDEPATH += E:/Qt/Libs/eigen-eigen-323c052e1731

SOURCES += main.cpp\
        mainwindow.cpp \
    constantwidthgen.cpp

HEADERS  += mainwindow.h \
    constantwidthgen.h

FORMS    += mainwindow.ui
