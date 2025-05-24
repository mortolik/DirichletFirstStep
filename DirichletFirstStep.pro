QT       += core gui
QT += datavisualization

RC_ICONS = favicon.ico
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    Dirichlet2step.cpp \
    Dirichle3StepModel.cpp \
    DirichletDisplayWidget.cpp \
    DirichletSolverModel.cpp \
    DirichletSolverModel2.cpp \
    DirichletWidget.cpp \
    main.cpp \
    mainwindow.cpp

HEADERS += \
    Dirichlet2step.hpp \
    Dirichle3StepModel.hpp \
    DirichletDisplayWidget.hpp \
    DirichletSolverModel.hpp \
    DirichletSolverModel2.hpp \
    DirichletWidget.hpp \
    mainwindow.h

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
