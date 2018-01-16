#include <QtWidgets\qapplication.h>
#include "new3d.h"
int main(int argc, char* argv[])
{
    QApplication a(argc, argv);
    NEW3D w;
    w.show();
    return a.exec();
}