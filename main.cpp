#include "sparsetracer.h"
#include <QtWidgets/QApplication>
#include <string>
#include <fstream>
#include <vtkAutoInit.h>
/*
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    SparseTracer w;
    w.show();

    return a.exec();
}*/

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	QSplashScreen *splash = new QSplashScreen;
	splash->setPixmap(QPixmap(":/new/tracer/Resource/Splash.png"));
	splash->show();
	splash->showMessage("Starting...", Qt::AlignRight, Qt::blue);

	a.processEvents();

	SparseTracer w;
	w.show();

	//splash->finish(w);
	delete splash;

	return a.exec();
}