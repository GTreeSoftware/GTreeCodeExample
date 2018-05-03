#ifndef CREATHDF5DIALOG_H
#define CREATHDF5DIALOG_H

#include <QDateTime>
#include <QEvent>
#include "ngtypes/ParamPack.h"
#include "Function/IO/imagereader.h"
#include "ui_creathdf5dialog.h"

class creathdf5dialog : public QDialog
{
	Q_OBJECT

public:
	creathdf5dialog(QWidget *parent = 0);
	~creathdf5dialog();
	
	QString saveH5;
	std::vector<std::string> tiffList;
	bool closeFlag = true;
	bool tiffFlag = false;
	double xResolution = 0;
	double yResolution = 0;
	double zResolution = 0;

	void SetParam(NGParamPack arg){ paramPack = arg; }

private slots:
    void on_okButton_clicked();
	void on_SelectButton_clicked();
	void on_pathButton_clicked();
	void on_cancelButton_clicked();
	void on_enterPressed();
	int getMaxLevelNum(int,int,int);

private:
	Ui::creathdf5dialog ui;
	QString fileList;
	QString saveDir;
	QDateTime selectTime;
	bool pathFlag = false;
	NGParamPack paramPack;
	
};

#endif // CREATHDF5DIALOG_H
