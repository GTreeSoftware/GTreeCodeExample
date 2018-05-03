#include "creathdf5dialog.h"

#include <QDialog>
#include <QString>
#include <QStringList>
#include <QFileDialog>
#include <QMessageBox>
#include <QObject>
#include <QDir>

creathdf5dialog::creathdf5dialog(QWidget *parent)
	: QDialog(parent)
{
	ui.setupUi(this);

	//set initia data
	this->setFixedSize(474, 420);
	this->setWindowTitle("Create HDF5 File");

	ui.pathEdit->setToolTip("Enter to confirm!");
	ui.okButton->setEnabled(false);
	ui.xSpinBox->setValue(1.0);
	ui.xSpinBox->setSingleStep(0.1);
	ui.ySpinBox->setValue(1.0);
	ui.ySpinBox->setSingleStep(0.1);
	ui.zSpinBox->setValue(1.0);
	ui.zSpinBox->setSingleStep(0.1);

	//QObject::connect(ui.pathEdit, SIGNAL(returnPressed()), this, SLOT(on_enterPressed()));
	QObject::connect(ui.pathEdit, SIGNAL(editingFinished()), this, SLOT(on_enterPressed()));
}

creathdf5dialog::~creathdf5dialog()
{

}

void creathdf5dialog::on_SelectButton_clicked()
{
	//get path of tif image and display current file info
	ui.infoBrowser->clear();
	tiffList.empty();
	fileList.clear();
	if (!pathFlag)
		ui.pathEdit->clear();

	QStringList qtiffList = QFileDialog::getOpenFileNames(this, "Choose TIFF Sequence", NULL, tr("TIFF (*.tif)"));
	if (qtiffList.isEmpty() || qtiffList.size() < 2)
	{
		QMessageBox::warning(this, "Warning", "Please select at least two images!");
		ui.okButton->setEnabled(false);
		closeFlag = true;
		return;
	}
	else
	{
		tiffFlag = true;
		std::for_each(qtiffList.begin(), qtiffList.end(), [&](const QString &file){tiffList.push_back(file.toStdString().c_str()); });
		std::for_each(qtiffList.begin(), qtiffList.end(), [&](QString &file){fileList += (file.section('/', 0, -1)) + "\n"; });
		if (!pathFlag)
		{
			QString tem(qtiffList[0]);
			saveDir = tem.replace(tem.section('/', -1), "");
			ui.pathEdit->setText(saveDir);
		}
		selectTime = QDateTime::currentDateTime();
		ui.infoBrowser->setText(fileList);
		QString temp1 = ui.nameEdit->text();
		if (temp1.isEmpty())
			ui.nameEdit->setText(selectTime.toString("yyyy_MM_dd_hh_mm") + ".h5");

		NGImageReader reader = ImageReader::New();
		reader->SetParam(paramPack);
		reader->GetImageInfo(tiffList[0]);
		int maxLevelNum = getMaxLevelNum(paramPack->xRangeMax_, paramPack->yRangeMax_, paramPack->zRangeMax_);
		int readerBpp = (reader->GetImageType() == DATATYPE::IMAGE8 ? 8 : 16);
		ui.infoLabel->setText(QString(tr("%1 bpp. MaxLevelNum: %2").arg(readerBpp).arg(maxLevelNum)));
		ui.okButton->setEnabled(true);
	}
}

void creathdf5dialog::on_pathButton_clicked()
{
	saveH5 = QFileDialog::getSaveFileName(this, "Create HDF5", saveDir, tr("HDF5 (*.h5)"));
	if (saveH5.isEmpty())
		return;
	else
	{
		QString tem(saveH5);
		saveDir = tem.replace(tem.section('/', -1), "");
		ui.pathEdit->setText(saveDir);
		ui.nameEdit->setText(saveH5.section('/', -1));
	}
}

void creathdf5dialog::on_okButton_clicked()
{
	if (saveH5.isEmpty())
	    saveH5 = saveDir + selectTime.toString("yyyy_MM_dd_hh_mm") + ".h5";

	xResolution = ui.xSpinBox->value();
	yResolution = ui.ySpinBox->value();
	zResolution = ui.zSpinBox->value();
	closeFlag = false;//The false delegate did not cancel the creat 
	this->close();
}

void creathdf5dialog::on_cancelButton_clicked()
{
	closeFlag = true;
	this->close();
}

int creathdf5dialog::getMaxLevelNum(int xR, int yR, int zR)
{
	int levelNum = 1;
	int maxLength = xR > yR ? (xR > zR ? xR : zR) : (yR > zR ? yR : zR);
	int minLength = xR < yR ? (xR < zR ? xR : zR) : (yR < zR ? yR : zR);
	//int maxLengthCp = maxLength;
	//int minLengthCp = minLength;
	while (true){
		if (maxLength < 256 || (minLength < 50 && maxLength < 2 * 256) || minLength < 10)
			break;
		maxLength /= 2;
		minLength /= 2;
		++levelNum;
	}
	return levelNum;
}

void creathdf5dialog::on_enterPressed()
{
	QString tem = ui.pathEdit->text();
	if (tem.isEmpty() || !ui.pathEdit->hasFocus())
		return;

	tem.replace('\\', '/');
	QDir testDir(tem);
	if (!testDir.exists())
	{
		QMessageBox::warning(this, "ERROR", "The path for save does not exits!");
		return;
	}

	if (testDir.isRoot())
	{
		QMessageBox::warning(this, "EROOR", "The path is root!");
		return;
	}
	else
	{
		if (!tem.endsWith('/'))
			tem += '/';
		pathFlag = true;
		ui.pathEdit->setText(tem);
		saveDir = tem;
		selectTime = QDateTime::currentDateTime();
		ui.nameEdit->setText(selectTime.toString("yyyy_MM_dd_hh_mm") + ".h5");
		ui.pathEdit->clearFocus();
	}
}
