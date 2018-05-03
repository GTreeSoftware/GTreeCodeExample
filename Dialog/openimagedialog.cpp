/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#include <tiffio.h>
#include <QFileDialog>
#include <QValidator>
#include <QMessageBox>
#include <QDoubleSpinBox>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif
#include "openimagedialog.h"
#include "ui_openimagedialog.h"
#include "Function/IO/MOSTDReader.h"
#include "Function/IO/imagereader.h"
#include "Function/IO/HDF5Reader.h"

OpenImageDialog::OpenImageDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::OpenImageDialog)
{
    ui->setupUi(this);
    this->setFixedSize(420,350);
    bitsPerPixel_ = 0;
    x_ = 0;
    y_ = 0;
    z_ = 0;
    xResolution_ = 1.0;
    yResolution_ = 1.0;
    zResolution_ = 1.0;
    ui->OKButton->setEnabled(false);
    //connect(ui->SomaPathEdit, SIGNAL(editingFinished()), this, UpdateSomaPath());
}

OpenImageDialog::~OpenImageDialog()
{
    delete ui;
}

const QString OpenImageDialog::GetOpenImageName() const
{
    return ui->FilePathEdit->text();
}

void OpenImageDialog::on_OpenButton_clicked()
{
    readImageName_ = QFileDialog::getOpenFileName(this, tr("Open Image File"), paramPack->defaultDir,
        tr("IMAGES (*.tif *.mostd *.h5)")); 
    if (readImageName_.isEmpty()) {
        printf("not open any file.\n");
        return;
    }
    if (readImageName_.section('.', -1).compare("tif") == 0) {
        paramPack->dataMode_ = READMODE::TIFF3D;
        QDialog *tmp = new QDialog(this);
        QVBoxLayout *vl = new QVBoxLayout(tmp);
        tmp->setLayout(vl);
        QHBoxLayout* hl1 = new QHBoxLayout();//layer
        QHBoxLayout *hl2 = new QHBoxLayout();//branch
        QLabel *inputLayerLabel = new QLabel(this); inputLayerLabel->setText("input Resolution");
        QDoubleSpinBox *inputXRes = new QDoubleSpinBox(tmp); inputXRes->setValue(1.0); inputXRes->setSingleStep(0.1);
        QDoubleSpinBox *inputYRes = new QDoubleSpinBox(tmp); inputYRes->setValue(1.0); inputYRes->setSingleStep(0.1);
        QDoubleSpinBox *inputZRes = new QDoubleSpinBox(tmp); inputZRes->setValue(1.0); inputZRes->setSingleStep(0.1);
        vl->addWidget(inputLayerLabel); vl->addLayout(hl1); vl->addLayout(hl2);
        hl1->addWidget(inputXRes); hl1->addWidget(inputYRes); hl1->addWidget(inputZRes);
        QPushButton *acceptButton = new QPushButton("accept", this);
        QPushButton *cancelButton = new QPushButton("cancel", this);
        hl2->addWidget(acceptButton);
        hl2->addWidget(cancelButton);
        QObject::connect(acceptButton, SIGNAL(clicked()), tmp, SLOT(accept()));
        QObject::connect(cancelButton, SIGNAL(clicked()), tmp, SLOT(reject()));
        if (tmp->exec() == QDialog::Accepted) {
            paramPack->xRes_ = inputXRes->value();
            paramPack->yRes_ = inputYRes->value();
            paramPack->zRes_ = inputZRes->value();
        }
    }
    else if (readImageName_.section('.', -1).compare("h5") == 0){
        paramPack->dataMode_ = READMODE::HDF5;
    }
    else{
        paramPack->dataMode_ = READMODE::MOSTD;
    }
    //paramPack->defaultPath = readFileName.section('/', 0, -2);
    SetDialogFileInfo(readImageName_);
    if(!ui->FilePathEdit->text().isEmpty()){
        ui->OKButton->setEnabled(true);
    }
    //automatically detect soma file
    QString autoStackFile = readImageName_.section('.',0,-2) + QString(".sta");
#ifdef _WIN32
    if(0 == _access(autoStackFile.toStdString().c_str(), 0) ){
#else
    if (0 == access(autoStackFile.toStdString().c_str(), 0)){
#endif
        ui->StackPathEdit->setText(autoStackFile);
        readStackName_ = autoStackFile;
    }
    paramPack->defaultDir = readImageName_.section('/', 0, -2);
    //2015-6-16
}

void OpenImageDialog::on_OpenStackButton_clicked()
{
    //if (readImageName_.isEmpty() || readImageName_.section('.', -1).compare("tif") == 0){
        readSwcList_ = QFileDialog::getOpenFileNames(this, tr("Open Trace File"), paramPack->defaultDir,
            tr("Trace (*.swc *.sta *.cbd)"));
    //}
    if (readSwcList_.isEmpty()) {
        printf("not open any file.\n");
        return;
    }
    paramPack->defaultDir = readSwcList_[0].section('/', 0, -2);
    //check validation
    int staNum = 0;
    int swcNum = 0;
	int cbdNum = 0;///////////////////////////////////////////////
    //bool mixFlag = false;
    for (auto &it : readSwcList_) {
        if (it.section('.', -1).compare("sta") == 0) ++staNum;
        if (staNum > 1) {
            QMessageBox::warning(this, "trace info warning", "You can import only one sta file.");
            return;
        }

        if (it.section('.', -1).compare("swc") == 0) ++swcNum;

		if (it.section('.', -1).compare("cbd") == 0) ++cbdNum;///////////////////////////////////////
    }
	if (staNum > 0 && swcNum > 0 || staNum&&cbdNum || swcNum&&cbdNum || staNum&&swcNum&&cbdNum) {
        QMessageBox::warning(this, "trace info warning", "You can not import different file.");//not two or more
        return;
    }
    //build text
    QString editString;
    for (int i = 0; i < readSwcList_.size(); ++i) 
        editString += readSwcList_[i] + ';';
    //
    ui->StackPathEdit->setText(editString);
	
	if (!ui->FilePathEdit->text().isEmpty() && !ui->StackPathEdit->text().isEmpty()){
		ui->OKButton->setEnabled(true);
	}
}

void OpenImageDialog::SetDialogFileInfo(const QString &name)
{
    if (name.isEmpty()) {
        return;
    }
    ui->FilePathEdit->setText(name);
    ui->getFileNameLable->setText(QString(name.section("/", -1)));
    GetImageInfo((const char*)name.toLocal8Bit());//2015-6-16
}

void OpenImageDialog::GetImageInfo(const char *path)
{
    switch (paramPack->dataMode_){
    case READMODE::MOSTD:{
        NGMOSTDReader reader = MOSTDReader::New();
        reader->SetParam(paramPack);
        reader->SetInputFileName(std::string(path));
        ui->getImageSizeLabel->setText(QString(tr("%1 slices, %2 * %3, %4 bpp").arg(paramPack->zRangeMax_).arg(
            paramPack->xRangeMax_).arg(paramPack->yRangeMax_).arg(reader->GetImageType() == DATATYPE::IMAGE8 ? 8 : 16)));
        ui->VoxelSizeLabel->setText(QString(tr("Voxel Size:%1 um x %2 um x %3 um").arg(paramPack->xRes_).arg(
            paramPack->yRes_).arg(paramPack->zRes_)));
        dataType_ = reader->GetImageType();
        break;
    }
    case READMODE::TIFF3D:{
        NGImageReader reader = ImageReader::New();
        reader->SetParam(paramPack);
        reader->GetImageInfo(std::string(path));
        ui->getImageSizeLabel->setText(QString(tr("%1 slices, %2 * %3, %4 bpp").arg(paramPack->zRangeMax_).arg(
            paramPack->xRangeMax_).arg(paramPack->yRangeMax_).arg(reader->GetImageType() == DATATYPE::IMAGE8 ? 8 : 16)));
        ui->VoxelSizeLabel->setText(QString(tr("Voxel Size:%1 um x %2 um x %3 um").arg(paramPack->xRes_).arg(
            paramPack->yRes_).arg(paramPack->zRes_)));
        dataType_ = reader->GetImageType();
        break;
    }
    case READMODE::HDF5:{
        NGHDF5Reader reader = HDF5Reader::New();
        reader->SetParam(paramPack);
        reader->SetInputFileName(std::string(path));
        ui->getImageSizeLabel->setText(QString(tr("%1 slices, %2 * %3, %4 bpp").arg(paramPack->zRangeMax_).arg(
            paramPack->xRangeMax_).arg(paramPack->yRangeMax_).arg(reader->GetImageType() == DATATYPE::IMAGE8 ? 8 : 16)));
        ui->VoxelSizeLabel->setText(QString(tr("Voxel Size:%1 um x %2 um x %3 um").arg(paramPack->xRes_).arg(
            paramPack->yRes_).arg(paramPack->zRes_)));
        dataType_ = reader->GetImageType();
        break;
    }
    default:
        break;
    }
}

const QString OpenImageDialog::GetOpenStackName() const
{
    return ui->StackPathEdit->text();
}
