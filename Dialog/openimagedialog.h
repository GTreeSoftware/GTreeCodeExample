//time : 2013-11-20
//this file declare a dialog to open image and configure
//author:zhouhang,WLNO BMP

#ifndef OPENIMAGEDIALOG_H
#define OPENIMAGEDIALOG_H

#include <QDialog>
#include <QString>
#include "ngtypes/ParamPack.h"

namespace Ui {
class OpenImageDialog;
}

class OpenImageDialog : public QDialog
{
    Q_OBJECT

public:
    explicit OpenImageDialog(QWidget *parent = 0);
    ~OpenImageDialog();
    const QString GetOpenImageName() const;
    const QString GetOpenStackName() const;
    void SetParam(NGParamPack arg){paramPack = arg;}
    DATATYPE GetDataType()const{ return dataType_; }

private slots:
    void on_OpenButton_clicked();

    void on_OpenStackButton_clicked();

private:
    Ui::OpenImageDialog *ui;
    QString readImageName_;
    QString readStackName_;
    QStringList readSwcList_;
    int bitsPerPixel_ ;
    int x_ ;
    int y_ ;
    int z_ ;
    double xResolution_ ;
    double yResolution_ ;
    double zResolution_ ;
    DATATYPE dataType_;
    NGParamPack paramPack;

    //function
    void SetDialogFileInfo(const QString&);//real time display current file info
    void GetImageInfo(const char*);
};

#endif // OPENIMAGEDIALOG_H
