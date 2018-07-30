<<<<<<< HEAD
/********************************************************************************
** Form generated from reading UI file 'openimagedialog.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OPENIMAGEDIALOG_H
#define UI_OPENIMAGEDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_OpenImageDialog
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QPushButton *OKButton;
    QPushButton *CancelButton;
    QGroupBox *FileInfoGroupBox;
    QGridLayout *gridLayout;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_2;
    QLabel *fileNameLabel;
    QLabel *getFileNameLable;
    QSpacerItem *horizontalSpacer;
    QHBoxLayout *horizontalLayout_3;
    QLabel *imageSizeLabel;
    QLabel *getImageSizeLabel;
    QSpacerItem *horizontalSpacer_2;
    QGroupBox *ResolutionRoupBox;
    QWidget *layoutWidget1;
    QVBoxLayout *verticalLayout_2;
    QLabel *ResolutionStyleLabel;
    QLabel *VoxelSizeLabel;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_5;
    QLabel *FilePathDisplayLabel;
    QLineEdit *FilePathEdit;
    QPushButton *OpenButton;
    QWidget *layoutWidget_2;
    QHBoxLayout *horizontalLayout_6;
    QLabel *SomaPathDisplayLabel;
    QLineEdit *StackPathEdit;
    QPushButton *OpenStackButton;

    void setupUi(QDialog *OpenImageDialog)
    {
        if (OpenImageDialog->objectName().isEmpty())
            OpenImageDialog->setObjectName(QStringLiteral("OpenImageDialog"));
        OpenImageDialog->resize(417, 347);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(OpenImageDialog->sizePolicy().hasHeightForWidth());
        OpenImageDialog->setSizePolicy(sizePolicy);
        horizontalLayoutWidget = new QWidget(OpenImageDialog);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(220, 310, 181, 25));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetMinimumSize);
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        OKButton = new QPushButton(horizontalLayoutWidget);
        OKButton->setObjectName(QStringLiteral("OKButton"));

        horizontalLayout->addWidget(OKButton);

        CancelButton = new QPushButton(horizontalLayoutWidget);
        CancelButton->setObjectName(QStringLiteral("CancelButton"));

        horizontalLayout->addWidget(CancelButton);

        FileInfoGroupBox = new QGroupBox(OpenImageDialog);
        FileInfoGroupBox->setObjectName(QStringLiteral("FileInfoGroupBox"));
        FileInfoGroupBox->setGeometry(QRect(10, 200, 391, 101));
        QFont font;
        font.setFamily(QStringLiteral("Arial Black"));
        FileInfoGroupBox->setFont(font);
        gridLayout = new QGridLayout(FileInfoGroupBox);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        fileNameLabel = new QLabel(FileInfoGroupBox);
        fileNameLabel->setObjectName(QStringLiteral("fileNameLabel"));
        QFont font1;
        font1.setFamily(QStringLiteral("Sans Serif"));
        fileNameLabel->setFont(font1);

        horizontalLayout_2->addWidget(fileNameLabel);

        getFileNameLable = new QLabel(FileInfoGroupBox);
        getFileNameLable->setObjectName(QStringLiteral("getFileNameLable"));
        getFileNameLable->setFont(font1);

        horizontalLayout_2->addWidget(getFileNameLable);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        imageSizeLabel = new QLabel(FileInfoGroupBox);
        imageSizeLabel->setObjectName(QStringLiteral("imageSizeLabel"));
        imageSizeLabel->setFont(font1);

        horizontalLayout_3->addWidget(imageSizeLabel);

        getImageSizeLabel = new QLabel(FileInfoGroupBox);
        getImageSizeLabel->setObjectName(QStringLiteral("getImageSizeLabel"));
        getImageSizeLabel->setFont(font1);

        horizontalLayout_3->addWidget(getImageSizeLabel);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_2);


        verticalLayout->addLayout(horizontalLayout_3);


        gridLayout->addLayout(verticalLayout, 0, 0, 1, 1);

        ResolutionRoupBox = new QGroupBox(OpenImageDialog);
        ResolutionRoupBox->setObjectName(QStringLiteral("ResolutionRoupBox"));
        ResolutionRoupBox->setGeometry(QRect(10, 90, 391, 101));
        ResolutionRoupBox->setFont(font);
        layoutWidget1 = new QWidget(ResolutionRoupBox);
        layoutWidget1->setObjectName(QStringLiteral("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(13, 31, 371, 61));
        verticalLayout_2 = new QVBoxLayout(layoutWidget1);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        ResolutionStyleLabel = new QLabel(layoutWidget1);
        ResolutionStyleLabel->setObjectName(QStringLiteral("ResolutionStyleLabel"));
        ResolutionStyleLabel->setFont(font1);

        verticalLayout_2->addWidget(ResolutionStyleLabel);

        VoxelSizeLabel = new QLabel(layoutWidget1);
        VoxelSizeLabel->setObjectName(QStringLiteral("VoxelSizeLabel"));
        VoxelSizeLabel->setSizeIncrement(QSize(3, 0));
        VoxelSizeLabel->setFont(font1);

        verticalLayout_2->addWidget(VoxelSizeLabel);

        layoutWidget = new QWidget(OpenImageDialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(10, 20, 391, 25));
        horizontalLayout_5 = new QHBoxLayout(layoutWidget);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(0, 0, 0, 0);
        FilePathDisplayLabel = new QLabel(layoutWidget);
        FilePathDisplayLabel->setObjectName(QStringLiteral("FilePathDisplayLabel"));

        horizontalLayout_5->addWidget(FilePathDisplayLabel);

        FilePathEdit = new QLineEdit(layoutWidget);
        FilePathEdit->setObjectName(QStringLiteral("FilePathEdit"));

        horizontalLayout_5->addWidget(FilePathEdit);

        OpenButton = new QPushButton(layoutWidget);
        OpenButton->setObjectName(QStringLiteral("OpenButton"));

        horizontalLayout_5->addWidget(OpenButton);

        layoutWidget_2 = new QWidget(OpenImageDialog);
        layoutWidget_2->setObjectName(QStringLiteral("layoutWidget_2"));
        layoutWidget_2->setGeometry(QRect(10, 50, 391, 25));
        horizontalLayout_6 = new QHBoxLayout(layoutWidget_2);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        horizontalLayout_6->setContentsMargins(0, 0, 0, 0);
        SomaPathDisplayLabel = new QLabel(layoutWidget_2);
        SomaPathDisplayLabel->setObjectName(QStringLiteral("SomaPathDisplayLabel"));

        horizontalLayout_6->addWidget(SomaPathDisplayLabel);

        StackPathEdit = new QLineEdit(layoutWidget_2);
        StackPathEdit->setObjectName(QStringLiteral("StackPathEdit"));

        horizontalLayout_6->addWidget(StackPathEdit);

        OpenStackButton = new QPushButton(layoutWidget_2);
        OpenStackButton->setObjectName(QStringLiteral("OpenStackButton"));

        horizontalLayout_6->addWidget(OpenStackButton);


        retranslateUi(OpenImageDialog);
        QObject::connect(OKButton, SIGNAL(clicked()), OpenImageDialog, SLOT(accept()));
        QObject::connect(CancelButton, SIGNAL(clicked()), OpenImageDialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(OpenImageDialog);
    } // setupUi

    void retranslateUi(QDialog *OpenImageDialog)
    {
        OpenImageDialog->setWindowTitle(QApplication::translate("OpenImageDialog", "OpenImage", 0));
        OKButton->setText(QApplication::translate("OpenImageDialog", "OK", 0));
        CancelButton->setText(QApplication::translate("OpenImageDialog", "Cancel", 0));
        FileInfoGroupBox->setTitle(QApplication::translate("OpenImageDialog", "FileInfo", 0));
        fileNameLabel->setText(QApplication::translate("OpenImageDialog", "File Name:", 0));
        getFileNameLable->setText(QApplication::translate("OpenImageDialog", "unknown File", 0));
        imageSizeLabel->setText(QApplication::translate("OpenImageDialog", "Image Size:", 0));
        getImageSizeLabel->setText(QApplication::translate("OpenImageDialog", "0 Slices, 0 * 0 , 0 bpp", 0));
        ResolutionRoupBox->setTitle(QApplication::translate("OpenImageDialog", "Resolution", 0));
        ResolutionStyleLabel->setText(QApplication::translate("OpenImageDialog", "Define : Voxel size (um)", 0));
        VoxelSizeLabel->setText(QApplication::translate("OpenImageDialog", "Voxel Size:", 0));
        FilePathDisplayLabel->setText(QApplication::translate("OpenImageDialog", "Image Path:", 0));
        OpenButton->setText(QApplication::translate("OpenImageDialog", "Open", 0));
        SomaPathDisplayLabel->setText(QApplication::translate("OpenImageDialog", "Stack Path:", 0));
        OpenStackButton->setText(QApplication::translate("OpenImageDialog", "Open", 0));
    } // retranslateUi

};

namespace Ui {
    class OpenImageDialog: public Ui_OpenImageDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OPENIMAGEDIALOG_H
=======
/********************************************************************************
** Form generated from reading UI file 'openimagedialog.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OPENIMAGEDIALOG_H
#define UI_OPENIMAGEDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_OpenImageDialog
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QPushButton *OKButton;
    QPushButton *CancelButton;
    QGroupBox *FileInfoGroupBox;
    QGridLayout *gridLayout;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_2;
    QLabel *fileNameLabel;
    QLabel *getFileNameLable;
    QSpacerItem *horizontalSpacer;
    QHBoxLayout *horizontalLayout_3;
    QLabel *imageSizeLabel;
    QLabel *getImageSizeLabel;
    QSpacerItem *horizontalSpacer_2;
    QGroupBox *ResolutionRoupBox;
    QWidget *layoutWidget1;
    QVBoxLayout *verticalLayout_2;
    QLabel *ResolutionStyleLabel;
    QLabel *VoxelSizeLabel;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_5;
    QLabel *FilePathDisplayLabel;
    QLineEdit *FilePathEdit;
    QPushButton *OpenButton;
    QWidget *layoutWidget_2;
    QHBoxLayout *horizontalLayout_6;
    QLabel *SomaPathDisplayLabel;
    QLineEdit *StackPathEdit;
    QPushButton *OpenStackButton;

    void setupUi(QDialog *OpenImageDialog)
    {
        if (OpenImageDialog->objectName().isEmpty())
            OpenImageDialog->setObjectName(QStringLiteral("OpenImageDialog"));
        OpenImageDialog->resize(417, 347);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(OpenImageDialog->sizePolicy().hasHeightForWidth());
        OpenImageDialog->setSizePolicy(sizePolicy);
        horizontalLayoutWidget = new QWidget(OpenImageDialog);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(220, 310, 181, 25));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetMinimumSize);
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        OKButton = new QPushButton(horizontalLayoutWidget);
        OKButton->setObjectName(QStringLiteral("OKButton"));

        horizontalLayout->addWidget(OKButton);

        CancelButton = new QPushButton(horizontalLayoutWidget);
        CancelButton->setObjectName(QStringLiteral("CancelButton"));

        horizontalLayout->addWidget(CancelButton);

        FileInfoGroupBox = new QGroupBox(OpenImageDialog);
        FileInfoGroupBox->setObjectName(QStringLiteral("FileInfoGroupBox"));
        FileInfoGroupBox->setGeometry(QRect(10, 200, 391, 101));
        QFont font;
        font.setFamily(QStringLiteral("Arial Black"));
        FileInfoGroupBox->setFont(font);
        gridLayout = new QGridLayout(FileInfoGroupBox);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        fileNameLabel = new QLabel(FileInfoGroupBox);
        fileNameLabel->setObjectName(QStringLiteral("fileNameLabel"));
        QFont font1;
        font1.setFamily(QStringLiteral("Sans Serif"));
        fileNameLabel->setFont(font1);

        horizontalLayout_2->addWidget(fileNameLabel);

        getFileNameLable = new QLabel(FileInfoGroupBox);
        getFileNameLable->setObjectName(QStringLiteral("getFileNameLable"));
        getFileNameLable->setFont(font1);

        horizontalLayout_2->addWidget(getFileNameLable);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        imageSizeLabel = new QLabel(FileInfoGroupBox);
        imageSizeLabel->setObjectName(QStringLiteral("imageSizeLabel"));
        imageSizeLabel->setFont(font1);

        horizontalLayout_3->addWidget(imageSizeLabel);

        getImageSizeLabel = new QLabel(FileInfoGroupBox);
        getImageSizeLabel->setObjectName(QStringLiteral("getImageSizeLabel"));
        getImageSizeLabel->setFont(font1);

        horizontalLayout_3->addWidget(getImageSizeLabel);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_2);


        verticalLayout->addLayout(horizontalLayout_3);


        gridLayout->addLayout(verticalLayout, 0, 0, 1, 1);

        ResolutionRoupBox = new QGroupBox(OpenImageDialog);
        ResolutionRoupBox->setObjectName(QStringLiteral("ResolutionRoupBox"));
        ResolutionRoupBox->setGeometry(QRect(10, 90, 391, 101));
        ResolutionRoupBox->setFont(font);
        layoutWidget1 = new QWidget(ResolutionRoupBox);
        layoutWidget1->setObjectName(QStringLiteral("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(13, 31, 371, 61));
        verticalLayout_2 = new QVBoxLayout(layoutWidget1);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        ResolutionStyleLabel = new QLabel(layoutWidget1);
        ResolutionStyleLabel->setObjectName(QStringLiteral("ResolutionStyleLabel"));
        ResolutionStyleLabel->setFont(font1);

        verticalLayout_2->addWidget(ResolutionStyleLabel);

        VoxelSizeLabel = new QLabel(layoutWidget1);
        VoxelSizeLabel->setObjectName(QStringLiteral("VoxelSizeLabel"));
        VoxelSizeLabel->setSizeIncrement(QSize(3, 0));
        VoxelSizeLabel->setFont(font1);

        verticalLayout_2->addWidget(VoxelSizeLabel);

        layoutWidget = new QWidget(OpenImageDialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(10, 20, 391, 25));
        horizontalLayout_5 = new QHBoxLayout(layoutWidget);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(0, 0, 0, 0);
        FilePathDisplayLabel = new QLabel(layoutWidget);
        FilePathDisplayLabel->setObjectName(QStringLiteral("FilePathDisplayLabel"));

        horizontalLayout_5->addWidget(FilePathDisplayLabel);

        FilePathEdit = new QLineEdit(layoutWidget);
        FilePathEdit->setObjectName(QStringLiteral("FilePathEdit"));

        horizontalLayout_5->addWidget(FilePathEdit);

        OpenButton = new QPushButton(layoutWidget);
        OpenButton->setObjectName(QStringLiteral("OpenButton"));

        horizontalLayout_5->addWidget(OpenButton);

        layoutWidget_2 = new QWidget(OpenImageDialog);
        layoutWidget_2->setObjectName(QStringLiteral("layoutWidget_2"));
        layoutWidget_2->setGeometry(QRect(10, 50, 391, 25));
        horizontalLayout_6 = new QHBoxLayout(layoutWidget_2);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        horizontalLayout_6->setContentsMargins(0, 0, 0, 0);
        SomaPathDisplayLabel = new QLabel(layoutWidget_2);
        SomaPathDisplayLabel->setObjectName(QStringLiteral("SomaPathDisplayLabel"));

        horizontalLayout_6->addWidget(SomaPathDisplayLabel);

        StackPathEdit = new QLineEdit(layoutWidget_2);
        StackPathEdit->setObjectName(QStringLiteral("StackPathEdit"));

        horizontalLayout_6->addWidget(StackPathEdit);

        OpenStackButton = new QPushButton(layoutWidget_2);
        OpenStackButton->setObjectName(QStringLiteral("OpenStackButton"));

        horizontalLayout_6->addWidget(OpenStackButton);


        retranslateUi(OpenImageDialog);
        QObject::connect(OKButton, SIGNAL(clicked()), OpenImageDialog, SLOT(accept()));
        QObject::connect(CancelButton, SIGNAL(clicked()), OpenImageDialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(OpenImageDialog);
    } // setupUi

    void retranslateUi(QDialog *OpenImageDialog)
    {
        OpenImageDialog->setWindowTitle(QApplication::translate("OpenImageDialog", "OpenImage", 0));
        OKButton->setText(QApplication::translate("OpenImageDialog", "OK", 0));
        CancelButton->setText(QApplication::translate("OpenImageDialog", "Cancel", 0));
        FileInfoGroupBox->setTitle(QApplication::translate("OpenImageDialog", "FileInfo", 0));
        fileNameLabel->setText(QApplication::translate("OpenImageDialog", "File Name:", 0));
        getFileNameLable->setText(QApplication::translate("OpenImageDialog", "unknown File", 0));
        imageSizeLabel->setText(QApplication::translate("OpenImageDialog", "Image Size:", 0));
        getImageSizeLabel->setText(QApplication::translate("OpenImageDialog", "0 Slices, 0 * 0 , 0 bpp", 0));
        ResolutionRoupBox->setTitle(QApplication::translate("OpenImageDialog", "Resolution", 0));
        ResolutionStyleLabel->setText(QApplication::translate("OpenImageDialog", "Define : Voxel size (um)", 0));
        VoxelSizeLabel->setText(QApplication::translate("OpenImageDialog", "Voxel Size:", 0));
        FilePathDisplayLabel->setText(QApplication::translate("OpenImageDialog", "Image Path:", 0));
        OpenButton->setText(QApplication::translate("OpenImageDialog", "Open", 0));
        SomaPathDisplayLabel->setText(QApplication::translate("OpenImageDialog", "Stack Path:", 0));
        OpenStackButton->setText(QApplication::translate("OpenImageDialog", "Open", 0));
    } // retranslateUi

};

namespace Ui {
    class OpenImageDialog: public Ui_OpenImageDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OPENIMAGEDIALOG_H
>>>>>>> 57f5b9a669c6df3771c0c49b0e4d450a07ba81e0
