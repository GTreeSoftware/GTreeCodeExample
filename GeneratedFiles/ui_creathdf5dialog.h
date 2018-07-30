<<<<<<< HEAD
/********************************************************************************
** Form generated from reading UI file 'creathdf5dialog.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CREATHDF5DIALOG_H
#define UI_CREATHDF5DIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_creathdf5dialog
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout1;
    QLabel *label;
    QLineEdit *nameEdit;
    QPushButton *pathButton;
    QWidget *horizontalLayoutWidget_3;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_3;
    QWidget *horizontalLayoutWidget_4;
    QHBoxLayout *horizontalLayout_3;
    QDoubleSpinBox *xSpinBox;
    QWidget *horizontalLayoutWidget_5;
    QHBoxLayout *horizontalLayout_4;
    QDoubleSpinBox *ySpinBox;
    QWidget *horizontalLayoutWidget_6;
    QHBoxLayout *horizontalLayout_5;
    QDoubleSpinBox *zSpinBox;
    QGroupBox *FileInfoGroupBox;
    QGridLayout *gridLayout;
    QTextBrowser *infoBrowser;
    QLabel *infoLabel;
    QWidget *horizontalLayoutWidget_7;
    QHBoxLayout *horizontalLayout_9;
    QPushButton *okButton;
    QPushButton *cancelButton;
    QWidget *horizontalLayoutWidget_8;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_4;
    QWidget *horizontalLayoutWidget_9;
    QHBoxLayout *horizontalLayout_8;
    QLabel *label_2;
    QSpacerItem *horizontalSpacer;
    QLineEdit *pathEdit;
    QPushButton *SelectButton;

    void setupUi(QDialog *creathdf5dialog)
    {
        if (creathdf5dialog->objectName().isEmpty())
            creathdf5dialog->setObjectName(QStringLiteral("creathdf5dialog"));
        creathdf5dialog->resize(474, 420);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(creathdf5dialog->sizePolicy().hasHeightForWidth());
        creathdf5dialog->setSizePolicy(sizePolicy);
        horizontalLayoutWidget = new QWidget(creathdf5dialog);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(20, 90, 441, 31));
        horizontalLayout1 = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout1->setSpacing(6);
        horizontalLayout1->setContentsMargins(11, 11, 11, 11);
        horizontalLayout1->setObjectName(QStringLiteral("horizontalLayout1"));
        horizontalLayout1->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(horizontalLayoutWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setEnabled(true);
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy1);

        horizontalLayout1->addWidget(label);

        nameEdit = new QLineEdit(horizontalLayoutWidget);
        nameEdit->setObjectName(QStringLiteral("nameEdit"));
        nameEdit->setEnabled(false);
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(nameEdit->sizePolicy().hasHeightForWidth());
        nameEdit->setSizePolicy(sizePolicy2);
        nameEdit->setClearButtonEnabled(false);

        horizontalLayout1->addWidget(nameEdit);

        pathButton = new QPushButton(horizontalLayoutWidget);
        pathButton->setObjectName(QStringLiteral("pathButton"));
        QSizePolicy sizePolicy3(QSizePolicy::Maximum, QSizePolicy::Maximum);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(pathButton->sizePolicy().hasHeightForWidth());
        pathButton->setSizePolicy(sizePolicy3);
        pathButton->setAutoDefault(false);

        horizontalLayout1->addWidget(pathButton);

        horizontalLayoutWidget_3 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_3->setObjectName(QStringLiteral("horizontalLayoutWidget_3"));
        horizontalLayoutWidget_3->setGeometry(QRect(20, 130, 71, 31));
        horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget_3);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        label_3 = new QLabel(horizontalLayoutWidget_3);
        label_3->setObjectName(QStringLiteral("label_3"));

        horizontalLayout_2->addWidget(label_3);

        horizontalLayoutWidget_4 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_4->setObjectName(QStringLiteral("horizontalLayoutWidget_4"));
        horizontalLayoutWidget_4->setGeometry(QRect(100, 130, 60, 31));
        horizontalLayout_3 = new QHBoxLayout(horizontalLayoutWidget_4);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        xSpinBox = new QDoubleSpinBox(horizontalLayoutWidget_4);
        xSpinBox->setObjectName(QStringLiteral("xSpinBox"));

        horizontalLayout_3->addWidget(xSpinBox);

        horizontalLayoutWidget_5 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_5->setObjectName(QStringLiteral("horizontalLayoutWidget_5"));
        horizontalLayoutWidget_5->setGeometry(QRect(170, 130, 60, 31));
        horizontalLayout_4 = new QHBoxLayout(horizontalLayoutWidget_5);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, 0, 0, 0);
        ySpinBox = new QDoubleSpinBox(horizontalLayoutWidget_5);
        ySpinBox->setObjectName(QStringLiteral("ySpinBox"));

        horizontalLayout_4->addWidget(ySpinBox);

        horizontalLayoutWidget_6 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_6->setObjectName(QStringLiteral("horizontalLayoutWidget_6"));
        horizontalLayoutWidget_6->setGeometry(QRect(240, 130, 60, 31));
        horizontalLayout_5 = new QHBoxLayout(horizontalLayoutWidget_6);
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(0, 0, 0, 0);
        zSpinBox = new QDoubleSpinBox(horizontalLayoutWidget_6);
        zSpinBox->setObjectName(QStringLiteral("zSpinBox"));

        horizontalLayout_5->addWidget(zSpinBox);

        FileInfoGroupBox = new QGroupBox(creathdf5dialog);
        FileInfoGroupBox->setObjectName(QStringLiteral("FileInfoGroupBox"));
        FileInfoGroupBox->setGeometry(QRect(20, 180, 431, 201));
        QFont font;
        font.setFamily(QStringLiteral("Arial Black"));
        FileInfoGroupBox->setFont(font);
        gridLayout = new QGridLayout(FileInfoGroupBox);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        infoBrowser = new QTextBrowser(FileInfoGroupBox);
        infoBrowser->setObjectName(QStringLiteral("infoBrowser"));
        QSizePolicy sizePolicy4(QSizePolicy::Expanding, QSizePolicy::Maximum);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(infoBrowser->sizePolicy().hasHeightForWidth());
        infoBrowser->setSizePolicy(sizePolicy4);
        QFont font1;
        font1.setFamily(QStringLiteral("System"));
        infoBrowser->setFont(font1);

        gridLayout->addWidget(infoBrowser, 1, 0, 1, 1);

        infoLabel = new QLabel(FileInfoGroupBox);
        infoLabel->setObjectName(QStringLiteral("infoLabel"));
        infoLabel->setTextFormat(Qt::PlainText);

        gridLayout->addWidget(infoLabel, 0, 0, 1, 1);

        horizontalLayoutWidget_7 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_7->setObjectName(QStringLiteral("horizontalLayoutWidget_7"));
        horizontalLayoutWidget_7->setGeometry(QRect(300, 380, 161, 41));
        horizontalLayout_9 = new QHBoxLayout(horizontalLayoutWidget_7);
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        horizontalLayout_9->setContentsMargins(0, 0, 0, 0);
        okButton = new QPushButton(horizontalLayoutWidget_7);
        okButton->setObjectName(QStringLiteral("okButton"));
        okButton->setEnabled(true);
        okButton->setAutoDefault(false);

        horizontalLayout_9->addWidget(okButton);

        cancelButton = new QPushButton(horizontalLayoutWidget_7);
        cancelButton->setObjectName(QStringLiteral("cancelButton"));
        cancelButton->setAutoDefault(false);

        horizontalLayout_9->addWidget(cancelButton);

        horizontalLayoutWidget_8 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_8->setObjectName(QStringLiteral("horizontalLayoutWidget_8"));
        horizontalLayoutWidget_8->setGeometry(QRect(20, 10, 320, 31));
        horizontalLayout_6 = new QHBoxLayout(horizontalLayoutWidget_8);
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        horizontalLayout_6->setContentsMargins(0, 0, 0, 0);
        label_4 = new QLabel(horizontalLayoutWidget_8);
        label_4->setObjectName(QStringLiteral("label_4"));

        horizontalLayout_6->addWidget(label_4);

        horizontalLayoutWidget_9 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_9->setObjectName(QStringLiteral("horizontalLayoutWidget_9"));
        horizontalLayoutWidget_9->setGeometry(QRect(20, 50, 441, 31));
        horizontalLayout_8 = new QHBoxLayout(horizontalLayoutWidget_9);
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        horizontalLayout_8->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(horizontalLayoutWidget_9);
        label_2->setObjectName(QStringLiteral("label_2"));
        QSizePolicy sizePolicy5(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy5);

        horizontalLayout_8->addWidget(label_2);

        horizontalSpacer = new QSpacerItem(12, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout_8->addItem(horizontalSpacer);

        pathEdit = new QLineEdit(horizontalLayoutWidget_9);
        pathEdit->setObjectName(QStringLiteral("pathEdit"));
        sizePolicy2.setHeightForWidth(pathEdit->sizePolicy().hasHeightForWidth());
        pathEdit->setSizePolicy(sizePolicy2);

        horizontalLayout_8->addWidget(pathEdit);

        SelectButton = new QPushButton(horizontalLayoutWidget_9);
        SelectButton->setObjectName(QStringLiteral("SelectButton"));
        sizePolicy.setHeightForWidth(SelectButton->sizePolicy().hasHeightForWidth());
        SelectButton->setSizePolicy(sizePolicy);
        SelectButton->setAutoDefault(false);

        horizontalLayout_8->addWidget(SelectButton);


        retranslateUi(creathdf5dialog);

        cancelButton->setDefault(false);


        QMetaObject::connectSlotsByName(creathdf5dialog);
    } // setupUi

    void retranslateUi(QDialog *creathdf5dialog)
    {
        creathdf5dialog->setWindowTitle(QApplication::translate("creathdf5dialog", "creathdf5dialog", 0));
        label->setText(QApplication::translate("creathdf5dialog", "Create File:", 0));
        pathButton->setText(QApplication::translate("creathdf5dialog", "Open", 0));
        label_3->setText(QApplication::translate("creathdf5dialog", "Resolution:", 0));
        FileInfoGroupBox->setTitle(QApplication::translate("creathdf5dialog", "FileInfo", 0));
        infoLabel->setText(QString());
        okButton->setText(QApplication::translate("creathdf5dialog", "OK", 0));
        cancelButton->setText(QApplication::translate("creathdf5dialog", "Cancel", 0));
        label_4->setText(QApplication::translate("creathdf5dialog", "Please input the path and name for the generate file:", 0));
        label_2->setText(QApplication::translate("creathdf5dialog", "Save Path:", 0));
        SelectButton->setText(QApplication::translate("creathdf5dialog", "Select Tiff", 0));
    } // retranslateUi

};

namespace Ui {
    class creathdf5dialog: public Ui_creathdf5dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CREATHDF5DIALOG_H
=======
/********************************************************************************
** Form generated from reading UI file 'creathdf5dialog.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CREATHDF5DIALOG_H
#define UI_CREATHDF5DIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_creathdf5dialog
{
public:
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout1;
    QLabel *label;
    QLineEdit *nameEdit;
    QPushButton *pathButton;
    QWidget *horizontalLayoutWidget_3;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_3;
    QWidget *horizontalLayoutWidget_4;
    QHBoxLayout *horizontalLayout_3;
    QDoubleSpinBox *xSpinBox;
    QWidget *horizontalLayoutWidget_5;
    QHBoxLayout *horizontalLayout_4;
    QDoubleSpinBox *ySpinBox;
    QWidget *horizontalLayoutWidget_6;
    QHBoxLayout *horizontalLayout_5;
    QDoubleSpinBox *zSpinBox;
    QGroupBox *FileInfoGroupBox;
    QGridLayout *gridLayout;
    QTextBrowser *infoBrowser;
    QLabel *infoLabel;
    QWidget *horizontalLayoutWidget_7;
    QHBoxLayout *horizontalLayout_9;
    QPushButton *okButton;
    QPushButton *cancelButton;
    QWidget *horizontalLayoutWidget_8;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_4;
    QWidget *horizontalLayoutWidget_9;
    QHBoxLayout *horizontalLayout_8;
    QLabel *label_2;
    QSpacerItem *horizontalSpacer;
    QLineEdit *pathEdit;
    QPushButton *SelectButton;

    void setupUi(QDialog *creathdf5dialog)
    {
        if (creathdf5dialog->objectName().isEmpty())
            creathdf5dialog->setObjectName(QStringLiteral("creathdf5dialog"));
        creathdf5dialog->resize(474, 420);
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(creathdf5dialog->sizePolicy().hasHeightForWidth());
        creathdf5dialog->setSizePolicy(sizePolicy);
        horizontalLayoutWidget = new QWidget(creathdf5dialog);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(20, 90, 441, 31));
        horizontalLayout1 = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout1->setSpacing(6);
        horizontalLayout1->setContentsMargins(11, 11, 11, 11);
        horizontalLayout1->setObjectName(QStringLiteral("horizontalLayout1"));
        horizontalLayout1->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(horizontalLayoutWidget);
        label->setObjectName(QStringLiteral("label"));
        label->setEnabled(true);
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy1);

        horizontalLayout1->addWidget(label);

        nameEdit = new QLineEdit(horizontalLayoutWidget);
        nameEdit->setObjectName(QStringLiteral("nameEdit"));
        nameEdit->setEnabled(false);
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(nameEdit->sizePolicy().hasHeightForWidth());
        nameEdit->setSizePolicy(sizePolicy2);
        nameEdit->setClearButtonEnabled(false);

        horizontalLayout1->addWidget(nameEdit);

        pathButton = new QPushButton(horizontalLayoutWidget);
        pathButton->setObjectName(QStringLiteral("pathButton"));
        QSizePolicy sizePolicy3(QSizePolicy::Maximum, QSizePolicy::Maximum);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(pathButton->sizePolicy().hasHeightForWidth());
        pathButton->setSizePolicy(sizePolicy3);
        pathButton->setAutoDefault(false);

        horizontalLayout1->addWidget(pathButton);

        horizontalLayoutWidget_3 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_3->setObjectName(QStringLiteral("horizontalLayoutWidget_3"));
        horizontalLayoutWidget_3->setGeometry(QRect(20, 130, 71, 31));
        horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget_3);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        label_3 = new QLabel(horizontalLayoutWidget_3);
        label_3->setObjectName(QStringLiteral("label_3"));

        horizontalLayout_2->addWidget(label_3);

        horizontalLayoutWidget_4 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_4->setObjectName(QStringLiteral("horizontalLayoutWidget_4"));
        horizontalLayoutWidget_4->setGeometry(QRect(100, 130, 60, 31));
        horizontalLayout_3 = new QHBoxLayout(horizontalLayoutWidget_4);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        xSpinBox = new QDoubleSpinBox(horizontalLayoutWidget_4);
        xSpinBox->setObjectName(QStringLiteral("xSpinBox"));

        horizontalLayout_3->addWidget(xSpinBox);

        horizontalLayoutWidget_5 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_5->setObjectName(QStringLiteral("horizontalLayoutWidget_5"));
        horizontalLayoutWidget_5->setGeometry(QRect(170, 130, 60, 31));
        horizontalLayout_4 = new QHBoxLayout(horizontalLayoutWidget_5);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, 0, 0, 0);
        ySpinBox = new QDoubleSpinBox(horizontalLayoutWidget_5);
        ySpinBox->setObjectName(QStringLiteral("ySpinBox"));

        horizontalLayout_4->addWidget(ySpinBox);

        horizontalLayoutWidget_6 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_6->setObjectName(QStringLiteral("horizontalLayoutWidget_6"));
        horizontalLayoutWidget_6->setGeometry(QRect(240, 130, 60, 31));
        horizontalLayout_5 = new QHBoxLayout(horizontalLayoutWidget_6);
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(0, 0, 0, 0);
        zSpinBox = new QDoubleSpinBox(horizontalLayoutWidget_6);
        zSpinBox->setObjectName(QStringLiteral("zSpinBox"));

        horizontalLayout_5->addWidget(zSpinBox);

        FileInfoGroupBox = new QGroupBox(creathdf5dialog);
        FileInfoGroupBox->setObjectName(QStringLiteral("FileInfoGroupBox"));
        FileInfoGroupBox->setGeometry(QRect(20, 180, 431, 201));
        QFont font;
        font.setFamily(QStringLiteral("Arial Black"));
        FileInfoGroupBox->setFont(font);
        gridLayout = new QGridLayout(FileInfoGroupBox);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        infoBrowser = new QTextBrowser(FileInfoGroupBox);
        infoBrowser->setObjectName(QStringLiteral("infoBrowser"));
        QSizePolicy sizePolicy4(QSizePolicy::Expanding, QSizePolicy::Maximum);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(infoBrowser->sizePolicy().hasHeightForWidth());
        infoBrowser->setSizePolicy(sizePolicy4);
        QFont font1;
        font1.setFamily(QStringLiteral("System"));
        infoBrowser->setFont(font1);

        gridLayout->addWidget(infoBrowser, 1, 0, 1, 1);

        infoLabel = new QLabel(FileInfoGroupBox);
        infoLabel->setObjectName(QStringLiteral("infoLabel"));
        infoLabel->setTextFormat(Qt::PlainText);

        gridLayout->addWidget(infoLabel, 0, 0, 1, 1);

        horizontalLayoutWidget_7 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_7->setObjectName(QStringLiteral("horizontalLayoutWidget_7"));
        horizontalLayoutWidget_7->setGeometry(QRect(300, 380, 161, 41));
        horizontalLayout_9 = new QHBoxLayout(horizontalLayoutWidget_7);
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        horizontalLayout_9->setContentsMargins(0, 0, 0, 0);
        okButton = new QPushButton(horizontalLayoutWidget_7);
        okButton->setObjectName(QStringLiteral("okButton"));
        okButton->setEnabled(true);
        okButton->setAutoDefault(false);

        horizontalLayout_9->addWidget(okButton);

        cancelButton = new QPushButton(horizontalLayoutWidget_7);
        cancelButton->setObjectName(QStringLiteral("cancelButton"));
        cancelButton->setAutoDefault(false);

        horizontalLayout_9->addWidget(cancelButton);

        horizontalLayoutWidget_8 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_8->setObjectName(QStringLiteral("horizontalLayoutWidget_8"));
        horizontalLayoutWidget_8->setGeometry(QRect(20, 10, 320, 31));
        horizontalLayout_6 = new QHBoxLayout(horizontalLayoutWidget_8);
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        horizontalLayout_6->setContentsMargins(0, 0, 0, 0);
        label_4 = new QLabel(horizontalLayoutWidget_8);
        label_4->setObjectName(QStringLiteral("label_4"));

        horizontalLayout_6->addWidget(label_4);

        horizontalLayoutWidget_9 = new QWidget(creathdf5dialog);
        horizontalLayoutWidget_9->setObjectName(QStringLiteral("horizontalLayoutWidget_9"));
        horizontalLayoutWidget_9->setGeometry(QRect(20, 50, 441, 31));
        horizontalLayout_8 = new QHBoxLayout(horizontalLayoutWidget_9);
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        horizontalLayout_8->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(horizontalLayoutWidget_9);
        label_2->setObjectName(QStringLiteral("label_2"));
        QSizePolicy sizePolicy5(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(label_2->sizePolicy().hasHeightForWidth());
        label_2->setSizePolicy(sizePolicy5);

        horizontalLayout_8->addWidget(label_2);

        horizontalSpacer = new QSpacerItem(12, 20, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout_8->addItem(horizontalSpacer);

        pathEdit = new QLineEdit(horizontalLayoutWidget_9);
        pathEdit->setObjectName(QStringLiteral("pathEdit"));
        sizePolicy2.setHeightForWidth(pathEdit->sizePolicy().hasHeightForWidth());
        pathEdit->setSizePolicy(sizePolicy2);

        horizontalLayout_8->addWidget(pathEdit);

        SelectButton = new QPushButton(horizontalLayoutWidget_9);
        SelectButton->setObjectName(QStringLiteral("SelectButton"));
        sizePolicy.setHeightForWidth(SelectButton->sizePolicy().hasHeightForWidth());
        SelectButton->setSizePolicy(sizePolicy);
        SelectButton->setAutoDefault(false);

        horizontalLayout_8->addWidget(SelectButton);


        retranslateUi(creathdf5dialog);

        cancelButton->setDefault(false);


        QMetaObject::connectSlotsByName(creathdf5dialog);
    } // setupUi

    void retranslateUi(QDialog *creathdf5dialog)
    {
        creathdf5dialog->setWindowTitle(QApplication::translate("creathdf5dialog", "creathdf5dialog", 0));
        label->setText(QApplication::translate("creathdf5dialog", "Create File:", 0));
        pathButton->setText(QApplication::translate("creathdf5dialog", "Open", 0));
        label_3->setText(QApplication::translate("creathdf5dialog", "Resolution:", 0));
        FileInfoGroupBox->setTitle(QApplication::translate("creathdf5dialog", "FileInfo", 0));
        infoLabel->setText(QString());
        okButton->setText(QApplication::translate("creathdf5dialog", "OK", 0));
        cancelButton->setText(QApplication::translate("creathdf5dialog", "Cancel", 0));
        label_4->setText(QApplication::translate("creathdf5dialog", "Please input the path and name for the generate file:", 0));
        label_2->setText(QApplication::translate("creathdf5dialog", "Save Path:", 0));
        SelectButton->setText(QApplication::translate("creathdf5dialog", "Select Tiff", 0));
    } // retranslateUi

};

namespace Ui {
    class creathdf5dialog: public Ui_creathdf5dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CREATHDF5DIALOG_H
>>>>>>> 57f5b9a669c6df3771c0c49b0e4d450a07ba81e0
