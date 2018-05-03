#ifndef TREECHECKER_H
#define TREECHECKER_H
#include <QDialog>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QSpinBox>
#include <QLineEdit>
#include <QPushButton>
#include <QCheckBox>
#include "ngtypes/ParamPack.h"
class TreeChecker :public QDialog
{
    Q_OBJECT
public:
    explicit TreeChecker(QWidget *parent = 0);
    ~TreeChecker();
    void SetParam(NGParamPack arg);
    QString TranslateLayerListIntoText(const std::vector<int>& layerList) ;
    bool TranslateTextIntoLayerList(const QString &str, std::vector<int>& layerList) ;
    void SetHasSoma(bool);
    void SetRadius(int arg);
    void SetLayerID(const std::vector<int>&);
    void SetBranchID(const std::vector<int>&);
    bool GetLayerID(std::vector<int>&) ;
    bool GetBranchID(std::vector<int>&) ;
    int GetRadius()const{ return inputRadius->value(); }
    bool HasSoma() const { return checkBox->isChecked(); }
 private slots:
    void UpdateTreeCheckerTips_Slot();
private:
    QVBoxLayout *vl;
    QHBoxLayout *hl1;
    QHBoxLayout *hl2;
    QHBoxLayout *hl3;
    QLabel *tips;
    QLabel *inputLayerLabel;
    QLabel *inputBranchLabel ;
    QLabel *inputRadiusLabel;
    QSpinBox *inputRadius;
    QLineEdit *inputLayerList ;
    QLineEdit *inputBranchList;
    QHBoxLayout *hl ;
    QPushButton *acceptButton ;
    QPushButton *cancelButton ;
    QCheckBox *checkBox ;

    NGParamPack paramPack;
};
#endif