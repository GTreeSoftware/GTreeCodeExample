
#include <QMessageBox>
#include "TreeChecker.h"
#include "ngtypes/ParamPack.h"
#include "ngtypes/tree.h"

TreeChecker::TreeChecker(QWidget *parent /*= 0*/) : QDialog(parent)
{
    vl = new QVBoxLayout(this);
    this->setLayout(vl);
    hl1 = new QHBoxLayout();//layer
    hl2 = new QHBoxLayout();//branch
    hl3 = new QHBoxLayout();//radius
    vl->addLayout(hl1); vl->addLayout(hl2); vl->addLayout(hl3);
    tips = new QLabel(this); 
    vl->addWidget(tips);
    inputLayerLabel = new QLabel(this); inputLayerLabel->setText("input layer");
    inputBranchLabel = new QLabel(this); inputBranchLabel->setText("input branch");
    inputRadiusLabel = new QLabel(this); inputRadiusLabel->setText("input radius");
    inputRadius = new QSpinBox(this); 
    inputLayerList = new QLineEdit(this);
    inputBranchList = new QLineEdit(this);
    connect(inputLayerList, SIGNAL(editingFinished()), this, SLOT(UpdateTreeCheckerTips_Slot()));
    inputLayerList->setMaximumHeight(inputRadius->height());
    
    inputBranchList->setMaximumHeight(inputRadius->height());
   
    hl1->addWidget(inputLayerLabel); hl1->addWidget(inputLayerList);
    hl2->addWidget(inputBranchLabel); hl2->addWidget(inputBranchList);
    hl3->addWidget(inputRadiusLabel); hl3->addWidget(inputRadius);
    QHBoxLayout *hl = new QHBoxLayout();
    vl->addLayout(hl);
    acceptButton = new QPushButton("accept", this);
    cancelButton = new QPushButton("cancel", this);
    hl->addWidget(acceptButton);
    hl->addWidget(cancelButton);
    checkBox = new QCheckBox("Has Soma", this); 
    vl->addWidget(checkBox);
    QObject::connect(acceptButton, SIGNAL(clicked()), this, SLOT(accept()));
    QObject::connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
    inputLayerList->setFocus();
}

TreeChecker::~TreeChecker()
{

}

QString TreeChecker::TranslateLayerListIntoText(const std::vector<int>& layerList) 
{
    QString tmp;
    if (!layerList.empty()) {
        for (size_t i = 0; i < layerList.size() - 1llu; ++i) {
            tmp += QString::number(layerList[i]) + ',';
        }
        tmp += QString::number(layerList.back());
    }
    return tmp;
}

bool TreeChecker::TranslateTextIntoLayerList(const QString &str, std::vector<int>& layerList) 
{
#define FAILURE_CALLBACK if (!ok) { \
    layerList.clear(); \
    QMessageBox::warning(this, "error in layer list", "You do not input layer nums.\n"); \
    return false; }
    //-----------------------------------------
    layerList.clear();
    if (str.isEmpty()) return false;
    bool ok;
    QStringList splitList = str.split(',');
    QStringList tmp;
    for (int i = 0; i < splitList.size(); ++i) {
        tmp = splitList[i].split('-');
        if (tmp.size() < 2) {
            layerList.push_back(splitList[i].toInt(&ok));
            FAILURE_CALLBACK;
        }
        else if (tmp.size() > 2) {
            FAILURE_CALLBACK;
        }
        else{
            int beg = tmp[0].toInt(&ok);
            FAILURE_CALLBACK;
            int end = tmp[1].toInt(&ok);
            FAILURE_CALLBACK;
            for (int k = beg; k <= end; ++k)
                layerList.push_back(k);
        }
    }
    return true;
}

void TreeChecker::SetRadius(int arg)
{
    inputRadius->setValue(arg);
}

void TreeChecker::SetLayerID(const std::vector<int>& layerDisplay_)
{
    inputLayerList->setText(TranslateLayerListIntoText(layerDisplay_));
}

void TreeChecker::SetBranchID(const std::vector<int>& branchDisplay_)
{
    //branch start from 0 , instead of 1
    std::vector<int> tmpBranchList = branchDisplay_; std::for_each(tmpBranchList.begin(), tmpBranchList.end(), [](int &arg){++arg; });
    inputBranchList->setText(TranslateLayerListIntoText(tmpBranchList));
    
}

void TreeChecker::SetParam(NGParamPack arg)
{
    paramPack = arg;
    int maxD = paramPack->activeTree->m_Topo.max_depth();
    tips->setText(tr("The max level depth of current tree is %1").arg(maxD-1));
}

void TreeChecker::SetHasSoma(bool hasSoma_)
{
    checkBox->setChecked(hasSoma_);
}

bool TreeChecker::GetLayerID(std::vector<int>& tmpLayer) 
{
    if (!TranslateTextIntoLayerList(inputLayerList->text(), tmpLayer)) return false;
    return true;
}

bool TreeChecker::GetBranchID(std::vector<int>& tmpBranch) 
{
    if (!TranslateTextIntoLayerList(inputBranchList->text(), tmpBranch)) return false;
    return true;
}

void TreeChecker::UpdateTreeCheckerTips_Slot()
{
    std::vector<int> tmpLayer;
    if (!GetLayerID(tmpLayer)) return;
    std::sort(tmpLayer.begin(), tmpLayer.end());
    auto level = tmpLayer.front();
    auto &topo = paramPack->activeTree->m_Topo;
    int branchNum = 0;
    for (auto iter = topo.begin() + 1; iter != topo.end(); ++iter) {
        if (topo.depth(iter) == level) {
            ++branchNum;
        }
    }
    tips->setText(tr("Max level is %1, top level branch num is %2").arg(paramPack->activeTree->m_Topo.max_depth()-1).arg(branchNum));
}