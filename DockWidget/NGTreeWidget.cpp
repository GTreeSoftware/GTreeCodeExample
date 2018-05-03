
#include <QContextMenuEvent>
#include <QMenu>
#include <QMessageBox>
#include "NGTreeWidget.h"


NGTreeWidget::NGTreeWidget(QWidget *par, WIDGETMODE arg) :QTreeWidget(par), widgetMode_(arg)
{
    CreateAction();
}


NGTreeWidget::~NGTreeWidget()
{
}

void NGTreeWidget::CreateAction()
{
    if (widgetMode_ == TREEWIDGETMODE) {
        pToggleTreeShow_ = new QAction("Toggle_Visible", this); //pToggleTreeShow_->setCheckable(true); pToggleTreeShow_->setChecked(true);
        connect(pToggleTreeShow_, SIGNAL(triggered()), this, SLOT(ToggleTreeVisible_Slot()));
        pActivate_ = new QAction("Activate", this);
        connect(pActivate_, SIGNAL(triggered()), this, SLOT(ActivateTree_Slot()));
        pCompare_ = new QAction("Compare", this);
        connect(pCompare_, SIGNAL(triggered()), this, SLOT(CompareTree_Slot()));
        pClearTraverseRecord_ = new QAction("Reset_Traverse", this);
        connect(pClearTraverseRecord_, SIGNAL(triggered()), this, SLOT(ResetTraverse_Slot()));
        pToggleShowTraverseFlag_ = new QAction("Toggle_TraFlag",this);
        connect(pToggleShowTraverseFlag_, SIGNAL(triggered()), this, SLOT(ToggleShowTraverseFlag_Slot()));

    }
    if (widgetMode_ == CHECKWIDGETMODE) {
        pGotoDiff_ = new QAction("Goto", this);
        connect(pGotoDiff_, SIGNAL(triggered()), this, SLOT(GotoDiff_Slot()));
        pToggleCheckState_ = new QAction("Toggle_Check", this);
        connect(pToggleCheckState_, SIGNAL(triggered()), this, SLOT(TriggeredCheckState_Slot()));
        pDeleteAnnotate_ = new QAction("Del_Annotate",this);
        connect(pDeleteAnnotate_, SIGNAL(triggered()), this, SLOT(DeleteAnnotate_Slot()));
        pToggleShowDiffArrow_ = new QAction("Show_Diff_Arrow", this);
        pToggleShowDiffArrow_->setCheckable(true); pToggleShowDiffArrow_->setChecked(false);
        connect(pToggleShowDiffArrow_, SIGNAL(toggled(bool)), this, SLOT(ToggleShowDiffArrow_Slot(bool)));
    }
}

void NGTreeWidget::contextMenuEvent(QContextMenuEvent *ev)
{
    choosedTreeWidgetItem = this->currentItem();
    if (!choosedTreeWidgetItem) 
        return;
    QMenu *menu = new QMenu(this);
    switch (widgetMode_)
    {
    case NGTreeWidget::TREEWIDGETMODE:
        menu->addAction(pToggleTreeShow_);
        menu->addAction(pActivate_);
        menu->addAction(pCompare_);
        menu->addAction(pClearTraverseRecord_);
        menu->addAction(pToggleShowTraverseFlag_);
        break;
    case NGTreeWidget::CHECKWIDGETMODE:
        menu->addAction(pGotoDiff_);
        menu->addAction(pToggleCheckState_);
        menu->addAction(pDeleteAnnotate_);
        menu->addAction(pToggleShowDiffArrow_);
        break;
    default:
        break;
    }
    menu->exec(QCursor::pos());
    delete menu;
    QWidget::contextMenuEvent(ev);
}

void NGTreeWidget::addRoot(const QStringList &itemStr)
{
    switch (widgetMode_)
    {
    case NGTreeWidget::TREEWIDGETMODE:
        if (itemStr.size() != 2lu) return;
        break;
    case NGTreeWidget::CHECKWIDGETMODE:
        if (itemStr.size() != 4lu) return;
        break;
    default:
        break;
    }
    QTreeWidgetItem *newIt = new QTreeWidgetItem(this);
    if (widgetMode_ == TREEWIDGETMODE) {
        newIt->setText(0, itemStr[0]);
        newIt->setText(1, itemStr[1]);
    }
    if (widgetMode_ == CHECKWIDGETMODE) {
        newIt->setText(0, itemStr[0]);
        newIt->setText(1, itemStr[1]);
        newIt->setText(2, itemStr[2]);
        newIt->setText(3, itemStr[3]);
    }
    this->addTopLevelItem(newIt);
}

void NGTreeWidget::ActivateTree_Slot()
{
    QString idStr = choosedTreeWidgetItem->text(0);
    idStr = idStr.mid(4);
    int id = idStr.toInt();
    emit ActiveTree_Signal(id);
}

void NGTreeWidget::CompareTree_Slot()
{
    if (choosedTreeWidgetItem->checkState(0) == Qt::Checked) {
        QMessageBox::information(this, "cannot compare", "Current tree is activated. Please choose another tree.");
        return;
    }
    QString idStr = choosedTreeWidgetItem->text(0);
    idStr = idStr.mid(4);
    int id = idStr.toInt();
    emit CompareTree_Signal(id);
}

void NGTreeWidget::ToggleTreeVisible_Slot()
{
    QString idStr = choosedTreeWidgetItem->text(0);
    idStr = idStr.mid(4);
    int id = idStr.toInt();
    emit ToggleTreeVisible_Signal(id);
}

void NGTreeWidget::GotoDiff_Slot()
{
    int id = this->indexOfTopLevelItem(choosedTreeWidgetItem);
    emit GotoDiff_Signal(id);
}

void NGTreeWidget::TriggeredCheckState_Slot()
{
    int index = this->indexOfTopLevelItem(choosedTreeWidgetItem);
    auto &stateList = paramPack_->activeTree->m_suspiciousCheckState_;
    if (choosedTreeWidgetItem->text(0) == "no") {
        choosedTreeWidgetItem->setText(0, "yes");
        stateList[index] = 1;
    }
    else if (choosedTreeWidgetItem->text(0) == "yes"){
        choosedTreeWidgetItem->setText(0, "tol");
        stateList[index] = 2;
    }
    else{
        choosedTreeWidgetItem->setText(0, "no");
        stateList[index] = 0;
    }
}

void NGTreeWidget::DeleteAnnotate_Slot()
{
    if (!paramPack_->activeTree) return;
    Vec3d pos;
    QString idStr = choosedTreeWidgetItem->text(1);
    pos(0)= idStr.toDouble();
    idStr = choosedTreeWidgetItem->text(2);
    pos(1) = idStr.toDouble();
    idStr = choosedTreeWidgetItem->text(3);
    pos(2) = idStr.toDouble();
    auto &sp = paramPack_->activeTree->m_suspiciousPositions_;
    auto &sc = paramPack_->activeTree->m_suspiciousCheckState_;
    auto it = std::find_if(sp.begin(), sp.end(), [&](const Vec3d &arg){ return (arg - pos).norm() < 2.0; });
    auto dist = std::distance(sp.begin(), it);
    sp.erase(it);
    sc.erase(sc.begin() + dist);
    int index = this->indexOfTopLevelItem(choosedTreeWidgetItem);
    this->takeTopLevelItem(index);
}

void NGTreeWidget::ToggleShowDiffArrow_Slot(bool arg)
{
    emit ToggleShowDiffArrow_Signal(arg);
}

void NGTreeWidget::ResetTraverse_Slot()
{
    //int index = this->indexOfTopLevelItem(choosedTreeWidgetItem);
    emit ResetTraverse_Signal();
}

void NGTreeWidget::ToggleShowTraverseFlag_Slot()
{
    emit ToggleShowTraverseFlag_Signal();
}
