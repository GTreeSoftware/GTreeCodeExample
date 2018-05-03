#ifndef NGTREEWIDGET_H
#define NGTREEWIDGET_H
#include <QTreeWidget>
#include <QAction>
#include "../ngtypes/ParamPack.h"
#include "../ngtypes/tree.h"
class QTreeWidgetItem;

class NGTreeWidget :
	public QTreeWidget
{
	Q_OBJECT
public:
	enum WIDGETMODE{ TREEWIDGETMODE, CHECKWIDGETMODE } widgetMode_;

	NGTreeWidget(QWidget *par = nullptr, WIDGETMODE = TREEWIDGETMODE);
	~NGTreeWidget();

	//enum ITEMTYPE{TREETYPE, CHECKTYPE};
	void addRoot(const QStringList&);
	void SetParam(NGParamPack arg){ paramPack_ = arg; }


public slots :
	void ActivateTree_Slot();
	void CompareTree_Slot();
	void ToggleTreeVisible_Slot();
	void GotoDiff_Slot();
	void TriggeredCheckState_Slot();
	void DeleteAnnotate_Slot();
	void ToggleShowDiffArrow_Slot(bool);
	void ResetTraverse_Slot();
	void ToggleShowTraverseFlag_Slot();

protected:
    void CreateAction();
    virtual void contextMenuEvent(QContextMenuEvent *ev);

signals:
    void ActiveTree_Signal(int);
    void CompareTree_Signal(int);
    void ToggleTreeVisible_Signal(int);
    void GotoDiff_Signal(int);
    void ToggleShowDiffArrow_Signal(bool);
    void ResetTraverse_Signal();
    void ToggleShowTraverseFlag_Signal();

private:
    QTreeWidgetItem *choosedTreeWidgetItem;
    QAction *pActivate_, *pCompare_, *pToggleTreeShow_, *pGotoDiff_, *pToggleCheckState_, *pDeleteAnnotate_, *pToggleShowDiffArrow_,
        *pClearTraverseRecord_, *pToggleShowTraverseFlag_;
    NGParamPack paramPack_;
};

#endif