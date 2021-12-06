// ControlNetworkDrawingDlg.cpp: 实现文件
//

#include "pch.h"
#include "ControlNetwork.h"
#include "ControlNetworkDrawingDlg.h"
#include "afxdialogex.h"


// ControlNetworkDrawingDlg 对话框

IMPLEMENT_DYNAMIC(ControlNetworkDrawingDlg, CDialogEx)

ControlNetworkDrawingDlg::ControlNetworkDrawingDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_CONTROLNETWORKDRAWING_DIALOG, pParent)
{

}

ControlNetworkDrawingDlg::~ControlNetworkDrawingDlg()
{
}

void ControlNetworkDrawingDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}


BEGIN_MESSAGE_MAP(ControlNetworkDrawingDlg, CDialogEx)
	ON_WM_PAINT()
END_MESSAGE_MAP()


// ControlNetworkDrawingDlg 消息处理程序


void ControlNetworkDrawingDlg::OnPaint()
{
	CPaintDC dc(this); // device context for painting
					   // TODO: 在此处添加消息处理程序代码
					   // 不为绘图消息调用 CDialogEx::OnPaint()

	CWnd* pWin = GetDlgItem(IDC_CTRLNET);//获取Picture控件的指针
	CRect rect;
	pWin->GetClientRect(rect);//把控件的长宽、坐标等信息保存在rect里
	CDC* pDC = pWin->GetDC();//获取该控件的画布

	CNA_Draw->ControlNetworkDrawing(pDC, rect);//画图
}

