#pragma once


// ControlNetworkDrawingDlg 对话框

class ControlNetworkDrawingDlg : public CDialogEx
{
	DECLARE_DYNAMIC(ControlNetworkDrawingDlg)

public:
	ControlNetworkDrawingDlg(CWnd* pParent = nullptr);   // 标准构造函数
	virtual ~ControlNetworkDrawingDlg();

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_CONTROLNETWORKDRAWING_DIALOG };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnPaint();

//实现
public:
	ControlNetworkAdjust* CNA_Draw;//控制网平差类对象
};
