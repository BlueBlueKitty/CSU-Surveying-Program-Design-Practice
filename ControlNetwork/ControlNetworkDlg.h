
// ControlNetworkDlg.h: 头文件
//

#pragma once


// CControlNetworkDlg 对话框
class CControlNetworkDlg : public CDialogEx
{
// 构造
public:
	CControlNetworkDlg(CWnd* pParent = nullptr);	// 标准构造函数

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_CONTROLNETWORK_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
private:
	ControlNetworkAdjust CNA_Cal;//定义控制网平差类对象
	ControlNetworkDrawingDlg CND;//定义控制网绘制对话框类
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	friend void ReadFileContent(CStdioFile& sf,CString& strFileContent);//从文件中读取数据到CString字符串中
	afx_msg void OnBnClickedReadfiledata();
	CString strData;
	CString strOutData;
	afx_msg void OnBnClickedCalapprocoord();
	afx_msg void OnBnClickedAdjustcal();
	afx_msg void OnBnClickedNetdrawing();
};
