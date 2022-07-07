import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox
import ExcelWindow
import dxfTxt
import json


class MyMainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.main_ui = ExcelWindow.Ui_MainWindow()
        self.main_ui.setupUi(self)
        with open("setting.json") as f:
            args = json.load(f)
        self.main_ui.kLineEdit.setText(str(args['k']))
        self.main_ui.stepLineEdit.setText(str(args['step']))
        self.main_ui.zMaxLineEdit.setText(str(args['zMax']))
        self.main_ui.xBackLineEdit.setText(str(args['xBack']))
        self.main_ui.yBackLineEdit.setText(str(args['yBack']))
        self.main_ui.zBackLineEdit.setText(str(args['zBack']))
        self.main_ui.dLineEdit.setText(str(args['delta']))
        self.main_ui.d2LineEdit.setText(str(args['delta2']))

    def chooseExcel(self):
        path = QtWidgets.QFileDialog.getOpenFileName(
            None, "选取dxf文件", "./", "ALL(*.*);;dxf文件(*.dxf)", "dxf文件(*.dxf)")  # 选择excel路径
        # print(path[0])
        self.main_ui.excelPath.setText(path[0])
        # self.main_ui.textEdit.setText(path)

    def start(self):
        try:
            path = self.main_ui.excelPath.text()
            if path == '':
                QMessageBox().about(None, "警告", '请先选择dxf文件')
                return
            name = path.split('/')[-1].split('.')[0]
            args = {'path': path,
                    'name': name,
                    # str = 'F:\\系统默认\\Desktop\\compression-12-2.xlsx'  # 文件路径
                    'k': float(self.main_ui.kLineEdit.text()),  # 0.2   # 倍率
                    # 0.25  # z 每次增加量
                    'step': float(self.main_ui.stepLineEdit.text()),
                    # 38     # z 最大值
                    'zMax': float(self.main_ui.zMaxLineEdit.text()),
                    'xBack': float(self.main_ui.xBackLineEdit.text()),
                    'yBack': float(self.main_ui.yBackLineEdit.text()),
                    'zBack': float(self.main_ui.zBackLineEdit.text()),
                    'delta': float(self.main_ui.dLineEdit.text()),
                    'delta2': float(self.main_ui.d2LineEdit.text())}
            with open("setting.json", "w") as f:
                f.write(json.dumps(args))
            dxfTxt.start(args)
            QMessageBox().about(None, "结果", f"生成成功，已保存到D:\\{name}.gcode")
        except:
            import traceback
            exc_type, exc_value, exc_traceback = sys.exc_info()
            error = str(repr(traceback.format_exception(
                exc_type, exc_value, exc_traceback)))  # 将异常信息转为字符串
            print(error)
            QMessageBox().about(None, "出错", error)


# pyinstaller -p D:\DevTools\Python38\Lib\site-packages\cvxopt\.lib -D   dxf2gcode.py
if __name__ == '__main__':
    myapp = QApplication(sys.argv)
    mainWindow = MyMainWindow()
    mainWindow.show()
    sys.exit(myapp.exec_())
