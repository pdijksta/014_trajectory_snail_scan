from PyQt5.QtCore import QObject, pyqtSignal

class WorkerBase(QObject):
    """
    Wraps a function call such that it can be used by QThread.
    """
    finished = pyqtSignal()
    progress = pyqtSignal(int)

    def __init__(self, *args, **kwargs):
        QObject.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self.outp = None
        self.error = False
        self.abort = False
        self.result_dict = {}

    def run(self):
        try:
            self.outp = self.func(*self.args, **self.kwargs)
        except Exception as e:
            self.error = True
            print('Worker failed with following error message:\n%s' % e)
        finally:
            self.finished.emit()

