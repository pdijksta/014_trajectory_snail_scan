from PyQt5.QtCore import QObject, pyqtSignal

#class Thread(QThread):
#    def __init__(self, parent):
#        QThread.__init__(self, parent=parent)
#        self.abort = False
#
#    def stop(self):
#        self.abort = True
#        #self.wait()


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
            print('Run status')
        except Exception as e:
            self.error = True
            print('%s failed with following error message:\n%s' % (self.function, e))
        finally:
            print('finished.emit')
            self.finished.emit()

