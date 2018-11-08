'''this module checks whether a specified directory exists and if it doesn't it creates the directory'''
import os

class isdir:
    @staticmethod
    def isdir(f):
        ''' Checks for directory and if needed makes it. Enter full path i.e. "C:\Users\owner" '''
        if not os.path.isdir(f):
            os.mkdir(f)

    
if __name__ == "__main__" :
    isdir.isdir("C:\Users\User\Desktop\isdir test\ ")
    