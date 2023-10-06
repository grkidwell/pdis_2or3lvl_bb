import os,sys
'''
PROJECT_ROOT = os.path.abspath(os.path.join(
                  os.path.dirname("__file__"), 
                  os.pardir)
)
'''
libdir=os.getcwd()
sys.path.append(libdir)
sys.path.append(libdir+r'\libs')
#print("ran code")