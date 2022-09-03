#from time import sleep
from threading import *


#WITHOUT extending the Thread class

#import threading
class ex:
    def B(self):
        lst = [2,1,'w',8.7,'abc']
        for x in lst:
            print("Child printing ",x)

myobj = ex()
t1=Thread(target=myobj.B)
t1.start()
t1.join()
print('done')
'''
#THREAD WITHOUT CLASS
def new1():
    for x in range(6):
        print('Child Executing...',current_thread().getName())
        sleep(1)

def new2():
    for x in range(6):
        print('Baby Executing...',current_thread().getName())
        sleep(1)

t1=Thread(target=new1)
print(current_thread().getName())
t2=Thread(target=new2)
t1.start()
t2.start()
t1.join()
t2.join()
print('Bye',current_thread().getName())






class Hello(Thread):
    def run(self):
        for i in range(5):
            print('Hello')
            sleep(1)
class Hi(Thread):
    def run(self):
        for i in range(5):
            print('Hi')
            sleep(1)

t1 = Hello()
t2 = Hi()

t1.start()

t2.start()
t1.join()
t2.join()
print('Bye')
'''
