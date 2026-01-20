#!/usr/bin/python

import queue
import threading
import time
import logging

exitFlag = 0

class myThread (threading.Thread):
    def __init__(self, threadID, queueLock, q, func):

        threading.Thread.__init__(self)
        self.threadID = threadID
        self.q = q
        self.func = staticmethod(func)

    def run(self):
        self.process_data()
    
    def process_data(self, *args):
        while not exitFlag:
            self.queueLock.acquire()
            if not self.q.empty():
                data = self.q.get()
                self.queueLock.release()
                self.func()
            else:
                self.queueLock.release()
            time.sleep(1)


def drucken(a, b):
    print(a + " und " + b)

queueLock = threading.Lock()
workQueue = queue.Queue(10)
threads = []
threadID = 1

# Create new threads
for threadID in range(10):
    thread = myThread(threadID, queueLock, workQueue, drucken, )
    thread.start()
    threads.append(thread)
    threadID += 1

# Fill the queue
queueLock.acquire()
for word in nameList:
   workQueue.put(word)
queueLock.release()

# Wait for queue to empty
while not workQueue.empty():
   pass

# Notify threads it's time to exit
exitFlag = 1

# Wait for all threads to complete
for t in threads:
   t.join()
