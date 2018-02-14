import numpy as np
import random
import threading
import time
from multiprocessing import Process, Queue, freeze_support, set_start_method
import matplotlib.pyplot as plt
import tkinter as tk
import tkinter.ttk as tt
from tkinter import *

def matrix_thread(A, B, C, K, N, n1, n2):
    for i in range(n1, n2):
        for j in range(0, K):
            for k in range(0, N):
                C[i, j] += A[i, k] * B[k, j]


def matrix_multy(q, A, B, K, N, n1, n2):
    C_new = 0
    list_of_C = []
    for i in range(n1, n2):
        for j in range(0, K):
            for l in range(0, N):
                C_new += A[i, l] * B[l, j]
            list_of_C.append([C_new, i, j])

    q.put(list_of_C)


if __name__ == '__main__':

    n=1
    seq =[]
    thred =[]
    mult =[]
    xstep = [100,200,300,400]
    for N in xstep:

        A = np.zeros((N, N))
        B = np.zeros((N, N))
        C = np.zeros((N, N))

        for i in range(0, N):
            for j in range(0, N):
                A[i][j] = random.randrange(-10, 10)

        for i in range(0, N):
            for j in range(0, N):
                B[i][j] = random.randrange(-10, 10)
        seq1 = 0
        thred1=0
        mult1=0
        for iter in range(0,n):
            time_start = time.time()
            matrix_thread(A, B, C, N, N, 0, N)
            seq1+=(time.time() - time_start)
        seq.append(seq1/n)

        for iter in range(0, n):
            list_of_threading = []
            list_of_Process = []
            list_of_Queue = []

            freeze_support()
            if (N==10 & n==0):
                set_start_method('spawn')


            for i in range(0, 4):
                list_of_threading.append(threading.Thread(target=matrix_thread, args=(A, B, C, N, N, int(i * N / 4), int(i * N / 4 + N / 4))))
                q = Queue()
                list_of_Queue.append(q)
                list_of_Process.append(Process(target=matrix_multy, args=(q, A, B, N, N, int(i * N / 4), int(i * N / 4 + N / 4))))



            time_start = time.time()
            for elem in list_of_threading:
                elem.start()
            for elem in list_of_threading:
                elem.join()
            thred1+=time.time() - time_start


            time_start = time.time()
            for elem in list_of_Process:
                elem.start()

            for elem in list_of_Queue:
                data = elem.get()
                for tt in data:
                    C[tt[1], tt[2]] = tt[0]
            mult1+=time.time() - time_start
        thred.append(thred1 / n)
        mult.append(mult1 / n)
    print(seq)
    print(thred)
    print(mult)
    plt.figure(1)
    plt.plot(xstep, seq)
    plt.plot(xstep,mult)
    plt.plot(xstep, thred)
    plt.title(r'$seq,\ mult,\ thred$')
    plt.show()
    for i in range(0,len(xstep)):
        mult[i] = mult[i]/seq[i]
        thred[i] = thred[i]/seq[i]
    plt.figure(2)
    plt.plot(xstep, seq)
    plt.plot(xstep,mult)
    plt.plot(xstep, thred)
    plt.title(r'$seq,\ mult,\ thred$')
    plt.show()