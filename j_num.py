from Get_global_value import BB
from Get_global_value import SE


def j_num(num_e):
    n = len(SE)
    ie = []  # list tye
    connection = []

    for i in range(n):
        if SE[i] == 1:
            ie.append(i)

    j_number = BB[ie[num_e]]

    connection.append(ie[num_e])  # connection是一个list
    while j_number != -1:
        connection.insert(0, j_number)  # 在首部插入 j_number
        j_number = BB[j_number]

    joint = connection   # return a list

    return joint



























