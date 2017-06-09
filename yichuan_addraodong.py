# -*- coding:UTF-8 -*-
#用三个基因片段表示三个变量，把三个基因片段当成一组解，每个变量用16位表示

#自变量编码为二进制，计算目标函数，计算适应度函数，自然选择，基因交换，突变
import math
import random
import matplotlib.pyplot as plt

# 定义全局变量
m_size = 1000 # 种群大小
gen_length = 16 # 基因长度
terminal_g = 150 # 迭代次数
p_gc = 0.90 #交叉概率
p_sc = 0.01 # 突变概率
# 种群初始化
theta1_m = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for i in range(m_size)]
theta2_m = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for i in range(m_size)]
theta3_m = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for i in range(m_size)]


# 二进制转十进制
def binary2decimal(b):
    value = 0
    for i in range(len(b)):
        value += b[i]* math.pow(2,i)
    value = value * 360.0/math.pow(2,gen_length)
    value -= 180
    # 从(0,360)转化(-180,180)
    return value;


# 计算目标函数值
def function_value(theta1,theta2,theta3):
    c1 = math.cos(theta1/180.0*math.pi)
    c2 = math.cos(theta2/180.0*math.pi)
    c3 = math.cos(theta3/180.0*math.pi)

    s1 = math.sin(theta1 / 180.0 * math.pi)
    s2 = math.sin(theta2 / 180.0 * math.pi)
    s3 = math.sin(theta3 / 180.0 * math.pi)

    p1 = 330*c1*s2-55*(c1*c2*s3-c1*c3*s2)
    p2 = 330*s1*s2-55*(c2*s1*s3-c3*s1*s2)
    p3 = 330*c2+55*c2*c3+55*s2*s3+275

    value = math.pow((-166.7099-p1),2) + math.pow((-288.75-p2),2) + math.pow((412.5-p3),2)
    g = math.pow(value,1.0/2)+ 0.1*(abs(theta1/180.0*math.pi) + \
        abs(theta2/180.0*math.pi) + abs(theta3/180.0*math.pi))

    return g

# 计算适应度函数值，最小化结果，则都转化成用最大值减去该值作为适应度函数
def calculate_p():
    result = {} # 存目标函数值
    p = [0] * m_size # 存轮盘赌占有的概率
    p_random = [0] * m_size # 存随机生成的概率值
    best_value = [10000.0] * 4 #存储最好的结果
    # 计算目标函数值,所有值都是正的
    for i in range(m_size):
        theta1 = binary2decimal(theta1_m[i])
        theta2 = binary2decimal(theta2_m[i])
        theta3 = binary2decimal(theta3_m[i])
        result[i] = function_value(theta1,theta2,theta3)
        # 更新最好结果
        if(result[i] < best_value[0]):
            best_value[0] = result[i]
            best_value[1] = theta1
            best_value[2] = theta2
            best_value[3] = theta3

        p_random[i] = random.random()

    # 2*(i-1)/(m_size-1)
    temp = sorted(result.iteritems(), key=lambda d:d[1], reverse = True)
    for i in range(len(temp)):
        result[temp[i][0]] = 2*(i + 1 - 1)/(m_size-1.0)
    # 求自然选择的概率
    p[0] = result[0]
    for i in range(1,m_size):
        p[i] = p[i - 1] + result[i]
    if (not (p[-1] == 0)):
        for i in range(1,m_size):
            p[i] = p[i] / p[-1]
    p_random.sort()
    return p,p_random,best_value


# 自然选择产生新种群
def nature_selection(p,p_random):
    new_index = 0
    old_index = 0
    new_theta1_m = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for i in range(m_size)]
    new_theta2_m = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for i in range(m_size)]
    new_theta3_m = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for i in range(m_size)]
    while (new_index < m_size):
        if (p_random[new_index] <= p[old_index]):
            new_theta1_m[new_index] = theta1_m[old_index]
            new_theta2_m[new_index] = theta2_m[old_index]
            new_theta3_m[new_index] = theta3_m[old_index]
            new_index += 1
        else:
            old_index += 1
    return new_theta1_m,new_theta2_m,new_theta3_m

# 个体交叉
def gen_cross():
    for i in range(m_size-1):
        # 三个变量是否需要交叉都是独立的
        p = random.random()
        if(p < p_gc):
            # 随机选择交叉点,和下一个基因交叉
            point = random.randint(0,gen_length)
            temp = theta1_m[i][point:]
            theta1_m[i][point:] = theta1_m[i+1][point:]
            theta1_m[i+1][point:] = temp
        p = random.random()
        if (p < p_gc):
            # 第二个变量是否需要交叉,和下一个基因交叉
            point = random.randint(0, gen_length)
            temp = theta2_m[i][point:]
            theta2_m[i][point:] = theta2_m[i + 1][point:]
            theta2_m[i + 1][point:] = temp
        p = random.random()
        if (p < p_gc):
            # 第二个变量是否需要交叉,和下一个基因交叉
            point = random.randint(0, gen_length)
            temp = theta3_m[i][point:]
            theta3_m[i][point:] = theta3_m[i + 1][point:]
            theta3_m[i + 1][point:] = temp

# 基因突变
def gen_suddenchange():
    for i in range(m_size):
        # 三个变量之间是否突变不相互影响
        p = random.random()
        if(p < 1):
            point = random.randint(0,gen_length-1)
            theta1_m[i][point] = (theta1_m[i][point] + 1) % 2
        p = random.random()
        if (p < 1):
            point = random.randint(0, gen_length - 1)
            theta2_m[i][point] = (theta2_m[i][point] + 1) % 2
        p = random.random()
        if (p < 1):
            point = random.randint(0, gen_length - 1)
            theta3_m[i][point] = (theta3_m[i][point] + 1) % 2


# 画图并输出结果
def output(best_value):
    x = [i + 1 for i in range(terminal_g)]
    y = [best_value[i][0] for i in range(terminal_g)]
    plt.plot(x, y)  # 画连线图
    plt.scatter(x, y)  # 画散点图
    plt.show()

    best_value.sort()
    print "最终结果如下(g,theta1,theta2.theta3)："
    print best_value[0]


# 十进制转二进制
def decimal2binary(value):
    value += 180
    value = int(value * math.pow(2,gen_length)/360)
    temp = []
    while(value/2 != 0):
        temp = [value%2] + temp
        value = value/2
    temp = [value%2] + temp
    temp = [0]*(gen_length-len(temp))+temp
    return temp


# 增加扰动
def raodong():
    for i in range(m_size):
        theta1 = binary2decimal(theta1_m[i])
        theta2 = binary2decimal(theta2_m[i])
        theta3 = binary2decimal(theta3_m[i])
        new_theta1 = (theta1 + random.uniform(-20,20))/400.0*360.0
        new_theta2 = (theta2 + random.uniform(-20, 20))/400.0*360.0
        new_theta3 = (theta3 + random.uniform(-20, 20))/400.0*360.0

        # 结果变好，接受扰动
        if(function_value(theta1,theta2,theta3) < function_value(new_theta1,new_theta2,new_theta3)):
            theta1_m[i] = decimal2binary(new_theta1)
            theta2_m[i] = decimal2binary(new_theta2)
            theta3_m[i] = decimal2binary(new_theta3)
        else:
            p = random.random()
            if(p < 0.1):
                theta1_m[i] = decimal2binary(new_theta1)
                theta2_m[i] = decimal2binary(new_theta2)
                theta3_m[i] = decimal2binary(new_theta3)


# 程序入口函数
def run():
    best_value = []
    for i in range(terminal_g):
        p,p_random,value = calculate_p()
        best_value.append(value) # 存最好的结果
        theta1_m,theta2_m,theta3_m = nature_selection(p,p_random) # 自然选择
        gen_cross() # 交叉
        gen_suddenchange() # 基因突变
        raodong() #添加

    output(best_value) # 输出结果

run()
