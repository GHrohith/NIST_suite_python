from math import floor, sqrt, erfc
from scipy.special import gammaincc
from numpy import matrix
from numpy.linalg import matrix_rank

def frequency(stream: str):
    n = len(stream)
    s_n = 0

    for bit in stream:
        s_n += 2 * int(bit) - 1

    s_obs = abs(s_n) / sqrt(n)

    p_value = erfc(s_obs / sqrt(2))
    return p_value

def blockFrequency(stream: str):
    n = len(stream)
    blockSize = max(20, 0.02*n)
    blockNums = floor(n / blockSize)

    X_obs_temp = 0

    for i in range(blockNums):
        pi_i  = stream[(i*blockSize):((i+1)*blockSize)].count("1") / blockSize
        print(pi_i)

        X_obs_temp += (pi_i - 0.5) ** 2
    
    X_obs = 4 * blockSize * X_obs_temp

    p_value = gammaincc(blockNums/2, X_obs/2)

    return p_value

def matrixRank(stream: str):
    n = len(stream)
    rowNum = 3
    colNum = 3
    blockNums = floor(n / (rowNum*colNum))

    Fm = Fm_1 = 0

    for i in range(blockNums):
        mat_arr = []

        for j in range(rowNum):
            mat_arr.append([int(bit) for bit in stream[(i*rowNum*colNum + j*colNum):(i*rowNum*colNum + (j+1)*colNum)]])

        mat = matrix(mat_arr)

        rank = matrix_rank(mat)

        if rank == rowNum:
            Fm += 1
        elif rank == rowNum - 1:
            Fm_1 += 1
    
    X_obs = ((Fm - 0.2888*blockNums)**2)/(0.2888*blockNums) + \
            ((Fm_1 - 0.5776*blockNums)**2)/(0.5776*blockNums) + \
            ((blockNums - Fm - Fm_1 - 0.1336*blockNums)**2)/(0.1336*blockNums)
    
    p_value = gammaincc(2/2, X_obs/2)
    return(p_value)

def main():
    str1 = "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000"
    str2 = "0110011010"
    str3 = "01011001001010101101"

    # print(frequency(str))
    # print(blockFrequency(str1))
    print(matrixRank(str3))

if __name__ == "__main__":
    main()