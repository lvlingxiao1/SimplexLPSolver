import numpy
import sys
from math import inf

'''
This Simplex LP solver is implemented according to the CLRS textbook.
The names of the variables also comes from the textbook. 
The main difference is, matrix A is with the same sign as in the slack form, instead of the negatives as
in the textbook. Same sign is more readable and easier to follow, though negatives might be somewhat more efficient.
All of A, b and c are implemented with dictionaries instead of arrays.
v is wrapped in c[-1] so that it can be modified in place.
In the slack form, A.keys() is equivalent to B, and A[i].keys() is equivalent to N, so both variables are
omitted. Additionally, b.keys() is equivalent to B, c.keys() is equivalent to N.
'''

INSTRUCTION = '''
Use -no flag to disable show-steps.
The input should be in a separate file with the following format.
The first line contains two number: number of variables n and number of constraints m
The input should be in standard form, which implies a maximization problem.
The next m lines contains the matrix A, with n elements in each row.
The next line contains the vector b, and last line contains c.

For example, for maximize  x1 + 2 x2 + 0.5 x3
               subject to  x1 + x2 + x3 <= 5
                           x1 <= 3
                           x2 <= 1
                           x3 <= 4
                           x1, x2, x3 >= 0
The input should be
3 4
1 1 1
1 0 0
0 1 0
0 0 1
5 3 1 4
1 2 0.5
Note that index starts from 0 instead of 1 in the algorithm, so they are actually x0, x1, x2
'''

SHOW_STEPS = True


def parse_input(s):
    s = s.split('\n')
    row = s[0].split()
    try:
        num_var = int(row[0])
        num_con = int(row[1])

        A = []
        for i in range(num_con):
            row = s[i + 1].split()
            A.append([])
            for j in range(num_var):
                A[i].append(float(row[j]))

        b = []
        row = s[num_con + 1].split()
        for i in range(num_con):
            b.append(float(row[i]))

        c = []
        row = s[num_con + 2].split()
        for i in range(num_var):
            c.append(float(row[i]))
        return A, b, c, num_var

    except ValueError:
        print('Input error.' + INSTRUCTION)
        sys.exit(1)


def to_slack_form(std_A, std_b, std_c, num_var):
    A = {}
    b = {}
    for index, row in enumerate(std_A):
        i = num_var + index
        A[i] = dict()
        for j, val in enumerate(row):
            A[i][j] = -val
        b[i] = std_b[index]

    c = {}
    for i in range(num_var):
        c[i] = std_c[i]

    c[-1] = 0  # c[-1]=v in the textbook, need a pointer for modifying in place

    return A, b, c


def print_slack(A, b, c):
    l = []
    for i in c:
        if i != -1:
            l.append(f'{c[i]} x{i}')
    print(f'z = {c[-1]} + {" + ".join(l)}')
    for i in A:
        l = []
        for j in A[i]:
            if A[i][j] != 0:
                l.append(f'{A[i][j]} x{j}')
        print(f'x{i} = {b[i]} + {" + ".join(l)}')


def can_improve(c):
    for i in c:
        if c[i] > 0 and i != -1:
            return i
    return -1


def simplex(A, b, c, num_var):
    A, b, c = to_slack_form(A, b, c, num_var)
    print_slack(A, b, c)
    # if any(i < 0 for i in B):
    #     get_first_feasible_solution(A, b, c)
    e = can_improve(c)
    while e != -1:
        tightest_bound = inf
        l = None
        for i in A:
            if A[i][e] < 0:
                bound = - b[i] / A[i][e]
                if bound < tightest_bound:
                    tightest_bound = bound
                    l = i
        if tightest_bound == inf:
            print('This LP is unbounded.')
            sys.exit(0)
        pivot(A, b, c, e, l)

        if SHOW_STEPS:
            print_slack(A, b, c)

        e = can_improve(c)

    print("\nSolution:")
    for i in range(num_var):
        if i in b:
            print(f'x{i}: {b[i]}')
        else:
            print(f'x{i}: 0')


def pivot(A, b, c, e, l):
    print(f'==>Pivot on x{e} with row of x{l}')
    b[e] = - b[l] / A[l][e]
    del b[l]

    A[e] = {}
    for i in A[l]:
        if i != e:
            A[e][i] = -A[l][i]/A[l][e]
    A[e][l] = 1/A[l][e]
    del A[l]

    for i in A:
        if i != e:
            b[i] += A[i][e]*b[e]
            for j in A[i]:
                if j != e:
                    A[i][j] += A[i][e]*A[e][j]
            A[i][l] = A[i][e] * A[e][l]
            del A[i][e]

    c[-1] += c[e]*b[e]
    for i in c:
        if i != e and i != -1:
            c[i] += c[e]*A[e][i]
    c[l] = c[e]*A[e][l]
    del c[e]


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage:\n\tpython simplex.py filename [-no]')
        # print(INSTRUCTION)
        EXAMPLE = '''3 4    # 3 variables and 4 constraints
                1 1 1
                1 0 0
                0 1 0
                0 0 1
                5 3 1 4     # this line contains vector b
                1 2 0.5     # this line contains vector c
                '''
        s = EXAMPLE
    else:
        f = open(sys.argv[1])
        s = f.read()
        if len(sys.argv) >= 3:
            if sys.argv[2].startswith('-no'):
                SHOW_STEPS = False
    A, b, c, num_var= parse_input(s)
    simplex(A, b, c, num_var)
