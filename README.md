This Simplex LP solver is implemented according to the CLRS textbook.

The names of the variables also comes from the textbook. 

The main difference is, matrix A is with the same sign as in the slack form, instead of the negatives as in the textbook. Same sign is more readable and easier to follow, though negatives might be somewhat more efficient.

All of A, b and c are implemented with dictionaries instead of arrays.

v is wrapped in c[-1] so that it can be modified in place.

In the slack form, A.keys() is equivalent to B, and A[i].keys() is equivalent to N, so both variables are omitted. Additionally, b.keys() is equivalent to B, c.keys() is equivalent to N.

```
Usage:
    python simplex.py filename [-no]
```
Use -no flag to disable show-steps.

The input should be in a separate file with the following format.

The first line contains two number: number of variables n and number of constraints m

The input should be in standard form, which implies a maximization problem.

The next m lines contains the matrix A, with n elements in each row.

The next line contains the vector b, and last line contains c.

For example, for 
```
maximize  	x1 + 2 x2 + 0.5 x3
subject to  x1 + x2 + x3 <= 5
            x1 <= 3
            x2 <= 1
            x3 <= 4
            x1, x2, x3 >= 0
```
The input should be
```
3 4
1 1 1
1 0 0
0 1 0
0 0 1
5 3 1 4
1 2 0.5
```
Note that variable index starts from 0 instead of 1 in the algorithm, so they will actually become x0, x1, x2
