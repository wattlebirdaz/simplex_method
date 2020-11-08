# Simplex method (Linear Programming) in C++

maximize c^Tx s.t. Ax \leq b and x \geq 0

## Usage

### solve feasible bounded problem where all elements in b is \geq 0


```c++
    tuple<Matrix, Vector, Vector> tupleValue = create_feasible_bounded_problem(2, 3, true); // A is 2*3 matrix
    Matrix A = get<0>(tupleValue);
    Vector b = get<1>(tupleValue);
    Vector c = get<2>(tupleValue);
    Matrix All = createMatrix(A, b, c);
    simplexMethod(All);
    showResult(All);
```

### solve feasible bounded problem where an element < 0 exists in b

```c++
    tuple<Matrix, Vector, Vector> tupleValue = create_feasible_bounded_problem(2, 3, false);
    Matrix A = get<0>(tupleValue);
    Vector b = get<1>(tupleValue);
    Vector c = get<2>(tupleValue);
    Matrix All = createMatrix(A, b, c);
    simplexMethod(All);

```


### solve unbounded problem

```c++
    tuple<Matrix, Vector, Vector> tupleValue = create_unbounded_problem(2, 3);
    Matrix A = get<0>(tupleValue);
    Vector b = get<1>(tupleValue);
    Vector c = get<2>(tupleValue);
    Matrix All = subProblem(A, b, c);
    try {
      simplexMethod(All); // throws error when unbounded
    } catch (char const* str) {
      cout << str << endl; // prints "unbounded"
    }
```

### solve infeasible problem

```c++
    tuple<Matrix, Vector, Vector> tupleValue = create_infeasible_problem(2, 3);
    Matrix A = get<0>(tupleValue);
    Vector b = get<1>(tupleValue);
    Vector c = get<2>(tupleValue);
    try {
      Matrix All = subProblem(A, b, c); // throws error when infeasible
    } catch (char const* str) {
      cout << str << endl; // prints "infeasible"
    }
```


## Time measurements

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
<colgroup>
<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<tbody>
<tr>
<td class="org-left">Processor Name</td>
<td class="org-left">Dual-Core Intel Core i5</td>
</tr>
<tr>
<td class="org-left">Processor Speed</td>
<td class="org-left">2 GHz</td>
</tr>
<tr>
<td class="org-left">Number of Processors</td>
<td class="org-left">1</td>
</tr>
<tr>
<td class="org-left">Total Number of Cores</td>
<td class="org-left">2</td>
</tr>
<tr>
<td class="org-left">L2 Cache (per Core)</td>
<td class="org-left">256 KB</td>
</tr>
<tr>
<td class="org-left">L3 Cache</td>
<td class="org-left">4 MB</td>
</tr>
<tr>
<td class="org-left">Memory</td>
<td class="org-left">8 GB</td>
</tr>
</tbody>
</table>


```shell
g++ simplex_method.cpp matrix.cpp
```

```
| => ./a.out
(i)-1
m: 10 n: 10 time: 156.3 [μs]
m: 10 n: 20 time: 179.1 [μs]
m: 10 n: 30 time: 222.9 [μs]
m: 10 n: 40 time: 233.3 [μs]
m: 10 n: 50 time: 368.9 [μs]
(i)-2
m: 10 n: 10 time: 137.8 [μs]
m: 20 n: 10 time: 444.8 [μs]
m: 30 n: 10 time: 1081.2 [μs]
m: 40 n: 10 time: 2076.1 [μs]
m: 50 n: 10 time: 2674.7 [μs]
(ii)-1
m: 10 n: 10^1 time: 134 [μs]
m: 10 n: 10^2 time: 610 [μs]
m: 10 n: 10^3 time: 4398 [μs]
m: 10 n: 10^4 time: 43958 [μs]
m: 10 n: 10^5 time: 370515 [μs]
m: 10 n: 10^6 time: 4.15787e+06 [μs]
m: 10 n: 10^7 time: 8.44057e+07 [μs]
(ii)-2
m: 10^1 n: 10 time: 637 [μs]
m: 10^2 n: 10 time: 8281 [μs]
m: 10^3 n: 10 time: 658616 [μs]
m: 10^4 n: 10 time: 1.38351e+08 [μs]
(iii)-1
m: 10 n: 10 time: 729.8 [μs]
m: 10 n: 20 time: 342.8 [μs]
m: 10 n: 30 time: 326.2 [μs]
m: 10 n: 40 time: 331.4 [μs]
m: 10 n: 50 time: 373.4 [μs]
(iii)-2
m: 10 n: 10 time: 213.7 [μs]
m: 20 n: 10 time: 1324.7 [μs]
m: 30 n: 10 time: 3755.7 [μs]
m: 40 n: 10 time: 6028.6 [μs]
m: 50 n: 10 time: 10970.5 [μs]
(iv)-unbounded
A
   9.856944e-01  -3.056640e-02  -2.177825e+00
   2.867853e-01  -1.153169e-01  -1.008329e+00
b
  -9.176064e-01  -1.296517e+00
c
   1.411685e+00  -2.035293e+00   3.045551e+00
Main problem
  -2.844165e-01   1.143644e-01   1.000000e+00   0.000000e+00  -9.917402e-01   1.285808e+00
   3.662850e-01   2.184994e-01   0.000000e+00   1.000000e+00  -2.159837e+00   1.882658e+00
  -2.277890e+00   2.383595e+00   0.000000e+00   0.000000e+00  -3.020395e+00   3.915993e+00
unbounded
(iv)-infeasible
A
  -4.695574e-01   1.687160e+00  -4.138367e-01
   3.711978e-01   1.932663e-01   1.779096e+00
b
  -2.069227e+00  -1.018689e+00
c
  -5.211084e-02  -5.421527e-01  -4.694724e-01
Sub problem
   1.000000e+00  -3.593085e+00   8.813335e-01  -2.129665e+00   0.000000e+00   2.129665e+00   0.000000e+00   4.406760e+00
   0.000000e+00  -1.527012e+00  -1.451947e+00  -7.905269e-01  -1.000000e+00   7.905269e-01   1.000000e+00   2.654468e+00
   0.000000e+00   1.527012e+00   1.451947e+00   7.905269e-01   1.000000e+00   2.094731e-01   0.000000e+00  -2.654468e+00
infeasible
```
