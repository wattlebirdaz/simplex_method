#include "matrix.hpp"
#include <random>
#include <chrono>
#include <cmath>

using std::cout;
using std::endl;
using std::random_device;
using std::mt19937;
using std::uniform_real_distribution;
using std::normal_distribution;
using std::tuple;
using std::make_tuple;
using std::get;
using std::abort;
using std::max;
using std::min;
using std::swap;
using std::pow;

// Transpose
Matrix T(Matrix A) {
  Matrix AT(A.getCol(), A.getRow());
  for (int i = 0; i < A.getRow(); i++) {
    for (int j = 0; j < A.getCol(); j++) {
      AT[j][i] = A[i][j];
    }
  }
  return AT;
}

Vector sign(Vector V) {
  Vector signs(V.getSize());
  for (int i = 0; i < V.getSize(); i++) {
    signs[i] = (V[i] < 0) ? -1 : (V[i] == 0) ? 0 : 1;
  }
  return signs;
}

Vector rand(int n) {
  Vector V(n);
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> d(0, 1);
  for (int i = 0; i < n; i++) {
    V[i] = d(gen);
  }
  return V;
}

Vector randn(int n) {
  Vector V(n);
  random_device rd; 
  mt19937 gen(rd()); 
  normal_distribution<double> d(0, 1);
  for (int i = 0; i < n; i++) {
    V[i] = d(gen);
  }
  return V;
}

Matrix randn(int n, int m) {
  Matrix A(n, m);
  random_device rd; 
  mt19937 gen(rd()); 
  normal_distribution<double> d(0, 1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      A[i][j] = d(gen);
    }
  }
  return A;
}

tuple<Matrix, Vector, Vector> create_feasible_bounded_problem(int m, int n, bool nonnegative_b = false) {
  Vector x = rand(n);
  Vector y = rand(m);
  Matrix A = randn(m, n);
  if (nonnegative_b) {
    Vector signs = sign(A * x);
    for (int i = 0; i < signs.getSize(); i++) {
      A[i] = signs[i] * A[i];
    }
  } else {
    while (true) {
      Vector signs = sign(A * x);
      int allsigns = 1;
      for (int i = 0; i < signs.getSize(); i++) {
        if (signs[i] < 0) {
          allsigns = -1;
          break;
        }
      }
      if (allsigns == -1) {
        break;
      } else {
        A = randn(m, n);
      }
    }
  }
  Vector b = A * x;
  Vector c = T(A) * y;
  return make_tuple(A, b, c);
}

tuple<Matrix, Vector, Vector> create_infeasible_problem(int m, int n) {
  double eps = 1e-15;
  Vector y = rand(m);
  Vector b(m);
  while (true) {
    b = randn(m);
    if (y * b <= -eps) break;
  }
  Matrix A = randn(n, m);
  Vector signs = sign(A * y);
  for (int i = 0; i < signs.getSize(); i++) {
    A[i] = signs[i] * A[i];
  }
  Vector c = randn(n);
  return make_tuple(T(A), b, c);
}

tuple<Matrix, Vector, Vector> create_unbounded_problem(int m, int n) {
  tuple<Matrix, Vector, Vector> tupleValue = create_infeasible_problem(n, m);
  Matrix A = get<0>(tupleValue);
  Vector b = get<1>(tupleValue);
  Vector c = get<2>(tupleValue);
  for (int i = 0; i < A.getRow(); i++) {
    A[i] = -1 * A[i];
  }
  return make_tuple(T(A), -c, -b);
}

template<class T> inline bool chmax(T& a, T b) { 
  if (a < b) {
    a = b;
    return true;
  }
  return false;
}

template<class T> inline bool chmin(T& a, T b) {
  if (a > b) {
    a = b;
    return true;
  }
  return false;
}

#define INF 1e5;

Matrix simplexMethod(Matrix &All) {
  // 誤差消去
  All.cleanup();
  // cout << All << endl;
  int m = All.getRow() - 1;
  int n = All.getCol() - 1 - m;
  // 目的関数の係数最大(実際は符号が反転しているので最小)の変数を見つける
  double S = INF;
  int SI = -1;
  
  // 最大係数規則
  for (int j = 0; j < m+n; j++) {
    SI = (chmin(S, All[m][j]) ? j : SI);
  }

  // ブランドの選択規則
  // for (int j = 0; j < m+n; j++) {
  //   if (All[m][j] < 0) {
  //     SI = (chmin(S, All[m][j]) ? j : SI);
  //     break;
  //   }
  // }
  
  // 最小の係数が0以上なら終了, 0より小さかったら続ける
  if (S < 0) {
    if (SI == -1) abort();

    // 変数の最小の可能増分を見つける
    double P = INF;
    int PI = -1;
    for (int i = 0; i < m; i++) {
      if (All[i][SI] > 0) {
        PI = (chmin(P, All[i][n+m] / All[i][SI]) ? i : PI);
      }
    }

    if (PI == -1) {
      cout << "Main problem" << endl;
      cout << All;
      throw "unbounded";
    };
    
    // 以下の行, 列をピボットにして掃き出しを行う
    int ROW = PI;
    int COL = SI;

    All[ROW] /= All[ROW][COL];
  
    for (int i = 0; i < m+1; i++) {
      if (i == ROW) continue;
      double d = All[i][COL];
      for (int j = 0; j < m+n+1; j++) {
        All[i][j] -= d * All[ROW][j];
      }
    } 
    simplexMethod(All);
  }
  
  return All;
}

tuple<Matrix, Vector> simplexMethod2(Matrix &All, Vector &Obj) {
  // 誤差消去
  All.cleanup();
  // cout << All << endl;
  int m = All.getRow() - 1;
  int n = All.getCol() - 1 - m;
  // 目的関数の係数最大(実際は符号が反転しているので最小)の変数を見つける
  double S = INF;
  int SI = -1;
  
  // 最大係数規則
  for (int j = 0; j < m+n; j++) {
    SI = (chmin(S, All[m][j]) ? j : SI);
  }

  // ブランドの選択規則
  // for (int j = 0; j < m+n; j++) {
  //   if (All[m][j] < 0) {
  //     SI = (chmin(S, All[m][j]) ? j : SI);
  //     break;
  //   }
  // }
  
  // 最小の係数が0以上なら終了, 0より小さかったら続ける
  if (S < 0) {
    if (SI == -1) abort();

    // 変数の最小の可能増分を見つける
    double P = INF;
    int PI = -1;
    for (int i = 0; i < m; i++) {
      if (All[i][SI] > 0) {
        PI = (chmin(P, All[i][n+m] / All[i][SI]) ? i : PI);
      }
    }

    if (PI == -1) {
      cout << "Main problem" << endl;
      cout << All;
      throw "unbounded";
    };
    
    // 以下の行, 列をピボットにして掃き出しを行う
    int ROW = PI;
    int COL = SI;

    All[ROW] /= All[ROW][COL];
  
    for (int i = 0; i < m+1; i++) {
      if (i == ROW) continue;
      double d = All[i][COL];
      for (int j = 0; j < m+n+1; j++) {
        All[i][j] -= d * All[ROW][j];
      }
    }
    Obj -= Obj[COL] * All[ROW];
    simplexMethod2(All, Obj);
  }
  
  return make_tuple(All, Obj);
}

Matrix subProblem(const Matrix &A, const Vector &b, const Vector &c) {
  int m = b.getSize();
  int n = c.getSize();
  int cnt = 0;
  for (int i = 0; i < m; i++) {
    if (b[i] < 0) cnt++;
  }
  Matrix Sub(m+1, n+m+cnt+1);
  for (int i = 0; i < Sub.getRow(); i++) {
    for (int j = 0; j < Sub.getCol(); j++) {
      if (i < m && j < n) Sub[i][j] = A[i][j];
      else if (i < m && j >= n && j < n+m) Sub[i][j] = (j-n == i) ? 1 : 0;
      else if (i < m && j == n+m+cnt) Sub[i][j] = b[i];
      else if (i == m && j < n+m) Sub[i][j] = 0;
      else if (i == m && j >= n+m  && j < n+m+cnt) Sub[i][j] = 1;
      else if (i == m && j == n+m+cnt) Sub[i][j] = 0;
    }
  }

  int flag = 0;
  for (int i = 0; i < m; i++) {
    if (b[i] < 0) {
      Sub[i][n+m+flag] = -1;
      flag++;
      Sub[i] = -Sub[i];
      Sub[m] -= Sub[i];
    }
  }

  Vector C(Sub.getCol());
  for (int i = 0; i < Sub.getCol(); i++) {
    if (i < n) C[i] = -c[i];
    else C[i] = 0;
  }

  simplexMethod2(Sub, C);

  if (Sub[m][n+m+cnt] < 0) {
    cout << "Sub problem" << endl;
    cout << Sub;
    throw "infeasible";
  };
  
  Matrix Main(m+1, n+m+1);
  // 一番右の行
  for (int i = 0; i < m+1; i++) {
    if (i != m) Main[i][n+m] = Sub[i][n+m+cnt];
    else if (i == m) Main[i][n+m] = C[n+m+cnt];
  }
  // それ以外
  for (int i = 0; i < m+1; i++) {
    for (int j = 0; j < n+m; j++) {
      if (i != m) Main[i][j] = Sub[i][j];
      else if (i == m) Main[i][j] = C[j];
    }
  }

  return Main;
}

// 入力行列を作る(bに負の要素がある場合はそれを考慮して人工変数を入れる)
Matrix createMatrix(const Matrix &A, const Vector &b, const Vector &c) {
  int m = b.getSize();
  int n = c.getSize();
  Matrix All(m+1, n+m+1);
  for (int i = 0; i < m+1; i++) {
    for (int j = 0; j < n+m+1; j++) {
      if (i < m && j < n) All[i][j] = A[i][j];
      else if (i < m && j >= n && j < n+m) All[i][j] = (j-n == i) ? 1 : 0;
      else if (i < m && j == n+m) All[i][j] = b[i];
      else if (i == m && j < n) All[i][j] = -c[j];
      else if (i == m && j >= n) All[i][j] = 0;
    }
  }
  return All;
}

// 結果出力
void showResult(const Matrix &All) {
  int m = All.getRow() - 1;
  int n = All.getCol() - 1 - m;
  Vector X(n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      if (All[i][j]  == 1.0) {
        X[j] = All[i][n+m];
      }
    }
  }
  cout << "X: " << X;
  cout << "f: " << All[m][n+m] << endl;
}

int main() {
  // (i)-1
  cout << "(i)-1" << endl;
  for (int i = 1; i <= 5; i++) {
    double T = 0.0;
    for (int j = 0; j < 10; j++) {
      tuple<Matrix, Vector, Vector> tupleValue = create_feasible_bounded_problem(10, i*10, true);
      Matrix A = get<0>(tupleValue);
      Vector b = get<1>(tupleValue);
      Vector c = get<2>(tupleValue);
      auto start = std::chrono::high_resolution_clock::now();
      Matrix All = createMatrix(A, b, c);
      simplexMethod(All);
      auto end = std::chrono::high_resolution_clock::now();
      double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
      T += elapsed;
    }
    cout << "m: 10 n: " << i*10 << " time: " << T/10 << " [μs]" << endl;
  }


  // (i)-2
  cout << "(i)-2" << endl;
  for (int i = 1; i <= 5; i++) {
    double T = 0.0;
    for (int j = 0; j < 10; j++) {
      tuple<Matrix, Vector, Vector> tupleValue = create_feasible_bounded_problem(i*10, 10, true);
      Matrix A = get<0>(tupleValue);
      Vector b = get<1>(tupleValue);
      Vector c = get<2>(tupleValue);
      auto start = std::chrono::high_resolution_clock::now();
      Matrix All = createMatrix(A, b, c);
      simplexMethod(All);
      auto end = std::chrono::high_resolution_clock::now();
      double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
      T += elapsed;
    }
    cout << "m: " << i*10 << " n: 10" << " time: " << T/10 << " [μs]" << endl;
  }


  // (ii)-1
  cout << "(ii)-1" << endl;
  for (int i = 1; i <= 7; i++) {
    tuple<Matrix, Vector, Vector> tupleValue = create_feasible_bounded_problem(10, pow(10, i), true);
    Matrix A = get<0>(tupleValue);
    Vector b = get<1>(tupleValue);
    Vector c = get<2>(tupleValue);
    auto start = std::chrono::high_resolution_clock::now();
    Matrix All = createMatrix(A, b, c);
    simplexMethod(All);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    cout << "m: 10" << " n: 10^" << i << " time: " << elapsed << " [μs]" << endl;
  }

  // (ii)-2
  cout << "(ii)-2" << endl;
  for (int i = 1; i <= 4; i++) {
    tuple<Matrix, Vector, Vector> tupleValue = create_feasible_bounded_problem(pow(10, i), 10, true);
    Matrix A = get<0>(tupleValue);
    Vector b = get<1>(tupleValue);
    Vector c = get<2>(tupleValue);
    auto start = std::chrono::high_resolution_clock::now();
    Matrix All = createMatrix(A, b, c);
    simplexMethod(All);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    cout << "m: 10^" << i << " n: 10" << " time: " << elapsed << " [μs]" << endl;
  }

  // (iii)-1
  cout << "(iii)-1" << endl;
  for (int i = 1; i <= 5; i++) {
    double T = 0.0;
    for (int j = 0; j < 10; j++) {
      tuple<Matrix, Vector, Vector> tupleValue = create_feasible_bounded_problem(10, i*10, false);
      Matrix A = get<0>(tupleValue);
      Vector b = get<1>(tupleValue);
      Vector c = get<2>(tupleValue);
      auto start = std::chrono::high_resolution_clock::now();
      Matrix All = subProblem(A, b, c);
      simplexMethod(All);
      auto end = std::chrono::high_resolution_clock::now();
      double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
      T += elapsed;
    }
    cout << "m: 10 n: " << i*10 << " time: " << T/10 << " [μs]" << endl;
  }


  // (iii)-2
  cout << "(iii)-2" << endl;
  for (int i = 1; i <= 5; i++) {
    double T = 0.0;
    for (int j = 0; j < 10; j++) {
      tuple<Matrix, Vector, Vector> tupleValue = create_feasible_bounded_problem(i*10, 10, false);
      Matrix A = get<0>(tupleValue);
      Vector b = get<1>(tupleValue);
      Vector c = get<2>(tupleValue);
      auto start = std::chrono::high_resolution_clock::now();
      Matrix All = subProblem(A, b, c);
      simplexMethod(All);
      auto end = std::chrono::high_resolution_clock::now();
      double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
      T += elapsed;
    }
    cout << "m: " << i*10 << " n: 10" << " time: " << T/10 << " [μs]" << endl;
  }


  cout << "(iv)-unbounded" << endl;
  {
    tuple<Matrix, Vector, Vector> tupleValue = create_unbounded_problem(2, 3);
    Matrix A = get<0>(tupleValue);
    Vector b = get<1>(tupleValue);
    Vector c = get<2>(tupleValue);
    cout << "A" << endl;
    cout << A;
    cout << "b" << endl;
    cout << b;
    cout << "c" << endl;
    cout << c;
    Matrix All = subProblem(A, b, c);
    try {
      simplexMethod(All);
    } catch (char const* str) {
      cout << str << endl;
    }
  }

  cout << "(iv)-infeasible" << endl;
  {
    tuple<Matrix, Vector, Vector> tupleValue = create_infeasible_problem(2, 3);
    Matrix A = get<0>(tupleValue);
    Vector b = get<1>(tupleValue);
    Vector c = get<2>(tupleValue);
    cout << "A" << endl;
    cout << A;
    cout << "b" << endl;
    cout << b;
    cout << "c" << endl;
    cout << c;
    try {
      Matrix All = subProblem(A, b, c);
    } catch (char const* str) {
      cout << str << endl;
    }
  }
}
