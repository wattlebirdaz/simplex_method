#include "matrix.hpp"

///////////////////////////////////////////////////////////////////////////////
//                                   Vector                                  //
///////////////////////////////////////////////////////////////////////////////

// not member /////////////////////////////////////////////////////////////////

const Vector operator+(const Vector &left, const Vector &right) {
  Vector v = left;
  return v += right;
}

const Vector operator-(const Vector &left, const Vector &right) {
  Vector v = left;
  return v -= right;
}

const Vector operator*(const Vector &left, double c) {
  Vector v = left;
  return v *= c;
}

const Vector operator*(double c, const Vector &right) {
  Vector v = right;
  return v *= c;
}

const Vector operator/(const Vector &left, double c) {
  Vector v = left;
  return v /= c;
}

// friend /////////////////////////////////////////////////////////////////////

std::ostream &operator<<(std::ostream &output, const Vector &a) {
  output.setf(std::ios_base::scientific);
  for (int i = 0; i < a.Dim; ++i) {
    output << std::setw(15) << a.ptr[i];
    if (!((i+1) % 20)) output << "\n";
  }
  if (a.Dim % 20) output << "\n";
  return output;
}


std::istream &operator>>(std::istream &input, Vector &a) {
  std::cout << "Input elements of "<< a.Dim << " dim" << " vector" << "\n";
  for (int i = 0; i < a.Dim; ++i) {
    input >> a.ptr[i];
  }
  return input;
}


double operator*(const Vector &left, const Vector &right) {
  if (left.Dim != right.Dim) {
    std::cout << "error: sizes do not match" << "\n";
    std::abort();
  }
  double a = 0.0;
  for (int i = 0; i < left.Dim; ++i) {
    a += left.ptr[i] * right.ptr[i];
  }
  return a;
}

bool operator==(const Vector &left, const Vector &right) {
  if (left.Dim != right.Dim) return false;
  for (int i = 0; i < left.Dim; ++i) {
    if (fabs(left.ptr[i] - right.ptr[i]) > NEARLY_ZERO) {
      return false;
    }
  }
  return true;
}

bool operator!=(const Vector &left, const Vector &right) {
  if (left == right) {
    return false;
  } else {
    return true;
  }
}

// public /////////////////////////////////////////////////////////////////////

Vector::Vector(int dim): Dim(dim) {
  new_vector();
  for (int i = 0; i < Dim; ++i) {
    ptr[i] = 0.0;
  }
}

Vector::Vector(const Vector &init) : Dim(init.Dim) {
  new_vector();
  for (int i = 0; i < Dim; ++i) {
    ptr[i] = init.ptr[i];
  }
}

Vector::Vector(const double *vec, int dim) : Dim(dim) {
  new_vector();
  for (int i = 0; i < Dim; ++i) {
    ptr[i] = vec[i];
  }
}

Vector::~Vector(){
  del_vector();
}

void Vector::setSize(int dim) {
  del_vector();
  Dim = dim;
  new_vector();
  for (int i = 0; i < Dim; ++i) {
    ptr[i] = 0.0;
  }
}

double Vector::norm() const {
  double a = 0.0;
  for (int i = 0; i < Dim; ++i) {
    a += ptr[i] * ptr[i];
  }
  return sqrt(a);
}

const Vector &Vector::normalize() {
  double a = norm();
  if (a < NEARLY_ZERO) return *this;
  for (int i = 0; i < Dim; ++i) {
    ptr[i] /= a;
  }
  return *this;
}

Vector &Vector::operator=(const Vector &right) {
  if (this != &right) {
    if (Dim != right.Dim) std::abort();
    for (int i = 0; i < Dim; ++i) {
      ptr[i] = right.ptr[i];
    }
  }
  return *this;
}

Vector &Vector::operator*=(double c) {
  for (int i = 0; i < Dim; ++i) {
    ptr[i] *= c;
  }
  return *this;
}

Vector &Vector::operator/=(double c) {
  if (abs(c) < NEARLY_ZERO) {
    std::cout << "error: divide by zero" << "\n";
    std::abort();
  }
  for (int i = 0; i < Dim; ++i) {
    ptr[i] /= c;
  }
  return *this;
}

Vector &Vector::operator+=(const Vector &right) {
  if (Dim != right.Dim) {
    std::cout << "error: sizes do not match" << "\n";
    std::abort();
  }
  for (int i = 0; i < Dim; ++i) {
    ptr[i] += right.ptr[i];
  }
  return *this;
}

Vector &Vector::operator-=(const Vector &right) {
  if (Dim != right.Dim) {
    std::cout << "error: sizes do not match" << "\n";
    std::abort();
  }
  for (int i = 0; i < Dim; ++i) {
    ptr[i] -= right.ptr[i];
  }
  return *this;
}

void Vector::cleanup() {
  int i;
  double max = 0.0;
  for (i = 0; i < Dim; ++i) {
    if (fabs(ptr[i]) > max) max = fabs(ptr[i]);
  }
  if (max > NEARLY_ZERO) {
    for (i = 0; i < Dim; ++i) {
      if (fabs(ptr[i]) / max < ZERO_TOLERANCE) {
        ptr[i] = 0.0;
      }
    }
  }
}

// private ////////////////////////////////////////////////////////////////////

void Vector::new_vector() {
  if (Dim == 0) {
    ptr = 0;
    return;
  }
  ptr = new double[Dim];
  if (ptr == 0) {
    std::cout << "error: failed to allocate memory" << "\n";
    std::abort();
  }
}

void Vector::del_vector() {
  delete [] ptr;
}

///////////////////////////////////////////////////////////////////////////////
//                                   Matrix                                  //
///////////////////////////////////////////////////////////////////////////////

// not member /////////////////////////////////////////////////////////////////

const Matrix operator+(const Matrix &left, const Matrix &right) {
  Matrix m = left;
  return m += right;
}

const Matrix operator-(const Matrix &left, const Matrix &right) {
  Matrix m = left;
  return m -= right;
}

// friend /////////////////////////////////////////////////////////////////////

std::istream& operator>>(std::istream& input, Matrix& a) {
  std::cout << "Input elements of "<< a.Row << "x" << a.Col << " matrix" << "\n";
  for (int i = 0; i < a.Row; ++i) {
    std::cout << "row:" << (i+1) << "\n";
    input >> a.ptr[i];
  }
  return input;
}

std::ostream& operator<<(std::ostream& output, const Matrix& a) {
  output.setf(std::ios_base::scientific);
  for (int i = 0; i < a.Row; ++i) {
    output << a.ptr[i];
  }
  return output;
}

const Vector operator*(const Matrix &a, const Vector &x) {
  if (a.Col != x.getSize()) {
    std::cout << "error: sizes do not match" << "\n";
    std::abort();
  }
  Vector y(a.Row);
  for (int i = 0; i < a.Row; ++i) {
    double sum = 0.0;
    for (int j = 0; j < a.Col; ++j) {
      sum += a.ptr[i].ptr[j] * x.ptr[j];
    }
    y.ptr[i] = sum;
  }
  y.cleanup();
  return y;
}

const Vector operator*(const Vector &x, const Matrix &a) {
  if (a.Row != x.getSize()) {
    std::cout << "error: sizes do not match" << "\n";
    std::abort();
  }
  Vector y(a.Col);
  for (int i = 0; i < a.Col; ++i) {
    double sum = 0.0;
    for (int j = 0; j < a.Row; ++j) {
      sum += x.ptr[j] * a.ptr[j].ptr[i];
    }
    y.ptr[i] = sum;
  }
  y.cleanup();
  return y;
}

const Matrix operator*(const Matrix &left, const Matrix &right) {
  if (left.Col != right.Row) {
    std::cout << "error: sizes do not match" << "\n";
    std::abort();
  }
  Matrix m(left.Row, right.Col);
  for (int i = 0; i < left.Row; ++i) {
    for (int j = 0; j < right.Col; ++j) {
      double sum = 0.0;
      for (int k = 0; k < left.Col; ++k) {
        sum += left.ptr[i].ptr[k] * right.ptr[k].ptr[j];
      }
      m.ptr[i].ptr[j] = sum;
    }
  }
  m.cleanup();
  return m;
}

bool operator==(const Matrix &left, const Matrix &right) {
  if (left.Row != right.Row || left.Col != right.Col) {
    return false;
  }
  for (int i = 0; i < left.Row; ++i) {
    if (left.ptr[i] != right.ptr[i]) {
      return false;
    }
  }
  return true;
}

bool operator!=(const Matrix &left, const Matrix &right) {
  if (left == right) {
    return false;
  } else {
    return true;
  }
}


// public /////////////////////////////////////////////////////////////////////

Matrix::Matrix(int row, int col): Row(row), Col(col) {
  new_matrix();
}

Matrix::Matrix(const Matrix &init): Row(init.Row), Col(init.Col) {
  new_matrix();
  for (int i = 0; i < Row; ++i) {
    ptr[i] = init.ptr[i];
  }
}

Matrix::~Matrix() {
  del_matrix();
}

void Matrix::setSize(int row, int col) {
  del_matrix();
  Row = row;
  Col = col;
  new_matrix();
}

Matrix &Matrix::operator=(const Matrix &right) {
  if (this != &right) {
    if ((Row != right.Row) || (Col != right.Col)) {
      std::abort();
    }
    for (int i = 0; i < Row; ++i) {
      ptr[i] = right.ptr[i];
    }
  }
  return *this;
}

Matrix &Matrix::operator+=(const Matrix &right) {
  if (Row != right.Row || Col != right.Col) {
    std::cout << "error: sizes do not match" << "\n";
    std::abort();
  }
  for (int i = 0; i < Row; ++i) {
    ptr[i] += right.ptr[i];
  }
  return *this;
}

Matrix &Matrix::operator-=(const Matrix &right) {
  if (Row != right.Row || Col != right.Col) {
    std::cout << "error: sizes do not match" << "\n";
    std::abort();
  }
  for (int i = 0; i < Row; ++i) {
    ptr[i] -= right.ptr[i];
  }
  return *this;
}

Matrix &Matrix::operator*=(const Matrix &right) {
  if ((Col != right.Row) || (Col != right.Col)) {
    std::cout << "error: sizes do not match" << "\n";
    std::abort();
  }

  Matrix m(Row, right.Col);
  for (int i = 0; i < Row; ++i) {
    for (int j = 0; j < right.Col; ++j) {
      double sum = 0.0;
      for (int k = 0; j < Col; ++k) {
        sum += ptr[i][k] * right.ptr[k][j];
      }
      m.ptr[i][j] = sum;
    }
  }
  m.cleanup();

  return *this = m;
}

void Matrix::cleanup() {
  int i, j;
  double max = 0.0;
  for (i = 0; i < Row; ++i) {
    for (j = 0; j < Col; ++j) {
      if (fabs(ptr[i][j]) > max) max = fabs(ptr[i][j]);
    }
  }
  if (max > NEARLY_ZERO) {
    for (i = 0; i < Row; ++i) {
      for (j = 0; j < Col; ++j) {
        if (fabs(ptr[i][j]) / max < ZERO_TOLERANCE) {
          ptr[i][j] = 0.0;
        }
      }
    }
  }
}

// private ////////////////////////////////////////////////////////////////////

void Matrix::new_matrix() {
  if (Row == 0 || Col == 0) {
    Row = 0;
    Col = 0;
    ptr = 0;
    return;
  }
  ptr = new Vector[Row];
  if (ptr == 0) {
    std::cout << "error: failed to allocate memory" << "\n";
    std::abort();
  }
  for (int i = 0; i < Row; ++i) {
    ptr[i].setSize(Col);
  }
}

void Matrix::del_matrix() {
  delete [] ptr;
}
