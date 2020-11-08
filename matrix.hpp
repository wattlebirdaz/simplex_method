#pragma once
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#define NDEBUG
#include <cassert>

#define NEARLY_ZERO 1.E-10
#define ZERO_TOLERANCE 1.E-12

class Vector;
class Matrix;

const Vector operator+(const Vector &, const Vector &);
const Vector operator-(const Vector &, const Vector &);
const Vector operator*(double, const Vector &);
const Vector operator*(const Vector &, double);
const Vector operator/(const Vector &, double);

class Vector {
  friend std::ostream &operator<<(std::ostream &, const Vector &);
  friend std::istream &operator>>(std::istream &, Vector &);
  friend double operator*(const Vector &, const Vector &);
  friend const Matrix operator*(const Matrix &, const Matrix &);
  friend const Vector operator*(const Matrix &, const Vector &);
  friend const Vector operator*(const Vector &, const Matrix &);
  friend bool operator==(const Vector &, const Vector &);
  friend bool operator!=(const Vector &, const Vector &);
  
 public:
  explicit Vector(int = 0);
  Vector(const Vector &);
  Vector(const double *, int);
  ~Vector();
  void setSize(int);
  int getSize() const { return Dim; }
  double norm() const;
  const Vector &normalize();
  Vector &operator=(const Vector &);
  double &operator[](int);
  const double &operator[](int) const;
  const Vector operator-() const {return -1.0*(*this);}
  Vector &operator*=(double);
  Vector &operator/=(double);
  Vector &operator+=(const Vector &);
  Vector &operator-=(const Vector &);
  void cleanup();

 private:
  double *ptr;
  int Dim;
  void new_vector();
  void del_vector();
};


const Matrix operator+(const Matrix &, const Matrix &);
const Matrix operator-(const Matrix &, const Matrix &);

class Matrix {
  friend std::ostream &operator<<(std::ostream &, const Matrix &);
  friend std::istream &operator>>(std::istream &, Matrix &);
  friend const Vector operator*(const Matrix &, const Vector &);
  friend const Vector operator*(const Vector &, const Matrix &);
  friend const Matrix operator*(const Matrix &, const Matrix &);
  friend bool operator==(const Matrix &, const Matrix &);
  friend bool operator!=(const Matrix &, const Matrix &);
  
 public:
  explicit Matrix(int = 0, int = 0);     /* default constructor */
  Matrix(const Matrix &);       /* copy constructor */
  ~Matrix();                    /* destructor */
  void setSize(int, int);
  int getRow() const {return Row;}
  int getCol() const {return Col;}
  Matrix &operator=(const Matrix &);
  Vector &operator[](int);
  const Vector &operator[](int) const;
  Matrix &operator+=(const Matrix &);
  Matrix &operator-=(const Matrix &);
  Matrix &operator*=(const Matrix &);
  void cleanup();

 private:
  Vector *ptr;
  int Row;
  int Col;

  /* utility functions */
  void new_matrix();            /* allocate memory */
  void del_matrix();            /* release memory */
};

inline double &Vector::operator[](int i) {
  assert(i >= 0 && i < Dim);
  return ptr[i];
}

inline const double &Vector::operator[](int i) const {
  assert(i >= 0 && i < Dim);
  return ptr[i];
}

inline Vector &Matrix::operator[](int i) {
  assert(i >= 0 && i < Row);
  return ptr[i];
}

inline const Vector &Matrix::operator[](int i) const {
  assert(i >= 0 && i < Row);
  return ptr[i];
}

