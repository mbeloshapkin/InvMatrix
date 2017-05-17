# QR Decomposition

import matrix, math

type QRDecomposition = object
  qr: Matrix[float]
  m,n: int
  rDiag: seq[float]

## sqrt(a^2 + b^2) without under/overflow.
proc hypot(a,b: float): float =
  if abs(a) > abs(b):
    result = b/a
    result = abs(a) * sqrt(1.0 + result * result)
  elif b != 0.0:
    result = a/b
    result = abs(b) * sqrt(1.0 + result * result)
  else:
    result = 0.0

proc newQRDecomposition*(a: Matrix[float]): QRDecomposition =
  var qr = a
  let m = a.m
  let n = a.n
  result.m = m
  result.n = n
  result.rDiag = newSeq[float](n)

  for k in 0..<n:
    var nrm = 0.0
    for i in k..<m:
      nrm = hypot(nrm, qr[i,k])

    if nrm != 0.0:
      # Form k-th Householder vector.
      if qr[k,k] < 0:
        nrm = -nrm
      for i in k..<m:
        qr[i,k] = qr[i,k] / nrm
      qr[k,k] = qr[k,k] + 1.0
      # Apply transformation to remaining columns.
      for j in k + 1..<n:
        var s = 0.0
        for i in k..<m:
          s += qr[i,k] * qr[i,j]
        s = (-s) / qr[k,k]
        for i in k..<m:
          qr[i,j] = qr[i,j] + s * qr[i,k]

    result.rDiag[k] = -nrm

  result.qr = qr

## Return the Householder vectors
## Lower trapezoidal matrix whose columns define the reflections
proc h(this: QRDecomposition): Matrix[float] =
  result = newMatrix[float](this.m, this.n)

  for i in 0..< this.m:
    for j in 0..< this.n:
      if i >= j: result[i,j] = this.qr[i,j]
      else: result[i,j] = 0.0

## Return the upper triangular factor
proc r(this: QRDecomposition): Matrix[float] =
  result = newMatrix[float](this.n, this.n)
  for i in 0..<this.n:
    for j in 0..<this.n:
      if i < j: result[i,j] = this.qr[i,j]
      elif i == j: result[i,j] = this.rDiag[i]
      else: result[i,j] = 0.0

## Generate and return the (economy-sized) orthogonal factor
proc q(this: QRDecomposition): Matrix[float] =
  result = newMatrix[float](this.m, this.n)
  for k in countdown(this.n - 1, 0):
    for i in 0..<this.m: result[i,k] = 0.0
    result[k,k] = 1.0

    for j in k..<this.n:
      if this.qr[k,k] != 0.0:
        var s = 0.0
        for i in k..<this.m: s += this.qr[i,k] * result[i,j]
        s = -s / this.qr[k,k]
        for i in k..<this.m: result[i,j] = result[i,j] + s * this.qr[i,k]

## Is the matrix full rank
## true if R, and hence A, has full rank.
proc isFullRank(this: QRDecomposition): bool =
  for j in 0..<this.n:
    if this.rDiag[j] == 0.0: return false
  result = true

## Least squares solution of A*X = B
## B is a Matrix with as many rows as A and any number of columns.
## X that minimizes the two norm of Q*R*X-B
proc solve*(this: QRDecomposition, b: Matrix): Matrix =
  doAssert(b.m == this.m, "Matrix row dimensions must agree.")
  doAssert(this.isFullRank, "Matrix is rank deficient.")

  let nx = b.n
  var x = b

  # Compute Y = transpose(Q)*B
  for k in 0..<this.n:
    for j in 0..<nx:
      var s = 0.0
      for i in k..<this.m: s += this.qr[i,k] * x[i,j]
      s = -s / this.qr[k,k]
      for i in k..<this.m: x[i,j] = x[i,j] + s * this.qr[i,k]

  # Solve R*X = Y
  for k in countdown(this.n - 1, 0):
    for j in 0..<nx: x[k,j] = x[k,j] / this.rDiag[k]
    for i in 0..<k:
      for j in 0..<nx: x[i,j] = x[i,j] - x[k,j] * this.qr[i,k]

  result = x.getMatrix(0, this.n - 1, 0, nx - 1)

     
