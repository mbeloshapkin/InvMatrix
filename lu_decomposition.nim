# Lower Upper decomposition by Alan Turing (1948)
# slightly modified by me for stable Kriging solution.
# https://en.wikipedia.org/wiki/LU_decomposition or google LU Decomposition
#
import matrix, math

type LUDecomposition* = object
  lu*: Matrix[float]
  piv*: seq[int]
  m, n, pivsign: int

proc newLUDecomposition*(a: Matrix[float]): LUDecomposition =
  result.lu = a
  result.m = a.m
  result.n = a.n
  result.piv = newSeq[int](a.m)
  #just make the code more readable
  let m = a.m
  let n = a.n

  for i in 0..<a.m:
    result.piv[i] = i

  result.pivsign = 1

  var luRowi = newVector[float](n)
  var luColj = newVector[float](m)

  for j in 0..<n:
  # Make a copy of the j-th column to localize references.
    for i in 0..<m:
      luColj[i] = result.lu[i,j]
  # Apply previous transformations.
    for i in 0..<m:
      luRowi = result.lu.getRow(i)
      let kmax = min(i,j)
      var s = 0.0
      for k in 0..<kmax:
        s = s + luRowi[k] * luColj[k]

      luColj[i] = luColj[i] - s
      luRowi[j] = luColj[i]
      result.lu.setRow(i, luRowi)

    ## Find pivot and exchange if necessary.
    var p = j
    for i in (j + 1)..<m:
      if abs(luColj[i]) > abs(luColj[p]):
        p = i

    if p != j:
      for k in 0..<n:
        let t = result.lu[p,k]
        result.lu[p,k] = result.lu[j,k]
        result.lu[j,k] = t

      let k2 = result.piv[p]
      result.piv[p] = result.piv[j]
      result.piv[j] = k2
      result.pivsign = -result.pivsign

    if result.lu[j,j] == 0.0:
      result.lu[j,j] = 0.0000000000001

    if j < m and result.lu[j,j] != 0.0:
      for ii in j + 1..<m:
        result.lu[ii,j] = result.lu[ii,j] / result.lu[j,j];

## Return lower triangular factor
proc getL(this: LUDecomposition): Matrix[float] =
  result = newMatrix[float](this.m, this.n)
  for i in 0..<this.m:
    for j in 0..<this.n:
      if i > j:
        result[i,j] = this.lu[i,j]
      else:
        if i == j: result[i,j] = 1.0
        else: result[i,j] = 0.0

proc getU(this: LUDecomposition): Matrix[float] =
  result = newMatrix[float](this.n, this.n)
  for i in 0..<this.n:
    for j in 0..<this.n:
      if i <= j: result[i,j] = this.lu[i,j]
      else: result[i,j] = 0.0

proc isNonSingular(this: LUDecomposition): bool =
  for j in 0..<this.n:
    if this.lu[j,j] == 0:
      result = false
      return
  result = true

##Solve A*X = B</summary>
##  B - A Matrix with as many rows as A and any number of columns.
## Returns X so that L*U*X = B(piv,:)
## exception Matrix row dimensions must agree.
## exception  Matrix is singular.
proc solve*(this: LUDecomposition, B: Matrix): Matrix =
  doAssert(B.m == this.m, "Matrix row dimensions must agree.")
  doAssert(this.isNonSingular, "Matrix is singular.")

  ## Copy right hand side with pivoting
  let nx = B.n

  var xMat = B.getMatrix(this.piv, 0, nx - 1);

  ## Solve L*Y = B(piv,:)
  for k in 0..<this.n:
    for i in k + 1..<this.n:
      for j in 0..< nx:
        xMat[i,j] = xMat[i,j] - xMat[k,j] * this.lu[i,k]

  ## Solve U*X = Y
  for k in countdown(this.n - 1,0): #(int k = n - 1; k >= 0; k--)
    for j in 0..<nx:
      xMat[k,j] = xMat[k,j] / this.lu[k,k]
    for i in 0..<k:
      for j in 0..<nx:
        xMat[i,j] = xMat[i,j] - xMat[k,j] * this.lu[i,k]

  result = xMat











