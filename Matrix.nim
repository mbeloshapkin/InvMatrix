# General Matrix

type
   Matrix*[T] = object
     width, height: int
     data: seq[T]

   Vector*[T] = seq[T]

method m*(this: Matrix): int = this.height # Number of rows

method n*(this: Matrix): int = this.width # Number of columns

proc newMatrix*[T](height, width: int): Matrix[T] =
  result.width = width
  result.height = height
  result.data = newSeq[T](width*height)

proc newMatrix*[T](height, width: int, data: seq[T]): Matrix[T] =
  result.width = width
  result.height = height
  result.data = data

proc newVector*[T](src: seq[T]): Vector[T] = src

proc newVector*[T](size: int): Vector[T] = newSeq[T](size)

proc `[]`*[T](mat: Matrix[T], rowIdx, colIdx: int): T =
  mat.data[rowIdx * mat.width + colIdx]

proc `[]`*[T](mat: var Matrix[T], rowIdx, colIdx: int): T =
  mat.data[rowIdx * mat.width + colIdx]


proc `[]=`*[T](mat: var Matrix[T], rowIdx, colIdx: int, val: T) =
   mat.data[rowIdx * mat.width + colIdx] = val

proc `*`*[T](this: Matrix[T], b: Matrix[T]): Matrix[T] =
  doAssert(b.m == this.n, "Matrix inner dimensions must agree.")
  result = newMatrix[T](this.m, b.n)

  var bColj = newVector[T](this.n)
  for j in 0..<b.n:
    for k in 0..<this.n:
      bColj[k] = b[k,j]

    for i in 0..<this.m:
      let aRowi = this.getRow(i)
      var s = 0.0
      for k in 0..<this.n:
        s += aRowi[k] * bColj[k]
      result[i,j] = s

proc getRow*[T](mat: Matrix[T], rowIdx: int): Vector[T] =
  result = newVector[T](mat.width)
  for ct in 0..<mat.width:
    result[ct] = mat.data[rowIdx * mat.width + ct]

proc setRow*[T](mat: var Matrix[T], rowIdx: int, val: Vector[T]) =
  for ct in 0..<mat.width:
    mat[rowIdx, ct] = val[ct]

proc getColumn*[T](mat: Matrix[T], colIdx: int): Vector[T] =
  result = newVector[T](mat.height)
  for ct in 0..<mat.height:
    result[ct] = mat.data[colIdx + ct*mat.width]

## Get a submatrix
## r - Array of row indices
## j0 - Initial column index
## j1 -  Final column index
proc getMatrix*[T](mat: Matrix[T], r: openArray[int], j0, j1: int): Matrix[T] = #Matrix[j1-j0+1, N, T] =   #array[M * (j1-j0+1), T] =
  result = newMatrix[T](r.len, j1-j0+1)
  for i in 0..<r.len:
    for j in  j0..j1:
      result[i,j-j0] = mat[r[i],j]

## Get a submatrix.
## i0 - Initial row index
## i1 - Final row index
## j0 - Initial column index
## j1 - Final column index
proc getMatrix*[T](mat: Matrix[T], i0, i1, j0, j1: int): Matrix[T] =
  result = newMatrix[T](i1 - i0 + 1, j1 - j0 + 1)
  for i in i0..i1:
    for j in j0..j1:
      result[i - i0, j - j0] = mat[i, j]
