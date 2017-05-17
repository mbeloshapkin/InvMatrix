# Test
import unittest, matrix, identity, lu_decomposition, qr_decomposition

proc printVector(v: Vector) =
  for x in v: echo x

proc printMatrix(m: Matrix) =
  for j in 0..<m.m:
    var r = m.getRow(j)
    echo r


suite "Vector and Matrix operations":
  test "Vector initialisation":
    let v = newVector(@[1.0, 3.0, 2.0, 8.0, -2.0])
    printVector(v)
  test "Matrix initialisation":
    let
      m = newMatrix[float](4, 3, @[
        1.0, 0.0, 2.0, -1.0,
        -1.0, 1.0, 3.0, 1.0,
        3.0, 2.0, 2.0, 4.0
      ])
    echo m[1,1]
  test "Get row and column":
    let
      m = newMatrix[float](4,3, @[
        1.0, 0.0, 2.0, -1.0,
        -1.0, 1.0, 3.0, 1.0,
        3.0, 2.0, 2.0, 4.0
      ])
    let v = getRow(m,2)
    echo "3d row"
    printVector(v)
    let v1 = m.getColumn(1)
    echo "2nd column"
    printVector(v1)

  test "LU Decomposition 2 by 2":
    var m = newMatrix[float](2,2)

    m[0,0] = 3.0
    m[0,1] = 1.0
    m[1,0] = 2.0
    m[1,1] = 4.0

    let b = newMatrix(2,1,@[1.0, 3.0])

    echo "M ="
    printMatrix(m)
    echo "B ="
    printMatrix(b)

    let luDec = newLUDecomposition(m)
    #echo "LU ="
    #printMatrix(luDec.lu)

    let x = luDec.solve(b)

    echo "X ="
    printMatrix(x)

    let mx = m * x
    echo "M * X ="
    printMatrix(mx)

  var m = newMatrix[float](3,3)

  m[0,0] = 1.0
  m[0,1] = 2.0
  m[0,2] = 3.0
  m[1,0] = 3.0
  m[1,1] = 5.0
  m[1,2] = 7.0
  m[2,0] = 1.0
  m[2,1] = 3.0
  m[2,2] = 4.0

  test "LU Decomposition 3 by 3":

    let b = newMatrix(3,1,@[3.0, 0.0, 1.0])

    echo "M ="
    printMatrix(m)
    echo "B ="
    printMatrix(b)

    let luDec = newLUDecomposition(m)
    #echo "LU ="
    #printMatrix(luDec.lu)

    let x = luDec.solve(b)

    echo "X ="
    printMatrix(x)

    let mx = m * x
    echo "M * X ="
    printMatrix(mx)

  test "Inverse matrix by LU Decomposition":
    echo "The matrix ="
    printMatrix(m)

    let luDec = newLUDecomposition(m)

    let id = newIdentity[float](m.m, m.m)  # matrix unity
    let mInv = luDec.solve(id)      # m * mInv = 1 equation

    echo "Inversed matrix ="
    printMatrix(mInv)

    echo "Product M * M_Inv ="
    let id1 = m * mInv
    printMatrix(id1) # That shall be identity

  test "Inverse matrix by QR Decomposition":
    echo "The matrix ="
    printMatrix(m)

    let qrDec = newQRDecomposition(m)

    let id = newIdentity[float](m.m, m.m)  # matrix unity
    let mInv = qrDec.solve(id)      # m * mInv = 1 equation

    echo "Inversed matrix ="
    printMatrix(mInv)

    echo "Product M * M_Inv ="
    let id1 = m * mInv
    printMatrix(id1) # That shall be identity
echo "Done"

