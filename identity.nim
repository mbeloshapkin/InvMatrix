# Identity matrix generator
# It just contains all zeroes except main diagonal which contains unities

import matrix

proc newIdentity*[T](height, width: int): Matrix[T] =
  result = newMatrix[T](height, width)
  for i in 0..<height:
    for j in 0..<width:
      if i == j: result[i,j] = 1.0
      else: result[i,j] = 0.0      
