# Debug utils

import matrix

proc echoMatrix*(m: Matrix) =
  for j in 0..<m.m:
    var r = m.getRow(j)
    echo r

proc echoVector*(title: string, v: Vector) = echo title, " = ", v
