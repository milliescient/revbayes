## name
mvMatrixElementScale
## title
Matrix Element Scaling Move
## description
Proposes a multiplicative update to one element of a matrix.
## details
This move picks one element of a matrix at random and
multiplies it with a randomly chosen factor, as `mvScale` does for a single
parameter.

By default the element is drawn from the whole matrix. The `row` and `col`
arguments confine the move to one row or one column, counting from 1, and
giving both pins a single element.
## authors
Nicolas Lartillot
June Walker
## see_also
mvMatrixElementSlide
mvMatrixElementSwap
mvScale
## example
    moves = VectorMoves()
    x ~ dnDecomposedInvWishart(df=10, kappa=2.0, dim=4)
    # any element of the matrix
    moves.append( mvMatrixElementScale(x, lambda=0.1, weight=1) )
    # only the elements of column 2
    moves.append( mvMatrixElementScale(x, lambda=0.1, col=2, weight=1) )

## references
