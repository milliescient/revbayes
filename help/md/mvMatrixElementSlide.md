## name
mvMatrixElementSlide
## title
Matrix Element Sliding-Window Move
## description
Proposes an additive update to one element of a matrix.
## details
This move picks one element of a matrix at random, draws
a random number from a uniform distribution and adds the draw to the element,
as `mvSlide` does for a single parameter.

By default the element is drawn from the whole matrix. The `row` and `col`
arguments confine the move to one row or one column, counting from 1, and
giving both pins a single element.
## authors
Nicolas Lartillot
June Walker
## see_also
mvMatrixElementScale
mvSlide
## example
    moves = VectorMoves()
    x ~ dnDecomposedInvWishart(df=10, kappa=2.0, dim=4)
    # any element of the matrix
    moves.append( mvMatrixElementSlide(x, delta=0.1, weight=1) )
    # only the elements of column 2
    moves.append( mvMatrixElementSlide(x, delta=0.1, col=2, weight=1) )

## references
